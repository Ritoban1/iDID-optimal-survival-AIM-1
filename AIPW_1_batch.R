# pass arguments to file
cat("Running on node:", system("hostname", intern = TRUE), "\n")
sessionInfo()
args = commandArgs(trailingOnly = T)

# takes first argument that was passed
seed = as.numeric(args[1])
# set seed (later name file according to seed)
set.seed(seed)
library(bridgedist)
library(rgenoud)
library(simsurv)
library(dplyr)
library(survival)
library(SuperLearner)
library(survSuperLearner)
n=1e3
t=1
tau=10
incr=0.01
setwd("~/Survival_DID/d_tilde_fixed/")
source("SL_functions.R")
source("Data_generator.R")

### d1_V
ref=readRDS("result_IPW_large_1.rds")
model_d1<-ref$model_d1
pred_d1=predict(model_d1,newdata=cbind.data.frame(1,V1,V2),type="response")
d1_V=as.numeric(ifelse(pred_d1>=0.5,1,0))

pi_LZ<- function(L_train, Z_train, V_train, L_test, Z_test, V_test){
  
  V_train_df <- as.data.frame(V_train)
  V_test_df   <- as.data.frame(V_test)
  
  #    Bind L to the covariates, then regress L against everything else (.)
  data_L  <- cbind(L_target = L_train, V_train_df)
  model_L <- glm(L_target ~ ., data = data_L, family = "binomial")
  
  # 3. Model Z (Using Training Data)
  data_Z  <- cbind(Z_target = Z_train, V_train_df)
  model_Z <- glm(Z_target ~ ., data = data_Z, family = "binomial")
  
  # 4. Predict Probabilities (Using New/Test Data)
  l_p <- predict(model_L, newdata = V_test_df, type = "response") # P(L=1 | V)
  z_p <- predict(model_Z, newdata = V_test_df, type = "response") # P(Z=1 | V)
  
  # 5. Calculate likelihood/joint probability on New Data
  #    (Fixed the bug where you referenced 'Z' instead of 'Z_test')
  l_out <- (Z_test * L_test * l_p * z_p) +                 # Case: L=1, Z=1
    ((1 - Z_test) * L_test * (1 - z_p) * l_p) +     # Case: L=1, Z=0
    (Z_test * (1 - L_test) * z_p * (1 - l_p)) +     # Case: L=0, Z=1
    ((1 - Z_test) * (1 - L_test) * (1 - z_p) * (1 - l_p)) # Case: L=0, Z=0
  
  return(l_out)
}

delta_A <- function(A_train, L_train, Z_train,
                    L_test, Z_test,
                    V_train, V_test) {
  
  # Ensure V inputs are data frames (handles matrices gracefully)
  # This ensures column names are consistent for predict()
  V_train_df <- as.data.frame(V_train)
  V_test_df  <- as.data.frame(V_test)
  
  # Helper: Fit on specific Train subset, Predict on Test set
  fit_predict_subset <- function(cond_L, cond_Z) {
    
    # Identify training rows
    subset_mask <- (L_train == cond_L & Z_train == cond_Z)
    
    df_train <- cbind(A = A_train[subset_mask], 
                      V_train_df[subset_mask, , drop = FALSE])
    
    # Fit model using "A ~ ." (A against Everything else in df_train)
    model <- glm(A ~ ., data = df_train, family = "binomial")
    
    # Predict on the Test Data
    return(predict(model, newdata = V_test_df, type = "response"))
  }
  
  # Calculate the 4 components
  p11 <- fit_predict_subset(cond_L = 1, cond_Z = 1)
  p10 <- fit_predict_subset(cond_L = 1, cond_Z = 0)
  p01 <- fit_predict_subset(cond_L = 0, cond_Z = 1)
  p00 <- fit_predict_subset(cond_L = 0, cond_Z = 0)
  
  # Calculate delta
  delta <- as.numeric(p11 - p10 - p01 + p00)
  p_orig<- (Z_test * L_test * p11) + 
    ((1 - Z_test) * L_test * p10) + 
    (Z_test * (1 - L_test) * p01) + 
    ((1 - Z_test) * (1 - L_test) * p00)
  
  return(list(delta_A=delta,mu_A=p_orig))
}
delta_N_tilde <- function(obs_N_train, delta_train, V_train, L_train, Z_train, 
                          V_test, L_test, Z_test, k.x.m_train,
                          sl_lib = c("SL.glm", "SL.lgb")) {
  
  # 1. Prepare Training Data ---------------------------------------------------
  make_df <- function(V, L, Z) {
    d <- as.data.frame(V)
    d$L <- L
    d$Z <- Z
    return(d)
  }
  
  # Y variable construction
  Y_train <- delta_train * obs_N_train / k.x.m_train
  
  # Full Training Data (for Delta model)
  X_train_full <- make_df(V_train, L_train, Z_train)
  
  # Subset Training Data (for Y model, only where delta/observed == 1)
  obs_idx <- which(delta_train == 1)
  X_train_obs <- X_train_full[obs_idx, , drop = FALSE]
  Y_train_obs <- Y_train[obs_idx]
  
  
  # 2. Fit Models (DO THIS ONCE ONLY) ------------------------------------------
  
  # Fit Y Model: E[Y | V, L, Z, Delta=1]
  sl_fit_Y <- SuperLearner(
    Y = Y_train_obs,
    X = X_train_obs,
    family = gaussian(),
    SL.library = sl_lib,
    method = "method.NNLS"
  )
  
  # Fit Delta Model: P(Delta=1 | V, L, Z)
  sl_fit_D <- SuperLearner(
    Y = delta_train,
    X = X_train_full,
    family = binomial(),
    SL.library = sl_lib,
    method = "method.CC_LS"
  )
  
  
  # 3. Create Counterfactual Test Data -----------------------------------------
  get_counterfactual_data <- function(v_data, l_val, z_val) {
    d <- as.data.frame(v_data)
    d$L <- l_val
    d$Z <- z_val
    return(d)
  }
  
  X_test_11 <- get_counterfactual_data(V_test, 1, 1)
  X_test_10 <- get_counterfactual_data(V_test, 1, 0)
  X_test_01 <- get_counterfactual_data(V_test, 0, 1)
  X_test_00 <- get_counterfactual_data(V_test, 0, 0)
  
  
  # 4. Predict on Test Data ----------------------------------------------------
  get_preds <- function(sl_model, newdata) {
    p <- predict(sl_model, newdata = newdata, onlySL = TRUE)$pred
    return(as.numeric(p))
  }
  
  # Predict Y component
  py_11 <- get_preds(sl_fit_Y, X_test_11)
  py_10 <- get_preds(sl_fit_Y, X_test_10)
  py_01 <- get_preds(sl_fit_Y, X_test_01)
  py_00 <- get_preds(sl_fit_Y, X_test_00)
  
  # Predict Delta component
  pd_11 <- get_preds(sl_fit_D, X_test_11)
  pd_10 <- get_preds(sl_fit_D, X_test_10)
  pd_01 <- get_preds(sl_fit_D, X_test_01)
  pd_00 <- get_preds(sl_fit_D, X_test_00)
  
  
  # 5. Compute mu_N and Delta_N_tilde ------------------------------------------
  
  # Calculate mu_N for each counterfactual scenario
  mu_N11 <- py_11 * pd_11
  mu_N10 <- py_10 * pd_10
  mu_N01 <- py_01 * pd_01
  mu_N00 <- py_00 * pd_00
  
  # Final Contrast Calculation
  delta_N_tilde <- mu_N11 - mu_N10 - mu_N01 + mu_N00
  
  # 6. Compute mu_N for the OBSERVED test data ---------------------------------
  # Pick the correct mu_N based on the actual L_test and Z_test values
  mu_N_orig <- (Z_test * L_test * mu_N11) + 
    ((1 - Z_test) * L_test * mu_N10) + 
    (Z_test * (1 - L_test) * mu_N01) + 
    ((1 - Z_test) * (1 - L_test) * mu_N00)
  
  return(list(delta_N_tilde = delta_N_tilde, mu_N = as.numeric(mu_N_orig)))
}
delta_D <- function(i_x_t_train, delta_train, V_train, L_train, Z_train, 
                    V_test, L_test, Z_test, k.x.m_train,
                    sl_lib = c("SL.glm", "SL.lgb")) {
  
  # 1. Prepare Training Data ---------------------------------------------------
  make_df <- function(V, L, Z) {
    d <- as.data.frame(V)
    d$L <- L
    d$Z <- Z
    return(d)
  }
  
  # Y variable construction
  Y_train <- delta_train * i_x_t_train / k.x.m_train
  
  # Full Training Data (for Delta model)
  X_train_full <- make_df(V_train, L_train, Z_train)
  
  # Subset Training Data (for Y model, only where delta/observed == 1)
  obs_idx <- which(delta_train == 1)
  X_train_obs <- X_train_full[obs_idx, , drop = FALSE]
  Y_train_obs <- Y_train[obs_idx]
  
  
  # 2. Fit Models (DO THIS ONCE ONLY) ------------------------------------------
  
  # Fit Y Model: E[Y | V, L, Z, Delta=1]
  sl_fit_Y <- SuperLearner(
    Y = Y_train_obs,
    X = X_train_obs,
    family = gaussian(),
    SL.library = sl_lib,
    method = "method.NNLS"
  )
  
  # Fit Delta Model: P(Delta=1 | V, L, Z)
  sl_fit_D <- SuperLearner(
    Y = delta_train,
    X = X_train_full,
    family = binomial(),
    SL.library = sl_lib,
    method = "method.CC_LS"
  )
  
  
  # 3. Create Counterfactual Test Data -----------------------------------------
  get_counterfactual_data <- function(v_data, l_val, z_val) {
    d <- as.data.frame(v_data)
    d$L <- l_val
    d$Z <- z_val
    return(d)
  }
  
  X_test_11 <- get_counterfactual_data(V_test, 1, 1)
  X_test_10 <- get_counterfactual_data(V_test, 1, 0)
  X_test_01 <- get_counterfactual_data(V_test, 0, 1)
  X_test_00 <- get_counterfactual_data(V_test, 0, 0)
  
  
  # 4. Predict on Test Data ----------------------------------------------------
  get_preds <- function(sl_model, newdata) {
    p <- predict(sl_model, newdata = newdata, onlySL = TRUE)$pred
    return(as.numeric(p))
  }
  
  # Predict Y component
  py_11 <- get_preds(sl_fit_Y, X_test_11)
  py_10 <- get_preds(sl_fit_Y, X_test_10)
  py_01 <- get_preds(sl_fit_Y, X_test_01)
  py_00 <- get_preds(sl_fit_Y, X_test_00)
  
  # Predict Delta component
  pd_11 <- get_preds(sl_fit_D, X_test_11)
  pd_10 <- get_preds(sl_fit_D, X_test_10)
  pd_01 <- get_preds(sl_fit_D, X_test_01)
  pd_00 <- get_preds(sl_fit_D, X_test_00)
  
  
  # 5. Compute mu_N and Delta_N_tilde ------------------------------------------
  
  # Calculate mu_N for each counterfactual scenario
  mu_D11 <- py_11 * pd_11
  mu_D10 <- py_10 * pd_10
  mu_D01 <- py_01 * pd_01
  mu_D00 <- py_00 * pd_00
  
  # Final Contrast Calculation
  delta_D <- mu_D11 - mu_D10 - mu_D01 + mu_D00
  
  # 6. Compute mu_N for the OBSERVED test data ---------------------------------
  # Pick the correct mu_N based on the actual L_test and Z_test values
  mu_D_orig <- (Z_test * L_test * mu_D11) + 
    ((1 - Z_test) * L_test * mu_D10) + 
    (Z_test * (1 - L_test) * mu_D01) + 
    ((1 - Z_test) * (1 - L_test) * mu_D00)
  
  return(list(delta_D = delta_D, mu_D = as.numeric(mu_D_orig)))
}
surv_fit<-function(X_train,delta_train,
                   A_train,A_test,L_train,L_test,
                   Z_train,Z_test,V_train,V_test,
                   cens_library,event_library,
                   eval.times,eps = 1e-06){
  
  cov_train=cbind.data.frame(A=A_train,L=L_train,Z=Z_train,V=V_train)
  cov_test=cbind.data.frame(A=A_test,L=L_test,Z=Z_test,V=V_test)
  
  fit <- survSuperLearner(
    time = X_train,
    event = delta_train,
    X = cov_train,
    newX = cov_test,
    new.times = eval.times,
    event.SL.library = cens_library,
    cens.SL.library = event_library
  )
  res=NULL
  res$h.hat<-fit$event.SL.predict
  if(any(res$h.hat == 0)) res$h.hat[res$h.hat == 0] <- min(res$h.hat[res$h.hat > 0]) # prevent 0
  
  res$k.hat <- fit$cens.SL.predict
  if(any(res$k.hat == 0)) res$k.hat[res$k.hat == 0] <- min(res$k.hat[res$k.hat > 0]) # prevent 0
  pred <- predict.survSuperLearner(fit,
                                   newdata = cov_train,
                                   new.times = eval.times)
  k.fitted <- pred$cens.SL.predict
  if(any(k.fitted == 0)) k.fitted[k.fitted == 0] <- min(k.fitted[k.fitted > 0]) # prevent 0
  res$k.x.m <- sapply(1:nrow(cov_train),
                      function(i) stepfun(eval.times, c(1,
                                                        k.fitted[i,]), right = FALSE)(X_train[i]-eps))
  return(res)
}

b_fit<-function(obs_N_train,X_train,delta_train,
                A_train,A_test,L_train,L_test,
                Z_train,Z_test,V_train,V_test,
                eval.times,fit.times,
                c_library,d_library,k.x.m_train){
  res<-NULL
  b.pred_temp=matrix(0,length(A_test),length(fit.times))
  Y_train=as.numeric(delta_train/k.x.m_train)*obs_N_train
  V_ALZ_train=cbind.data.frame(A=A_train,V=V_train,L=L_train,Z=Z_train)
  V_ALZ_test=cbind.data.frame(A=A_test,V=V_test,L=L_test,Z=Z_test)
  idx.d <- delta_train == 1
  for (u in 1:length(fit.times)){
    idx.c <- X_train> fit.times[u]
    crit.c <- (length(unique(Y_train[idx.c & idx.d])) == 1) | (sum(idx.c & idx.d) < 10)
    crit.d <- (length(unique(delta_train[idx.c])) == 1) | (sum(idx.c) < 10)
    
    if (!crit.c) {
      fit.c <- SuperLearner(
        Y = Y_train[idx.c & idx.d],          # I(X>u) = 1 and Delta = 1
        X = V_ALZ_train[idx.c & idx.d,],
        family = gaussian(),
        newX = V_ALZ_test,
        method = "method.NNLS",
        SL.library = c_library
      )
      c.pred=fit.c$SL.predict
    }
    
    if (!crit.d) {
      # fit the d model
      fit.d <- SuperLearner(
        Y = delta_train[idx.c],              # I(X>u) = 1
        X = V_ALZ_train[idx.c,],
        family = binomial(),
        newX = V_ALZ_test,
        method = "method.CC_LS",
        SL.library = d_library
      )
      d.pred <- fit.d$SL.predict
    }
    b.pred_temp[,u]=c.pred* d.pred
  }
  res$b.hat<- t(sapply(1:length(A_test), function(i) {
    approx(x = fit.times, y = b.pred_temp[i,], xout = eval.times, rule = 2)$y
  }))
  return(res)
}

kfolds=5


V=cbind.data.frame(V1,V2)
t_fits<-landmark_times
i_x_t=X > t_fits[t]

df <- data.frame(delta, A, Z, L)
df <- df %>%
  group_by(delta, A, Z, L) %>%
  mutate(idx = sample(rep(1:kfolds, length = n()))) %>%
  ungroup()

idx <- df$idx

#### integral part
u_grid_prime <- seq(0, tau, by = incr)
u_grid <- sort(unique(c(u_grid_prime,X)))
fit.times = quantile(u_grid_prime, probs = seq(0, 0.95, by = 0.05))
empt <- matrix(0, nrow = n, ncol = length(u_grid))
eval.times=u_grid

pi_LZ_fit=rep(0,times=n)
delta_A_fit=rep(0,times=n)
k.hat=matrix(0,n,length(eval.times))
h.hat=matrix(0,n,length(eval.times))
b.hat=matrix(0,n,length(eval.times))
delta_N_fit=rep(0,times=n)
delta_D_fit=rep(0,times=n)
mu_A_fit=rep(0,times=n)
mu_N_fit=rep(0,times=n)
mu_D_fit=rep(0,times=n)

for (fold in 1:kfolds) {
  # training and test indexes
  idx_test <- idx == fold
  idx_trn <- idx != fold
  
  pi_LZ_fit[idx_test]=as.numeric(pi_LZ(L[idx_trn], Z[idx_trn], 
                                       V[idx_trn,], L[idx_test], Z[idx_test], V[idx_test,]))
  delta_A_fun<-delta_A(A[idx_trn], L[idx_trn], Z[idx_trn], 
                       L[idx_test], Z[idx_test],
                       V[idx_trn,], V[idx_test,])
  delta_A_fit[idx_test]=delta_A_fun$delta_A
  mu_A_fit[idx_test]=delta_A_fun$mu_A
  
  surv_res<-surv_fit(X[idx_trn],delta[idx_trn],
                     A[idx_trn],A[idx_test],L[idx_trn],L[idx_test],
                     Z[idx_trn],Z[idx_test],V[idx_trn,],V[idx_test,],
                     cens_library=c("survSL.rfsrc", "survSL.coxph"),
                     event_library=c("survSL.rfsrc", "survSL.coxph"),
                     eval.times,eps = 1e-06)
  
  k.hat[idx_test,]=surv_res$k.hat
  h.hat[idx_test,]=surv_res$h.hat
  delta_N_res<-delta_N_tilde(obs_N[idx_trn,t], delta[idx_trn], 
                             V[idx_trn,], L[idx_trn], Z[idx_trn], 
                             V[idx_test,],L[idx_test], Z[idx_test],
                             k.x.m_train=surv_res$k.x.m,
                             sl_lib = c("SL.glm", "SL.lgb"))
  delta_N_fit[idx_test]<-delta_N_res$delta_N_tilde
  mu_N_fit[idx_test]<-delta_N_res$mu_N
  
  delta_D_res<-delta_D(i_x_t[idx_trn], delta[idx_trn], 
                       V[idx_trn,], L[idx_trn], Z[idx_trn], 
                       V[idx_test,],L[idx_test], Z[idx_test],
                       k.x.m_train=surv_res$k.x.m,
                       sl_lib = c("SL.glm", "SL.lgb"))
  delta_D_fit[idx_test]<-delta_D_res$delta_D
  mu_D_fit[idx_test]<-delta_D_res$mu_D
  
  b_res<-b_fit(obs_N[idx_trn,t],X[idx_trn],delta[idx_trn],
               A[idx_trn],A[idx_test],L[idx_trn],L[idx_test],
               Z[idx_trn],Z[idx_test],V[idx_trn,],V[idx_test,],
               eval.times,fit.times,
               c_library=c("SL.glm", "SL.lgb"),
               d_library=c("SL.glm", "SL.lgb"),
               k.x.m_train=surv_res$k.x.m)
  
  b.hat[idx_test,]=b_res$b.hat
}

nu <- length(eval.times)
eps = 1e-06
int.f <- sapply(1:n, function(i) {
  vals <- diff(-log(k.hat[i,])) * ((b.hat[i,-nu] + b.hat[i,-1])/2) #trapezoidal rule
  if(any(eval.times[-1] > X[i])) vals[eval.times[-1] > X[i]] <- 0 #only integrate until x
  sum(vals)
})

b.x <- sapply(1:n, function(i) stepfun(eval.times, c(b.hat[i,1],b.hat[i,]), right = FALSE)(X[i]))

k.x <- sapply(1:n, function(i) stepfun(eval.times, c(1,k.hat[i,]), right = FALSE)(X[i]))

k.x.m <- sapply(1:n, function(i) stepfun(eval.times, c(1,k.hat[i,]), right = FALSE)(X[i]-eps))
h.x <- sapply(1:n, function(i) stepfun(eval.times, c(1,h.hat[i,]), right = FALSE)(X[i]))
h.x.v.t <- sapply(1:n, function(i) {
  sf <- stepfun(eval.times, c(1, h.hat[i,]), right = FALSE)
  return(sf(max(t_fits[t], X[i])))
})

uvt <- sapply(1:nu, function(k) {max(eval.times[k], t_fits[t])})
int.h <- sapply(1:n, function(i) {
  huvt <- stepfun(eval.times, c(1,h.hat[i,]), right = FALSE)(uvt)
  vals <- diff(1/k.hat[i,]) * ((huvt[-nu] / h.hat[i,-nu] + huvt[-1] / h.hat[i,-1])/2)
  if(any(eval.times[-1] > X[i])) vals[eval.times[-1] > X[i]] <- 0
  sum(vals)
})

part_i=(delta_N_fit/delta_A_fit)*(2*A-1)
part_ii=((2*Z-1)*(2*L-1)/(pi_LZ_fit*delta_A_fit)*(delta/k.x.m*obs_N[,t]-mu_N_fit))*(2*A-1)
part_iii=((2*Z-1)*(2*L-1)/(pi_LZ_fit*delta_A_fit)*(delta_N_fit/delta_A_fit*(A-mu_A_fit)))*(2*A-1)
part_iv_i=((2*Z-1)*(2*L-1)/(pi_LZ_fit*delta_A_fit)*((1-delta)*b.x))*(2*A-1)
part_iv_ii=((2*Z-1)*(2*L-1)/(pi_LZ_fit*delta_A_fit)*(int.f))*(2*A-1)

aipw_n<-part_i+part_ii-part_iii+part_iv_i-part_iv_ii

part_i=(delta_D_fit/delta_A_fit)*(2*A-1)
part_ii=((2*Z-1)*(2*L-1)/(pi_LZ_fit*delta_A_fit)*(delta/k.x.m*i_x_t-mu_D_fit))*(2*A-1)
part_iii=((2*Z-1)*(2*L-1)/(pi_LZ_fit*delta_A_fit)*(delta_D_fit/delta_A_fit*(A-mu_A_fit)))*(2*A-1)
part_iv_i=((2*Z-1)*(2*L-1)/(pi_LZ_fit*delta_A_fit)*((1-delta) * h.x.v.t)/(h.x *k.x))*(2*A-1)
part_iv_ii=((2*Z-1)*(2*L-1)/(pi_LZ_fit*delta_A_fit)*(int.h))*(2*A-1)

aipw_d<-part_i+part_ii-part_iii+part_iv_i-part_iv_ii

##
V_hat_Nobs<-function(d){
  d_V=as.numeric(c(cbind(1, V1, V2) %*% d) > 0)
  constraint_val=mean((as.numeric(d_V==A)-as.numeric(d1_V==A))*aipw_d)
  if (constraint_val <= 0) {
    # 1e6 is a large penalty, + squared deviation pulls it back
    return(1e6 + (constraint_val)^2) 
  }else{
    return(mean(as.numeric(d_V==A)*aipw_n))
  }
}
V_opt_Nobs<- genoud(V_hat_Nobs, default.domains = 1000, nvars = 3,
                    pop.size = 6000, max = FALSE, print.level = 0,
                    wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
                    P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
                    P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
                    starting.values = NULL, hard.generation.limit = FALSE,
                    solution.tolerance = 1e-04, optim.method = "Nelder-Mead")
final_d_Nobs <- V_opt_Nobs$par

V1_eval=rnorm(1e6)
V2_eval=rnorm(1e6)

final_d_Nstar=ref$d_hat
estimated_d_Nstar=as.numeric(c(cbind(1, V1_eval, V2_eval) %*% final_d_Nstar) > 0)

estimated_d_Nobs=as.numeric(c(cbind(1, V1_eval, V2_eval) %*% final_d_Nobs) > 0)

d_V=rep(1,n)
d1_V=rep(0,n)

V_N<-mean((as.numeric(d_V==A)-as.numeric(d1_V==A))*aipw_n)
Var_N<-mean(aipw_n^2*(as.numeric(d_V==A)-as.numeric(d1_V==A))^2)

res<-list(d=as.numeric(final_d_Nobs),accuracy=mean(estimated_d_Nobs==estimated_d_Nstar),V_10=V_N,Var_10=Var_N)

setwd("~/Survival_DID/d_tilde_fixed/n_1000/Results_AIPW_1")
saveRDS(res,file=paste0("AIPW_1_",seed,".rds"))
