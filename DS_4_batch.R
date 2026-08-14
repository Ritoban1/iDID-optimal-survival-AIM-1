# pass arguments to file
args = commandArgs(trailingOnly = T)

# takes first argument that was passed
seed = as.numeric(args[1])
# set seed (later name file according to seed)
set.seed(seed)

library(bridgedist)
library(rgenoud)
library(simsurv)
library(survival)
library(nloptr)

n=1e3
tau=10
t_e=4
setwd("~/Survival_DID/d_tilde_fixed/")
source("Data_generator.R")

### pi(L,Z|V)
prob_L=mean(L)
model_Z<-glm(Z~V1+V2,family="binomial")
pi_LZ=as.numeric((Z*model_Z$fitted.values+(1-Z)*(1-model_Z$fitted.values))*prob_L)

### d1_V
ref=readRDS("result_IPW_large_4.rds")
model_d1<-ref$model_d1
pred_d1=predict(model_d1,newdata=cbind.data.frame(1,V1,V2),type="response")
d1_V=as.numeric(ifelse(pred_d1>=0.5,1,0))

####
A11=A[L==1 & Z==1]
A10=A[L==1 & Z==0]
A01=A[L==0 & Z==1]
A00=A[L==0 & Z==0]

V1_11=V1[L==1 & Z==1]
V1_10=V1[L==1 & Z==0]
V1_01=V1[L==0 & Z==1]
V1_00=V1[L==0 & Z==0]

V2_11=V2[L==1 & Z==1]
V2_10=V2[L==1 & Z==0]
V2_01=V2[L==0 & Z==1]
V2_00=V2[L==0 & Z==0]

model_A<-function(A_res,V1_cov,V2_cov){
  model<-glm(A_res ~ V1_cov + V2_cov, family = "binomial")
  p<-predict(model,data.frame(V1_cov=V1,V2_cov=V2),type="response")
  return(p)
}

delta_A=as.numeric(model_A(A11,V1_11,V2_11)-model_A(A10,V1_10,V2_10)-
                     model_A(A01,V1_01,V2_01)+model_A(A00,V1_00,V2_00))



# u grid
u_grid <- seq(0, tau, by = 0.01)

#interpolate N_E
ts <- landmark_times
N_E <- obs_N
N_E <- cbind(0, N_E)
ts <- union(c(0), ts)
dN_E <- matrix(0, nrow = n, ncol = length(u_grid))

for (t in 1:(length(ts)-1)) {
  for (i in 1:n) {
    dN_E[i,(u_grid >= ts[t]) & (u_grid < ts[t+1])] <- (N_E[i,t+1] - N_E[i,t])/(sum((u_grid >= ts[t]) & (u_grid < ts[t+1])))
  }
}

####
data_surv=data.frame(A,V1,V2,L,Z)
cens_model <- coxph(Surv(X, (1-delta)) ~ A + V1 + V2 + L + Z)
cum_hazard_C <- predict(cens_model, type = "expected")
prob_C_ge_D <- exp(-cum_hazard_C)
fit_individual <- survfit(cens_model, newdata = data_surv)
prob_at_times <- summary(fit_individual, times = u_grid)$surv
prob_at_times=t(prob_at_times)

###
DS_N <- sapply(1:n, function(i){
  (2 * Z[i] - 1) * (2 * L[i] - 1)* (2 * A[i]-1)*1/(pi_LZ[i] * delta_A[i])*cumsum(dN_E[i,]/(prob_at_times[i,]))
}) 

DS_D=(2 * Z - 1) * (2 * L - 1)* (2 * A-1)*delta/ (prob_C_ge_D*pi_LZ * delta_A)
IPW_D=DS_D*(as.numeric(X>t_e))


#mu_hat <- rowMeans(mu_hat, na.rm = TRUE) 
l_t_u=which(u_grid %in% t_e)
DS_N=DS_N[l_t_u,]

###

V_hat_Nobs<-function(d){
  d_V=as.numeric(c(cbind(1, V1, V2) %*% d) > 0)
  constraint_val=mean((as.numeric(d_V==A)-as.numeric(d1_V==A))*IPW_D)
  if (constraint_val <= 0) {
    # 1e6 is a large penalty, + squared deviation pulls it back
    return(1e6 + (constraint_val)^2) 
  }else{
    return(mean(as.numeric(d_V==A)*DS_N))
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
res<-list(d=as.numeric(final_d_Nobs),accuracy=mean(estimated_d_Nobs==estimated_d_Nstar))

setwd("~/Survival_DID/d_tilde_fixed/n_1000/Results_DS_4")
saveRDS(res,file=paste0("DS_4_",seed,".rds"))
