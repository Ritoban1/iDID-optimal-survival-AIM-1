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
t=1
n=1e3
setwd("~/Survival_DID/d_tilde_fixed/")
source("Data_generator.R")

### pi(L,Z|V)
prob_L=mean(L)
model_Z<-glm(Z~V1+V2,family="binomial")
pi_LZ=as.numeric((Z*model_Z$fitted.values+(1-Z)*(1-model_Z$fitted.values))*prob_L)

### d1_V
ref=readRDS("result_IPW_large_1.rds")
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

##

### censoring
cens_model <- coxph(Surv(X, (1-delta)) ~ A + V1 + V2 + L + Z)
cum_hazard_C <- predict(cens_model, type = "expected")
prob_C_ge_D <- exp(-cum_hazard_C)

###



comm=(2 * Z - 1) * (2 * L - 1)* (2 * A-1)*delta/ (prob_C_ge_D*pi_LZ * delta_A)
IPW_D=comm*(as.numeric(X>t))
IPW_N=comm*obs_N[,t]

##
V_hat_Nobs<-function(d){
  d_V=as.numeric(c(cbind(1, V1, V2) %*% d) > 0)
  constraint_val=mean((as.numeric(d_V==A)-as.numeric(d1_V==A))*IPW_D)
  if (constraint_val <= 0) {
    # 1e6 is a large penalty, + squared deviation pulls it back
    return(1e6 + (constraint_val)^2) 
  }else{
    return(mean(as.numeric(d_V==A)*IPW_N))
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

setwd("~/Survival_DID/d_tilde_fixed/n_1000/Results_IPW_1")
saveRDS(res,file=paste0("IPW_1_",seed,".rds"))



