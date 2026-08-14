library(bridgedist)
library(rgenoud)
library(simsurv)
library(survival)
library(nloptr)
set.seed(1000)
###
n=1e6
t=1
source("Data_generator.R")


### pi(L,Z|V)
prob_L=mean(L)
model_Z<-glm(Z~V1+V2,family="binomial")
pi_LZ=as.numeric((Z*model_Z$fitted.values+(1-Z)*(1-model_Z$fitted.values))*prob_L)

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

### d1_V
model_d1<-glm(A~V1+V2,family="binomial")
pred_d1=predict(model_d1,type="response")
d1_V=as.numeric(ifelse(pred_d1>=0.5,1,0))


##
V_hat_Nstar<-function(d){
  d_V=as.numeric(c(cbind(1, V1, V2) %*% d) > 0)
  constraint_val=mean((as.numeric(d_V==A)-as.numeric(d1_V==A))*(2 * Z - 1) * 
                        (2 * L - 1)* (2 * A-1)*(as.numeric(T_event>t))/ 
                        (1/4* delta_A))
  if (constraint_val <= 0) {
    # 1e6 is a large penalty, + squared deviation pulls it back
    return(1e6 + (constraint_val)^2) 
  }else{
    return(mean(as.numeric(d_V==A)*(2 * Z - 1) * (2 * L - 1)* (2*A-1)*N_star[, t]/ 
                  (1/4* delta_A)))
  }
}

V_opt_Nstar<- genoud(V_hat_Nstar, default.domains = 1000, nvars = 3,
                     max = FALSE, print.level = 0,
                     wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
                     P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
                     P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
                     starting.values = NULL, hard.generation.limit = FALSE,
                     solution.tolerance = 1e-04, optim.method = "Nelder-Mead")
final_d_Nstar <- V_opt_Nstar$par

V_10=mean((as.numeric(A==1)-as.numeric(A==0))*(2 * Z - 1) * (2 * L - 1)* (2*A-1)*N_star[, t]/ 
            (1/4* delta_A))

d_V=as.numeric(c(cbind(1, V1, V2) %*% final_d_Nstar) > 0)

Value=mean((as.numeric(A==d_V)-as.numeric(A==d1_V))*(2 * Z - 1) * (2 * L - 1)* (2*A-1)*N_star[, t]/ 
             (1/4* delta_A))
#estimated_d_Nstar=as.numeric(c(cbind(1, V1, V2) %*% final_d_Nstar) > 0)
saveRDS(list(d_hat=final_d_Nstar,V_10=V_10,Value=Value,d1_V=d1_V,model_d1=model_d1),"result_IPW_large_1.rds")

