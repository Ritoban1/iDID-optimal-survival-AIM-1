####
landmark_times = seq(1,5,by=1)
#####

library(simsurv)


simulate_weibull_ph_times <- function(X, beta, weibull_shape, weibull_scale) {
  
  # n is determined from the design matrix
  n <- nrow(X)
  
  # Calculate the linear predictor
  # This is the 'eta' term
  eta <- X %*% beta
  
  # Generate U ~ Unif(0, 1)
  U <- runif(n, 0, 1)
  
  # Calculate times using the inversion formula
  # T = scale * ( (-log(U)) / exp(eta) )^(1 / shape)
  times <- weibull_scale * ( (-log(U)) / exp(eta) )^(1 / weibull_shape)
  
  return(times)
}
V1 <- rnorm(n)
V2 <- rnorm(n)
Z <- rbinom(n, size = 1, prob = 0.5)
L <- rbinom(n, size = 1, prob = 0.5)

U0 <- rbridge(n)
U1 <- rbridge(n)

A0 <- rbinom(n, size = 1, prob = plogis(2 - 7 * Z + 0.2 * U0 + 2 * V1))
A1 <- rbinom(n, size = 1, prob = plogis(-1.5 + 5 * Z - 0.15 * U1 + 1.5 * V2))
A <- L * A1 + (1 - L) * A0

###
beta_C=c(-3,0.3,0.1,0.1,0.2,0.1)

X_C <- cbind(
  intercept = 1,
  A = A,
  L = L,
  Z = Z,
  V1 =V1,
  V2 = V2)

C <- simulate_weibull_ph_times(X_C, beta = beta_C, weibull_shape = 1 , weibull_scale = 2) 

### Death
lambda_0 <- 0.2 + ((-0.1 + 0.2*V1 -0.2*V2)/10 * A0) + (-0.01 * V1) + (0.02 * V2) + (0.03 * U0) + (0.01 * Z) 

if(any(lambda_0 <= 0)) lambda_0[lambda_0 <= 0] <- as.numeric(quantile(lambda_0,probs = 0.25))
T_event_0 <- rexp(n, rate = lambda_0)

lambda_1 <- 0.25 + ((-0.1 + 0.2*V1 -0.2*V2)/10  * A1) + (-0.01 * V1) + (0.02 * V2) + (0.03 * U1) + (0.01 * Z) 

if(any(lambda_1 <= 0)) lambda_1[lambda_1 <= 0] <- as.numeric(quantile(lambda_1,probs = 0.25))

T_event_1 <- rexp(n, rate = lambda_1)

T_event= T_event_1*L + (1-L)*T_event_0
X=apply(cbind(C, T_event),1,min)
delta=ifelse(C>= T_event,1,0)

# Initialize cumulative counts (at t=0, count is 0)
current_N_0 <- matrix(0,n,length(landmark_times))
current_N_1 <-  matrix(0,n,length(landmark_times))
prev_time <- rep(0, n)

# Recurrent Event Baseline Parameters (Weibull)
wb_scale_rec <- 1.0
wb_shape_rec <- 1.1
time=1
for (t_k in landmark_times) {
  eff_end <- t_k
  # Effective start of interval
  eff_start <-prev_time
  # Duration of risk in this specific interval
  duration <- pmax(0, eff_end - eff_start)
  # 2. Calculate Expected Increment in this Interval
  # E[dN] = (Weibull_Cum(eff_end) - Weibull_Cum(eff_start)) + (Additive_Rate * duration)
  
  # A. Weibull Part (Baseline Cumulative Hazard)
  cum_base_end <- (eff_end / wb_scale_rec)^wb_shape_rec
  cum_base_start <- (eff_start / wb_scale_rec)^wb_shape_rec
  diff_base <- cum_base_end - cum_base_start
  
  # B. Additive Part (Covariates + Time Trend)
  # We assume the "Time Trend" effect depends on the landmark time t_k itself
  period_effect <- 0.1 * t_k
  
  # Additive Rate for Potential Outcome 0
  rate_0 <- period_effect + ((-0.1 + 0.2*V1 -0.2*V2)  * A0) + (0.05 * V1) + (0.03 * V2) + (0.05 * U0) + (0.02 * Z)
  
  # Additive Rate for Potential Outcome 1
  rate_1 <- period_effect + 0.2 + ((-0.1 + 0.2*V1 -0.2*V2) * A1) + (0.05 * V1) + (0.03 * V2) + (0.05 * U1) + (0.02 * Z)
  
  # Realized Rate
  
  # Expected additive contribution
  diff_additive_0 <- rate_0 * duration
  diff_additive_1 <- rate_1 * duration
  
  # Total expected increment
  expected_dN_0 <- diff_base +  diff_additive_0
  expected_dN_1 <- diff_base +  diff_additive_1
  
  # Safety: ensure non-negative
  if(any(expected_dN_0 < 0)) expected_dN_0[expected_dN_0 < 0] <- 0
  if(any(expected_dN_1 < 0)) expected_dN_1[expected_dN_1 < 0] <- 0
  
  # 3. Generate Incremental Count dN
  dN_0 <- rpois(n, lambda = expected_dN_0)
  dN_1 <- rpois(n, lambda = expected_dN_1)
  # Update Cumulative Count N(t)
  if(time==1){
    current_N_0[,time] <- dN_0
    current_N_1[,time] <- dN_1
  }else{
    current_N_0[,time] <- current_N_0[,(time-1)] + dN_0
    current_N_1[,time] <- current_N_1[,(time-1)] + dN_1
  }
  prev_time <- t_k
  time=time+1
}

star_obs<-function(current_N){
  N_star <- matrix(0, n, length(landmark_times))
  obs_N  <- matrix(0, n, length(landmark_times))
  
  # Pre-calculate the index of the last valid landmark for every subject
  # Find the largest index 'j' such that landmark_times[j] <= T_event
  # We can use findInterval
  idx_star <- findInterval(T_event, landmark_times)
  idx_obs  <- findInterval(X, landmark_times)
  
  for (k in 1:length(landmark_times)) {
    t_k <- landmark_times[k]
    
    # --- 1. N_star (Stopped at T_event) ---
    # If current landmark t_k <= T_event, we haven't stopped yet: use count at t_k
    # If current landmark t_k > T_event, we stopped earlier: use count at idx_star
    
    # Alive part: direct lookup
    alive_mask <- t_k <= T_event
    N_star[alive_mask, k] <- current_N[alive_mask, k]
    
    # Stopped part: use the last valid count before death
    stopped_mask <- !alive_mask
    # Identify which index to pull from for the stopped subjects
    # If idx_star is 0 (died before first landmark), count is 0.
    valid_stops <- stopped_mask & (idx_star > 0)
    
    if (any(valid_stops)) {
      # We need to pull specific column for specific rows. 
      # Since R matrices are column-major, simple indexing is tricky.
      # Easier to loop or use matrix indexing.
      # Here we use the vector of indices idx_star[valid_stops]
      rows <- which(valid_stops)
      cols <- idx_star[rows]
      # Linear index in matrix
      lin_idx <- cbind(rows, cols)
      N_star[rows, k] <- current_N[lin_idx]
    }
    # Note: N_star initialized to 0, so cases with idx_star == 0 are already handled.
    
    
    # --- 2. obs_N (Stopped at X) ---
    # Same logic using X and idx_obs
    
    # At risk part
    risk_mask <- t_k <= X
    obs_N[risk_mask, k] <- current_N[risk_mask, k]
    
    # Censored/Event part
    stopped_mask_obs <- !risk_mask
    valid_stops_obs <- stopped_mask_obs & (idx_obs > 0)
    
    if (any(valid_stops_obs)) {
      rows <- which(valid_stops_obs)
      cols <- idx_obs[rows]
      lin_idx <- cbind(rows, cols)
      obs_N[rows, k] <- current_N[lin_idx]
    }
  }
  return(list(N_star=N_star,obs_N=obs_N))
}

current_N=current_N_1*L + (1-L)*current_N_0
star_obs_obs=star_obs(current_N)

N_star=star_obs_obs$N_star
obs_N=star_obs_obs$obs_N






