make_folds <- function(data,seed,k=5L) {
  set.seed(seed)
  group <- interaction(data$L,data$Z,data$delta,drop=TRUE)
  f <- integer(nrow(data))
  for(ix in split(seq_len(nrow(data)),group)) f[ix] <- {
    labels <- rep(seq_len(k), length.out = length(ix))
    labels[sample.int(length(ix))]
  }
  if (length(unique(f)) < k) {
    ord <- sample.int(nrow(data))
    f[ord] <- rep(seq_len(k), length.out = nrow(data))
  }
  f
}
cell_predictions <- function(fit,w,cfg) {
  nd <- do.call(rbind,lapply(0:3,function(j) {
    data.frame(L=j%/%2,Z=j%%2,w[cfg$wcols],check.names=FALSE)
  }))
  matrix(predict_ml(fit,nd),nrow=nrow(w),ncol=4)
}
cell_delta <- function(p) p[,4]-p[,3]-p[,2]+p[,1]
cell_observed <- function(p,d) p[cbind(seq_len(nrow(d)),1+2*d$L+d$Z)]
glm_probability <- function(y,x,newx) {
  if(!length(y)) stop('Every period-instrument cell needs training observations for cellwise logistic fitting.')
  if(length(unique(y))==1L) return(rep(y[1],nrow(newx)))
  fit <- glm(y~.,data=data.frame(y=y,x),family=binomial(),control=glm.control(maxit=100))
  as.numeric(predict(fit,newdata=newx,type='response'))
}
get_base <- function(sim,seed,cfg) {
  parametric <- cfg$propensity == 'logistic'
  cached('base',cfg,function() {
    d <- sim$data; n <- nrow(d); folds <- if(is.null(cfg$fold_id)) make_folds(d,seed+200000L,cfg$outer_folds) else cfg$fold_id
    out <- list(pi=numeric(n),muA=numeric(n),delta=numeric(n),folds=folds,
                fold=list(),diagnostics=list(delta_clipped=0L,pi_clipped=0L))
    for(f in seq_len(cfg$outer_folds)) {
      tr <- which(folds!=f); te <- which(folds==f); train <- d[tr,]; test <- d[te,]
      s <- seed+f*10000L
      if(parametric) {
        lp <- glm_probability(train$L,train[cfg$wcols],test[cfg$wcols])
        zp <- glm_probability(train$Z,train[c('L',cfg$wcols)],test[c('L',cfg$wcols)])
        pa <- matrix(0,length(te),4)
        for(j in 0:3) {
          use <- train$L==j%/%2 & train$Z==j%%2
          pa[,j+1] <- glm_probability(train$A[use],train[use,cfg$wcols],test[cfg$wcols])
        }
      } else {
        fitL <- fit_ml(train$L,train[cfg$wcols],'binomial',s+1L,cfg)
        lp <- predict_ml(fitL,test[cfg$wcols])
        fitZ <- fit_ml(train$Z,train[c('L',cfg$wcols)],'binomial',s+2L,cfg)
        zp <- predict_ml(fitZ,test[c('L',cfg$wcols)])
        fitA <- fit_ml(train$A,train[cfg$xcols],'binomial',s+3L,cfg)
        pa <- cell_predictions(fitA,test,cfg)
      }
      surv <- fit_survival(train,test,'censor',s+4L,cfg)
      out$diagnostics$pi_clipped <- out$diagnostics$pi_clipped +
        sum(lp<cfg$prob_floor | lp>1-cfg$prob_floor | zp<cfg$prob_floor | zp>1-cfg$prob_floor)
      lp <- pmin(1-cfg$prob_floor,pmax(cfg$prob_floor,lp))
      zp <- pmin(1-cfg$prob_floor,pmax(cfg$prob_floor,zp))
      out$pi[te] <- ifelse(test$L==1,lp,1-lp)*ifelse(test$Z==1,zp,1-zp)
      if(any(!is.finite(pa))) stop('Nonfinite treatment probabilities')
      pa <- pmin(pmax(pa,0),1)
      raw <- cell_delta(pa)
      out$diagnostics$delta_clipped <- out$diagnostics$delta_clipped+sum(cfg$first_stage_sign*raw<cfg$delta_floor)
      out$diagnostics$min_delta_raw <- min(c(out$diagnostics$min_delta_raw,raw))
      out$delta[te] <- cfg$first_stage_sign*pmax(cfg$delta_floor,cfg$first_stage_sign*raw)
      out$muA[te] <- cell_observed(pa,test)
      out$fold[[f]] <- list(train=tr,test=te,K=surv)
      if(cfg$verbose) message('  Base fold ',f,'/',cfg$outer_folds,if(parametric) ' parametric' else ' ML')
    }
    out$diagnostics$min_delta <- min(out$delta)
    out$diagnostics$min_pi <- min(out$pi)
    out$diagnostics$min_K <- min(vapply(out$fold,function(z) min(c(1,z$K$test)),numeric(1)))
    out$diagnostics$K_test_floored <- sum(vapply(out$fold,function(z) sum(z$K$diagnostics$test_floored),numeric(1)))
    out$diagnostics$K_training_floored <- sum(vapply(out$fold,function(z) sum(z$K$diagnostics$train_oob_floored),numeric(1)))
    out
  })
}
get_death <- function(sim,base,seed,cfg) {
  cached('death',cfg,function() {
    lapply(seq_along(base$fold),function(f) {
      ix <- base$fold[[f]]
      fit_survival(sim$data[ix$train,],sim$data[ix$test,],'death',seed+f*10000L+5L,cfg)
    })
  })
}
b_knots <- function(train, landmarks, cfg) {
  end <- min(cfg$tau, max(train$X))
  n_knots <- max(cfg$b_points, ceiling(3 * nrow(train)^(1/4)))
  sort(unique(c(seq(0, end, length.out = n_knots), landmarks[landmarks <= end])))
}

b_risk <- function(train, u, end) {
  if (u == end) which(train$X >= end) else which(train$X > u)
}

get_completion <- function(sim,base,seed,cfg) {
  # e(V) and d(u,Q) do not depend on the count landmark. Cache per DATASET/fold,
  # so t=1 and t=4 share predictions without sharing training/evaluation data.
  cached('completion',cfg,function() {
    qcols <- cfg$qcols
    xcols <- cfg$xcols
    out <- list(fold=list(),diagnostics=list(b_d_small_risksets=0L,
                  b_d_constant_fits=0L,b_d_zero_completion_risksets=0L,
                  mu_e_constant_fits=0L))
    for(f in seq_along(base$fold)) {
      ix <- base$fold[[f]]
      train <- sim$data[ix$train,]; test <- sim$data[ix$test,]
      # e(V) for the two marginal outcome means excludes A. It is shared
      # across N/Y and t, but is NOT the A-conditional d used for b below.
      efit <- fit_ml(train$delta,train[xcols],'binomial',
                     seed+f*10000L+500000L,cfg)
      if(identical(efit$type,'constant')) out$diagnostics$mu_e_constant_fits <- out$diagnostics$mu_e_constant_fits+1L
      E <- cell_predictions(efit,test,cfg)
      knots <- b_knots(train,sim$b_times,cfg)
      D <- matrix(NA_real_,nrow(test),length(knots))
      for(j in seq_along(knots)) {
        risk <- b_risk(train,knots[j],max(knots))
        if(!length(risk)) stop('Empty risk set for b at time ',knots[j])
        if(length(risk)<10L) out$diagnostics$b_d_small_risksets <- out$diagnostics$b_d_small_risksets+1L
        if(!any(train$delta[risk]==1L)) out$diagnostics$b_d_zero_completion_risksets <- out$diagnostics$b_d_zero_completion_risksets+1L
        dfit <- fit_ml(train$delta[risk],train[risk,qcols,drop=FALSE],
                       'binomial',seed+f*10000L+600000L+j,cfg)
        if(identical(dfit$type,'constant')) out$diagnostics$b_d_constant_fits <- out$diagnostics$b_d_constant_fits+1L
        D[,j] <- pmin(1,pmax(0,predict_ml(dfit,test[qcols])))
      }
      out$fold[[f]] <- list(knots=knots,D=D,E=E)
      if(cfg$verbose) message('  AIPW completion fold ',f,'/',length(base$fold),' (shared across landmarks)')
    }
    out
  })
}
fit_b_counts <- function(train,train_RN,test,knots,D,seed,cfg) {
  # c(u,t,Q) = E[N(t)/K(X-) | delta=1, X>u, Q]. Completion includes
  # completed ZERO counts. This is not a positive-count hurdle regression.
  qcols <- cfg$qcols
  B <- matrix(NA_real_,nrow(test),length(knots))
  diag <- list(b_small_risksets=0L,b_constant_fits=0L,
               b_c_small_risksets=0L,b_c_empty_risksets=0L,
               b_completed_rows_min=Inf)
  for(j in seq_along(knots)) {
    risk <- b_risk(train,knots[j],max(knots))
    if(!length(risk)) stop('Empty risk set for b at time ',knots[j])
    if(length(risk)<10L) diag$b_small_risksets <- diag$b_small_risksets+1L
    complete <- risk[train$delta[risk]==1L]
    diag$b_completed_rows_min <- min(diag$b_completed_rows_min,length(complete))
    if(length(complete)<10L) diag$b_c_small_risksets <- diag$b_c_small_risksets+1L
    if(!length(complete)) {
      # d is identically zero here; c is unidentified but its product is zero.
      if(any(D[,j]!=0)) stop('Nonzero completion prediction in an all-censored risk set')
      B[,j] <- 0
      diag$b_c_empty_risksets <- diag$b_c_empty_risksets+1L
      next
    }
    cfit <- fit_ml(train_RN[complete],train[complete,qcols,drop=FALSE],
                   'gaussian',seed+20L+j,cfg)
    if(identical(cfit$type,'constant')) diag$b_constant_fits <- diag$b_constant_fits+1L
    B[,j] <- pmax(0,predict_ml(cfit,test[qcols]))*D[,j]
  }
  list(B=B,diagnostics=diag)
}
interpolate_b <- function(B,times,q) {
  # Linear interpolation applies ONLY to the separately estimated b regression.
  # Survival functions are always evaluated as exact step functions.
  n <- nrow(B); q <- rep_len(q,n)
  j <- pmax(1L,pmin(length(times)-1L,findInterval(q,times)))
  w <- pmin(1,pmax(0,(q-times[j])/(times[j+1L]-times[j])))
  (1-w)*B[cbind(seq_len(n),j)]+w*B[cbind(seq_len(n),j+1L)]
}
get_scores <- function(sim,t,method,seed,cfg) {
    base <- get_base(sim,seed,cfg)
    d <- sim$data; n <- nrow(d); col <- match(t,sim$landmark_times)
    if(is.na(col)) stop('Unsupported landmark')
    ans <- list(N=numeric(n),Y=numeric(n),diagnostics=base$diagnostics,folds=base$folds)
    if(method=='AIPW') {
      death <- get_death(sim,base,seed,cfg)
      completion <- get_completion(sim,base,seed,cfg)
      ans$diagnostics <- c(ans$diagnostics,completion$diagnostics)
      ans$diagnostics$min_H <- min(vapply(death,function(z) min(c(1,z$test)),numeric(1)))
      ans$diagnostics$H_test_floored <- sum(vapply(death,function(z) sum(z$diagnostics$test_floored),numeric(1)))
      ans$diagnostics$mu_completed_rows_min <- Inf
      ans$diagnostics$mu_N_constant_fits <- 0L
      ans$diagnostics$mu_Y_constant_fits <- 0L
    }
    ans$diagnostics$b_small_risksets <- 0L
    ans$diagnostics$b_constant_fits <- 0L
    ans$diagnostics$b_c_small_risksets <- 0L
    ans$diagnostics$b_c_empty_risksets <- 0L
    ans$diagnostics$b_completed_rows_min <- Inf
    ans$diagnostics$b_max_knots <- 0L
    ans$diagnostics$b_endpoint_rows_min <- Inf
    ans$diagnostics$b_tail_extrapolation_folds <- 0L
    for(f in seq_len(cfg$outer_folds)) {
      ix <- base$fold[[f]]; tr <- ix$train; te <- ix$test
      train <- d[tr,]; test <- d[te,]; K <- ix$K
      kx <- curve_lookup(K$test,K$times,test$X,left=TRUE)
      w <- (2*test$Z-1)*(2*test$L-1)/(base$pi[te]*base$delta[te])
      RN <- test$delta*sim$obs_N[te,col]/kx
      RY <- test$delta*as.integer(test$X>=t)/kx
      if(method=='IPW1') {
        ans$N[te] <- w*RN; ans$Y[te] <- w*RY
        next
      }
      if(method=='IPW2') {
        ans$N[te] <- w * weighted_events(sim$events[te], K, t, cfg$event_grid)
        ans$Y[te] <- w*as.integer(test$X>=t)/curve_lookup(K$test,K$times,t,left=TRUE)
        next
      }
      if(method!='AIPW') stop('Unknown method: ',method)
      s <- seed+f*10000L+as.integer((t*100) %% 100000L)
      train_RN <- train$delta*sim$obs_N[tr,col]/K$train_oob_at_X
      train_RY <- train$delta*as.integer(train$X>=t)/K$train_oob_at_X
      xcols <- cfg$xcols
      complete <- which(train$delta==1L)
      if(!length(complete)) stop('No completed training histories for outcome means')
      fitN <- fit_ml(train_RN[complete],train[complete,xcols,drop=FALSE],
                     'gaussian',s+10L,cfg)
      fitY <- fit_ml(train_RY[complete],train[complete,xcols,drop=FALSE],
                     'gaussian',s+11L,cfg)
      ans$diagnostics$mu_completed_rows_min <- min(ans$diagnostics$mu_completed_rows_min,length(complete))
      ans$diagnostics$mu_N_constant_fits <- ans$diagnostics$mu_N_constant_fits+as.integer(identical(fitN$type,'constant'))
      ans$diagnostics$mu_Y_constant_fits <- ans$diagnostics$mu_Y_constant_fits+as.integer(identical(fitY$type,'constant'))
      # Both responses can exceed one after weighting; neither is binomial.
      # Keep completed zero outcomes, including deaths before t for Y.
      e <- completion$fold[[f]]$E
      pn <- pmax(0,cell_predictions(fitN,test,cfg))*e
      py <- pmax(0,cell_predictions(fitY,test,cfg))*e
      muN <- cell_observed(pn,test); muY <- cell_observed(py,test)
      tauN <- cell_delta(pn)/base$delta[te]; tauY <- cell_delta(py)/base$delta[te]

      knots <- completion$fold[[f]]$knots
      ans$diagnostics$b_tail_extrapolation_folds <- ans$diagnostics$b_tail_extrapolation_folds + as.integer(max(knots)<cfg$tau)
      ans$diagnostics$b_max_knots <- max(ans$diagnostics$b_max_knots,length(knots))
      ans$diagnostics$b_endpoint_rows_min <- min(ans$diagnostics$b_endpoint_rows_min,sum(train$X>=cfg$tau))
      bfit <- fit_b_counts(train,train_RN,test,knots,completion$fold[[f]]$D,s,cfg)
      B <- bfit$B
      for(nm in names(bfit$diagnostics)) {
        ans$diagnostics[[nm]] <- if(nm=='b_completed_rows_min') {
          min(ans$diagnostics[[nm]],bfit$diagnostics[[nm]])
        } else ans$diagnostics[[nm]]+bfit$diagnostics[[nm]]
      }
      use <- which(K$times<cfg$tau)
      times <- K$times[use]; kp <- K$test[,use,drop=FALSE]
      kprev <- if(length(use)) cbind(1,K$test[,-ncol(K$test),drop=FALSE])[,use,drop=FALSE] else kp
      hazard <- pmax(1-kp/kprev,0)
      Ydag <- outer(test$X,times,'>') |
        (outer(test$X,times,'==') & (test$J==0L))
      intN <- numeric(length(te))
      for(j in seq_along(times)) {
        bj <- interpolate_b(B,knots,times[j])
        intN <- intN+bj*Ydag[,j]*hazard[,j]
      }
      augN <- (1-test$delta)*interpolate_b(B,knots,test$X)-intN

      H <- death[[f]]
      hu <- curve_eval(H$test,H$times,times,left=FALSE)
      ht <- curve_lookup(H$test,H$times,t,left=TRUE)
      numerator <- hu
      if(any(times<t)) numerator[,times<t] <- ht
      BY <- numerator/(hu*kp)
      intY <- rowSums(BY*Ydag*hazard)
      hx <- curve_lookup(H$test,H$times,test$X,left=FALSE)
      hnumx <- ifelse(test$X<t,ht,hx)
      kpostx <- curve_lookup(K$test,K$times,test$X,left=FALSE)
      augY <- (1-test$delta)*hnumx/(hx*kpostx)-intY
      ans$N[te] <- tauN+w*(RN-muN-tauN*(test$A-base$muA[te])+augN)
      ans$Y[te] <- tauY+w*(RY-muY-tauY*(test$A-base$muA[te])+augY)
      if(cfg$verbose) message('  AIPW outcome fold ',f,'/',cfg$outer_folds,' t=',t)
    }
    if(any(!is.finite(ans$N)) || any(!is.finite(ans$Y))) stop('Nonfinite effect scores')
    ans
}
