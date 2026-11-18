reference_path <- function() file.path(getwd(),'Reference','reference.rds')
load_reference <- function() {
  path <- reference_path()
  if(!file.exists(path)) stop('Run Setup_reference.R once before submitting simulation jobs')
  ref <- readRDS(path)
  if(ref$n!=1000000L || !identical(ref$dgp_md5,unname(tools::md5sum('Data_generator.R')))) {
    stop('Reference does not match this DGP or the required sample size')
  }
  ref
}
reference_policy <- function(ref,w) as.integer(predict_ml(ref$model,w[c('V1','V2')])>=.5)

normal_quadrature <- function(order) {
  j <- seq_len(order-1L); beta <- j/sqrt(4*j*j-1)
  J <- matrix(0,order,order)
  J[cbind(j,j+1L)] <- J[cbind(j+1L,j)] <- beta
  e <- eigen(J,symmetric=TRUE)
  ix <- order(e$values)
  nodes <- 2*e$values[ix]
  weight <- 4*e$vectors[1,ix]^2*dnorm(nodes)/(pnorm(2)-pnorm(-2))
  list(nodes=nodes,weights=weight)
}
reference_truth <- function(model,order=1024L) {
  q <- normal_quadrature(order)
  answer <- matrix(0,2,2,dimnames=list(c('1','4'),c('N','Y')))
  for(start in seq(1L,order,by=32L)) {
    cols <- start:min(order,start+31L)
    w <- expand.grid(V1=q$nodes,V2=q$nodes[cols])
    weight <- as.vector(outer(q$weights,q$weights[cols]))
    d <- as.integer(predict_ml(model,w)>=.5)
    for(t in c(1,4)) {
      cate <- idid_true_cate(w$V1,w$V2,t)
      answer[as.character(t),] <- answer[as.character(t),]+
        c(sum(weight*d*cate$tau_N),sum(weight*d*cate$tau_Y))
    }
  }
  answer
}
setup_reference <- function(cfg,force=FALSE) {
  if(file.exists(reference_path()) && !force) {
    message('Using existing frozen reference: ',reference_path())
    return(invisible(load_reference()))
  }
  message('Generating independent reference sample n=1,000,000')
  sim <- simulate_valid_idid(cfg$reference_n,seed=cfg$reference_seed)
  x <- sim$data[c('V1','V2')]; y <- sim$data$A
  rm(sim); invisible(gc())
  message('Fitting frozen ML behavioral reference')
  model <- fit_ml(y,x,'binomial',cfg$reference_fit_seed,cfg)
  rm(x,y); invisible(gc())
  message('Integrating reference policy under the new DGP')
  truth <- reference_truth(model,cfg$truth_quadrature)
  check <- reference_truth(model,cfg$truth_quadrature%/%2L)
  ref <- list(model=model,n=cfg$reference_n,seed=cfg$reference_seed,
              fit_seed=cfg$reference_fit_seed,threshold=.5,truth=truth,
              quadrature_order=cfg$truth_quadrature,quadrature_check=check,
              quadrature_difference=truth-check,
              dgp_md5=unname(tools::md5sum('Data_generator.R')),
              nuisance_settings=cfg,created=as.character(Sys.time()))
  atomic_rds(ref,reference_path(),compress=TRUE)
  write.csv(data.frame(t=c(1,4),ref_N=truth[,'N'],ref_Y=truth[,'Y'],
                      N_quadrature_difference=(truth-check)[,'N'],
                      Y_quadrature_difference=(truth-check)[,'Y']),
            file.path('Reference','truth.csv'),row.names=FALSE)
  message('Frozen reference saved; maximum quadrature difference: ',max(abs(truth-check)))
  invisible(ref)
}
true_contrast <- function(theta,t,ref) {
  value <- idid_true_policy_contrast(theta=theta,t=t)
  c(value_contrast=unname(value[1]-ref$truth[as.character(t),'N']),
    survival_contrast=unname(value[2]-ref$truth[as.character(t),'Y']))
}
