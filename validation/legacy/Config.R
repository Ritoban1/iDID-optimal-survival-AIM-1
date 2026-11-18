# Shared settings for every n, method, and landmark.
get_config <- function() {
  fast <- identical(Sys.getenv('IDID_FAST'), '1')
  x <- list(outer_folds=5L,inner_folds=3L,rf_trees=500L,rf_nodesize=15L,
            rf_sample_min=1000L,rf_sample_power=.75,
            lgb_rounds=200L,prob_floor=.01,delta_floor=.05,surv_floor=.02,
            b_points=21L,b_estimation='two_part',mu_estimation='two_part',
            threads=1L,tau=10,replications=500L,
            genoud_population=2000L,genoud_generations=40L,genoud_wait=6L,
            theta_bound=3,reference_n=1000000L,reference_seed=700000001L,
            reference_fit_seed=700000002L,truth_seed=700000003L,
            truth_quadrature=1024L,test_mode=fast)
  if(fast) {
    x$rf_trees <- 80L; x$lgb_rounds <- 60L; x$b_points <- 11L
    x$genoud_population <- 400L; x$genoud_generations <- 12L; x$genoud_wait <- 3L
  }
  code <- c('Config.R','Data_generator.R','ML_functions.R','Score_functions.R')
  h <- tools::md5sum(code[file.exists(code)])
  x$signature <- paste(c('idid-cluster-v2-two-part',unname(h),unlist(x)),collapse='|')
  x
}

dataset_seed <- function(n,replication,evaluation=FALSE) {
  as.integer(n*1000+replication+if(evaluation) 100000000 else 1000000)
}
result_root <- function() path.expand(Sys.getenv('IDID_RESULT_ROOT',getwd()))
cache_root <- function() path.expand(Sys.getenv('IDID_CACHE_ROOT',file.path(getwd(),'Cache')))
atomic_rds <- function(x,path,compress=FALSE) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp <- tempfile('.writing_',tmpdir=dirname(path))
  on.exit(unlink(tmp),add=TRUE)
  saveRDS(x,tmp,compress=compress)
  if(!file.rename(tmp,path)) stop('Could not save ',path)
  invisible(path)
}
cached <- function(path,cfg,compute) {
  if(file.exists(path)) {
    old <- tryCatch(readRDS(path),error=function(e) NULL)
    if(is.list(old) && identical(old$signature,cfg$signature)) return(old$value)
  }
  value <- compute()
  atomic_rds(list(signature=cfg$signature,value=value),path)
  value
}
