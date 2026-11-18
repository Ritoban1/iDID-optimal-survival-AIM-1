# Every wrapper starts here; common code keeps all n/method implementations equal.
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
options(rf.cores=1L)
for(file in c('Config.R','Data_generator.R','ML_functions.R','Score_functions.R',
              'Optimize_functions.R','Reference_functions.R')) source(file)

job_result_path <- function(n,t,method,fixed,seed) {
  if(fixed) {
    folder <- file.path(result_root(),paste0('n_',n),'fixed',paste0('Results_value_',t))
    name <- paste0('Value_',t,'_',seed,'.rds')
  } else {
    parent <- if(method=='IPW_parametric') 'IPW' else method
    folder <- file.path(result_root(),paste0('n_',n),'learned',parent,paste0('Results_',method,'_',t))
    name <- paste0(method,'_',t,'_',seed,'.rds')
  }
  file.path(folder,name)
}

run_job <- function(n,t,method='AIPW',fixed=FALSE,seed) {
  cfg <- get_config()
  stopifnot(n %in% c(1000L,1500L,2500L,5000L),t %in% c(1L,4L),
            method %in% c('AIPW','DS','IPW','IPW_parametric'),
            seed>=1L,seed<=cfg$replications,!fixed || method=='AIPW')
  path <- job_result_path(n,t,method,fixed,seed)
  start <- proc.time()[3]
  warnings <- character()
  meta <- list(n=n,t=t,method=method,fixed=fixed,seed=seed,test_mode=cfg$test_mode)
  result <- tryCatch(withCallingHandlers({
    required <- c('SuperLearner','lightgbm','randomForestSRC','survival','rgenoud')
    missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
    if(length(missing)) stop('Missing R packages: ',paste(missing,collapse=', '))
    ref <- load_reference()
    reference_id <- unname(tools::md5sum(reference_path()))
    result_code <- c('Optimize_functions.R','Reference_functions.R','Run_functions.R')
    result_code_id <- paste(unname(tools::md5sum(result_code)),collapse=':')
    signature <- paste0(cfg$signature,'|result_code=',result_code_id,'|reference=',reference_id)
    if(file.exists(path)) {
      old <- tryCatch(readRDS(path),error=function(e) NULL)
      if(is.list(old) && identical(old$status,'ok') && identical(old$signature,signature)) {
        message('Already complete: ',path)
        return(invisible(old))
      }
    }
    train_constraint <- train_value <- optimal_train_value <- NA_real_
    optimal_train_constraint <- optimal_train_constraint_se <- NA_real_
    optimizer <- NULL; train_diagnostics <- NULL
    theta <- c(.5,-1,1)
    if(!fixed) {
      train_seed <- dataset_seed(n,seed,FALSE)
      train <- simulate_valid_idid(n,seed=train_seed)
      key <- file.path(paste0('n_',n),paste0('training_',seed))
      message('Learning ',method,', n=',n,', t=',t,', replication=',seed)
      scores <- get_scores(train,t,method,train_seed,key,cfg)
      dref <- reference_policy(ref,train$data)
      optimizer <- learn_policy(train$data,scores,dref,train_seed+t*1000L,cfg)
      theta <- optimizer$theta
      train_constraint <- mean((policy(theta,train$data)-dref)*scores$Y)
      training_value <- mean((policy(theta,train$data)-dref)*scores$N)
      train_value <- training_value
      optimal_difference <- policy(c(.5,-1,1),train$data)-dref
      optimal_train_value <- mean(optimal_difference*scores$N)
      optimal_survival_score <- optimal_difference*scores$Y
      optimal_train_constraint <- mean(optimal_survival_score)
      optimal_train_constraint_se <- sqrt(mean((optimal_survival_score-optimal_train_constraint)^2)/n)
      train_diagnostics <- scores$diagnostics
    }
    actual <- true_contrast(theta,t,ref)
    optimal <- true_contrast(c(.5,-1,1),t,ref)
    estimate <- if(fixed) NA_real_ else training_value
    se <- lower <- upper <- evaluation_constraint <- NA_real_
    achieved_covered <- optimal_covered <- NA
    diagnostics <- if(fixed) list() else train_diagnostics
    inference <- method=='AIPW'
    if(inference) {
      eval_seed <- dataset_seed(n,seed,TRUE)
      evaluation <- simulate_valid_idid(n,seed=eval_seed)
      key <- file.path(paste0('n_',n),paste0('evaluation_',seed))
      message('Independent AIPW evaluation, n=',n,', t=',t,', replication=',seed)
      eval_scores <- get_scores(evaluation,t,'AIPW',eval_seed,key,cfg)
      difference <- policy(theta,evaluation$data)-reference_policy(ref,evaluation$data)
      contribution <- difference*eval_scores$N
      estimate <- mean(contribution)
      variance <- mean((contribution-estimate)^2)
      se <- sqrt(variance/n)
      lower <- estimate-qnorm(.975)*se; upper <- estimate+qnorm(.975)*se
      evaluation_constraint <- mean(difference*eval_scores$Y)
      achieved_covered <- lower<=actual[1] && actual[1]<=upper
      optimal_covered <- lower<=optimal[1] && optimal[1]<=upper
      diagnostics <- eval_scores$diagnostics
      if(!fixed) diagnostics$train_delta_clipped <- train_diagnostics$delta_clipped
    }
    diagnostics$warning_count <- length(unique(warnings))
    c(meta,list(status='ok',theta=as.numeric(theta),inference=inference,
      accuracy=policy_accuracy(theta),true_value=unname(actual[1]),
      true_survival=unname(actual[2]),optimal_value=unname(optimal[1]),
      estimated_value=estimate,se=se,ci_lower=lower,ci_upper=upper,
      covered_achieved=achieved_covered,covered_optimal=optimal_covered,
      train_constraint=train_constraint,evaluation_constraint=evaluation_constraint,
      train_value=train_value,optimal_train_value=optimal_train_value,
      optimal_train_constraint=optimal_train_constraint,
      optimal_train_constraint_se=optimal_train_constraint_se,
      elapsed=unname(proc.time()[3]-start),diagnostics=diagnostics,
      warnings=unique(warnings),optimizer=optimizer,signature=signature,
      reference_id=reference_id,configuration=cfg,
      R_version=R.version.string,
      package_versions=vapply(required,function(p) as.character(packageVersion(p)),character(1))))
  },warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')
  }),error=function(e) {
    c(meta,list(status='failed',error=conditionMessage(e),warnings=unique(warnings),
               elapsed=unname(proc.time()[3]-start),signature=cfg$signature))
  })
  atomic_rds(result,path,compress=TRUE)
  if(!identical(result$status,'ok')) stop('Job failed; diagnostic result saved at ',path,': ',result$error)
  message('Saved ',path,'; elapsed ',round(result$elapsed,2),' seconds')
  invisible(result)
}
