# Run from any directory after installing ididPolicy and optional paper dependencies.
# Rscript validation/compare_legacy.R [sample_size] [replication]
library(ididPolicy)
args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args) >= 1L) as.integer(args[1L]) else 1000L
replication <- if (length(args) >= 2L) as.integer(args[2L]) else 1L
source_arg <- grep("^--file=", commandArgs(), value = TRUE)
home <- dirname(normalizePath(sub("^--file=", "", source_arg[1L])))
output <- file.path(home, sprintf("comparison_n%d_rep%d.csv", n, replication))
tmp <- tempfile("idid-legacy-")
dir.create(tmp)
old_wd <- getwd()
main <- function() {
  on.exit({setwd(old_wd); unlink(tmp, recursive = TRUE)}, add = TRUE)
  Sys.setenv(IDID_FAST = "0", IDID_CACHE_ROOT = file.path(tmp, "cache"),
             IDID_RESULT_ROOT = file.path(tmp, "results"))
  setwd(file.path(home, "legacy"))
  source("Run_functions.R", local = .GlobalEnv)
  cfg <- get_config()
  reference <- load_reference()
  ref <- function(w) reference_policy(reference, w)
  opt <- function(w) as.integer(.5 - w$V1 + w$V2 > 0)
  times <- c(1, 4)
  labels <- c(IPW = "IPW1", DS = "IPW2", AIPW = "AIPW", IPW_parametric = "IPW1")
  make_events <- function(sim) {
    increments <- sim$obs_N - cbind(0, sim$obs_N[, -ncol(sim$obs_N), drop = FALSE])
    lapply(seq_len(nrow(increments)), function(i) rep(sim$landmark_times, increments[i, ]))
  }
  train_seed <- dataset_seed(n, replication, FALSE)
  eval_seed <- dataset_seed(n, replication, TRUE)
  train <- simulate_valid_idid(n, train_seed)
  eval <- simulate_valid_idid(n, eval_seed)
  control <- idid_control("paper", b_extra_times = train$landmark_times)
  train_events <- make_events(train)
  eval_events <- make_events(eval)
  message("Refitting package scores for training and independent evaluation")
  train_fit <- idid_fit_scores(train$data, train_events, times, train$tau,
                              c("V1", "V2"), control = control, seed = train_seed)
  eval_fit <- idid_fit_scores(eval$data, eval_events, times, eval$tau,
                             c("V1", "V2"), methods = "AIPW",
                             control = control, seed = eval_seed)
  par_fit <- idid_fit_scores(train$data, train_events, times, train$tau,
                            c("V1", "V2"), methods = "IPW1", propensity = "logistic",
                            censor = "cox", control = control, seed = train_seed)
  key <- file.path(paste0("n_", n), paste0("training_", replication))
  eval_key <- file.path(paste0("n_", n), paste0("evaluation_", replication))
  rows <- list()
  max_error <- function(a, b) max(abs(a - b))
  add <- function(method, t, arm, diffs) {
    rows[[length(rows) + 1L]] <<- data.frame(n = n, replication = replication,
      method = method, time = t, arm = arm, metric = names(diffs), max_abs_diff = unname(diffs))
  }
  for (j in seq_along(times)) {
    t <- times[j]
    for (method in names(labels)) {
      m <- unname(labels[method])
      pf <- if (method == "IPW_parametric") par_fit else train_fit
      legacy_scores <- get_scores(train, t, method, train_seed, key, cfg)
      stopifnot(identical(as.integer(pf$fold_id), as.integer(legacy_scores$folds)))
      legacy_result <- run_job(n, t, method, FALSE, replication)
      learned <- idid_learn(pf, t, ref, m, seed = train_seed + t * 1000L)
      decisions <- predict(learned, train$data)
      pkg_pcd <- policy_accuracy(unname(learned$coefficients))
      pkg_truth <- true_contrast(unname(learned$coefficients), t, reference)
      train_value <- idid_value(pf, t, learned, ref, m)$estimate[1]
      diffs <- c(count_scores = max_error(pf$scores[[m]]$N[, j], legacy_scores$N),
        survival_scores = max_error(pf$scores[[m]]$Y[, j], legacy_scores$Y),
        coefficients = max_error(learned$coefficients, legacy_result$theta),
        decisions = max_error(decisions, policy(legacy_result$theta, train$data)),
        constraint = max_error(learned$survival_advantage, legacy_result$train_constraint),
        PCD = max_error(pkg_pcd, legacy_result$accuracy),
        achieved_truth = max_error(pkg_truth[1], legacy_result$true_value),
        training_value = max_error(train_value, legacy_result$train_value))
      if (m == "AIPW") {
        value <- idid_value(eval_fit, t, learned, ref)
        diffs <- c(diffs, evaluation_value = max_error(value$estimate[1], legacy_result$estimated_value),
          SE = max_error(value$se[1], legacy_result$se),
          CI = max_error(c(value$lower[1], value$upper[1]),
                         c(legacy_result$ci_lower, legacy_result$ci_upper)),
          optimal_coverage_flag = max_error(
            value$lower[1] <= legacy_result$optimal_value && legacy_result$optimal_value <= value$upper[1],
            legacy_result$covered_optimal))
      }
      add(method, t, "learned", diffs)
    }
    legacy_scores <- get_scores(eval, t, "AIPW", eval_seed, eval_key, cfg)
    legacy_result <- run_job(n, t, "AIPW", TRUE, replication)
    fixed <- idid_value(eval_fit, t, opt, ref)
    add("AIPW", t, "fixed", c(count_scores = max_error(eval_fit$scores$AIPW$N[, j], legacy_scores$N),
      survival_scores = max_error(eval_fit$scores$AIPW$Y[, j], legacy_scores$Y),
      evaluation_value = max_error(fixed$estimate[1], legacy_result$estimated_value),
      SE = max_error(fixed$se[1], legacy_result$se),
      CI = max_error(c(fixed$lower[1], fixed$upper[1]), c(legacy_result$ci_lower, legacy_result$ci_upper)),
      optimal_coverage_flag = max_error(
        fixed$lower[1] <= legacy_result$optimal_value && legacy_result$optimal_value <= fixed$upper[1],
        legacy_result$covered_optimal)))
  }
  report <- do.call(rbind, rows)
  write.csv(report, output, row.names = FALSE)
  print(aggregate(max_abs_diff ~ metric, report, max), row.names = FALSE)
  if (any(!is.finite(report$max_abs_diff)) || any(report$max_abs_diff > 1e-10))
    stop("Differences exceeded tolerance. See ", output)
  cat("PASS:", nrow(report), "comparisons, absolute tolerance 1e-10.\n")
}
main()
