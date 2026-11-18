# Install ididPolicy, randomForestSRC, lightgbm, and rgenoud before running.
library(ididPolicy)

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here <- if (length(script)) dirname(normalizePath(sub("^--file=", "", script[1L]))) else "examples"
source(file.path(here, "simulate_data.R"))

training <- simulate_valid_idid(n = 1000L, seed = 2000001L)
evaluation <- simulate_valid_idid(n = 1000L, seed = 101000001L)
landmarks <- c(1, 4)
covariates <- c("V1", "V2")

observed_event_times <- function(sim) {
  previous <- cbind(0, sim$obs_N[, -ncol(sim$obs_N), drop = FALSE])
  increments <- sim$obs_N - previous
  lapply(seq_len(nrow(increments)), function(i) rep(sim$landmark_times, increments[i, ]))
}

control <- idid_control("lightgbm", b_extra_times = training$landmark_times)
training_fit <- idid_fit_scores(
  training$data, observed_event_times(training), times = landmarks, tau = training$tau,
  covariates = covariates, control = control, seed = 2000001L
)
evaluation_fit <- idid_fit_scores(
  evaluation$data, observed_event_times(evaluation), times = landmarks, tau = evaluation$tau,
  covariates = covariates, methods = "AIPW", control = control, seed = 101000001L
)

reference <- function(w) rep(0L, nrow(w))
for (t in landmarks) {
  for (method in names(training_fit$scores)) {
    rule <- idid_learn(training_fit, t, reference, method = method, seed = 2000001L + 1000L * t)
    cat("\nLandmark", t, "| method", method, "\n")
    print(rule$coefficients)
    if (method == "AIPW") print(idid_value(evaluation_fit, t, rule, reference))
  }
}
