library(ididPolicy)
# At u=1, the left limit is 1. At u=2.3 the left limit is .8.
# Events after the requested landmark are excluded, including same-subject events.
histories <- list(c(.5, 1, 2.3, 2.3, 3.2), numeric())
curve <- list(times = c(1, 2.9), test = rbind(c(.8, .6), c(.9, .7)))
expected <- c(2 + 2/.8, 0)
stopifnot(identical(ididPolicy:::weighted_events(histories, curve, 2.5, NULL), expected),
          identical(ididPolicy:::weighted_events(histories, curve, 2.5,
                                                     c(.5, 1, 2.3, 3.2)), expected))
# Three baseline covariates, renamed columns, noninteger landmarks and custom
# wrappers in an unattached environment. A negative first-stage orientation
# should be preserved when requested.
set.seed(102)
n <- 240L
tau <- 4.25
s <- data.frame(w1 = runif(n, -1, 1), w2 = rnorm(n), w3 = rnorm(n),
                period = rep(c(0L, 1L), each = n/2), iv = rep(c(0L, 1L), n/2))
s$trt <- rbinom(n, 1, ifelse(s$period == 0L, .6, ifelse(s$iv == 1L, .15, .85)))
end <- rexp(n, .25)
censor <- rexp(n, .16)
s$stop <- pmin(end, censor, tau)
s$dead <- as.integer(end <= pmin(censor, tau))
s$complete <- as.integer(censor >= pmin(end, tau))
events <- lapply(s$stop, function(x) runif(rpois(1, x), max = x))
fields <- c(A = "trt", L = "period", Z = "iv", X = "stop", delta = "complete", J = "dead")
covars <- c("w1", "w2", "w3")
env <- new.env(parent = baseenv())
env$SL.my_mean <- SuperLearner::SL.mean
control <- idid_control(sl_library = list(binomial = "SL.mean", gaussian = "SL.my_mean"),
                        first_stage_sign = -1, b_points = 5L,
                        outer_folds = 3L, genoud_population = 60L,
                        genoud_generations = 5L, genoud_wait = 2L)
times <- c(1.37, tau, .48)
fit <- idid_fit_scores(s, events, times, tau, covars, columns = fields,
                       control = control, learner_env = env)
stopifnot(identical(fit$times, times), all(fit$first_stage <= -.05),
          identical(dim(fit$scores$IPW2$N), c(n, 3L)),
          length(unique(fit$fold_id)) == 3L)
for (element in fit$scores) {
  stopifnot(all(is.finite(element$N)), all(is.finite(element$Y)))
}
reference <- function(w) as.integer(w$w1 > 0)
fixed <- function(w) as.integer(w$w3 > 0)
value <- idid_value(fit, times[1L], fixed, reference)
brute <- (fixed(s) - reference(s)) * fit$scores$AIPW$N[, 1]
stopifnot(all.equal(value$estimate[1], mean(brute)),
          all.equal(value$se[1], sqrt(mean((brute - mean(brute))^2)/n)),
          is.na(idid_value(fit, times[1L], fixed, reference, "IPW2")$se[1]))
invalid <- try(idid_fit_scores(s, events, tau + .01, tau, covars, columns = fields), silent = TRUE)
stopifnot(inherits(invalid, "try-error"))
# A no-censoring fold produces an empty jump grid for K.
nc <- s
nc$stop <- pmin(end, tau)
nc$dead <- as.integer(end <= tau)
nc$complete <- 1L
fit_nc <- idid_fit_scores(nc, events, times[1L], tau, covars, methods = "AIPW",
                          columns = fields, control = control, learner_env = env)
stopifnot(all(is.finite(fit_nc$scores$AIPW$N)))
# Verify that both survival engines can be selected independently.
if (requireNamespace("randomForestSRC", quietly = TRUE)) {
  ctl <- idid_control(sl_library = "SL.mean", outer_folds = 3L,
                      rf_trees = 60L, b_points = 5L, first_stage_sign = -1)
  mixed <- idid_fit_scores(s, events, times[1L], tau, covars, methods = "AIPW",
                           columns = fields, censor = "cox",
                           death = "randomForestSRC", control = ctl)
  stopifnot(identical(unname(mixed$survival_methods), c("cox", "randomForestSRC")),
            all(is.finite(mixed$scores$AIPW$Y)))
}
if (requireNamespace("rgenoud", quietly = TRUE)) {
  # With a zero survival score, all candidate rules are feasible in this API check.
  fit$scores$AIPW$Y[] <- 0
  for (cols in list(covars, "w1")) {
    rule <- idid_learn(fit, times[1L], reference, rule_covariates = cols)
    stopifnot(length(rule$coefficients) == length(cols) + 1L,
              abs(tail(rule$coefficients, 1L)) == 1,
              length(predict(rule, s)) == n)
  }
}
cat("PASS: fractional landmarks, exact IPW2 integral, custom learners, multiple covariates, Cox/RF, and rule/value API.\n")
