# Run after installing ididPolicy. The independent samples below illustrate the API.
library(ididPolicy)
set.seed(119)
tau <- 4.5
times <- c(.75, 2.35, tau)
simulate_example <- function(n) {
  x <- data.frame(w1 = runif(n, -1, 1), w2 = runif(n, -1, 1), w3 = runif(n, -1, 1),
                  L = rbinom(n, 1, .5), Z = rbinom(n, 1, .5))
  probability <- with(x, plogis(ifelse(L == 0, 2.4 - 4.8 * Z, -2.4 + 4.8 * Z) + .25 * w1))
  x$A <- rbinom(n, 1, probability)
  rate <- with(x, exp(-2.3 + .1 * w2 - .2 * A))
  death <- rexp(n, rate)
  censor <- rexp(n, .07)
  x$X <- pmin(death, censor, tau)
  x$delta <- as.integer(censor >= pmin(death, tau))
  x$J <- as.integer(death <= pmin(censor, tau))
  intensity <- with(x, exp(-.15 + .2 * w3 + .5 * A * w1))
  number <- rpois(n, x$X * intensity)
  events <- Map(function(m, u) sort(runif(m, max = u)), number, x$X)
  list(data = x, events = events)
}

training <- simulate_example(600)
evaluation <- simulate_example(600)
covariates <- c("w1", "w2", "w3")
control <- idid_control(
  sl_library = list(binomial = c("SL.mean", "SL.idid_cell_glm"),
                    gaussian = c("SL.mean", "SL.idid_cell_glm")))
reference <- function(w) rep(0L, nrow(w))
training_fit <- idid_fit_scores(training$data, training$events, times, tau, covariates,
                                censor = "cox", death = "cox", control = control)
rules <- lapply(names(training_fit$scores), function(method) {
  idid_learn(training_fit, t = 2.35, reference = reference, method = method)
})
names(rules) <- names(training_fit$scores)
evaluation_fit <- idid_fit_scores(evaluation$data, evaluation$events, times, tau, covariates,
                                  methods = "AIPW", control = control, seed = 2L)
for (method in names(rules)) {
  cat("\nRule learned with", method, "\n")
  print(rules[[method]]$coefficients)
  print(idid_value(evaluation_fit, t = 2.35, rule = rules[[method]], reference = reference))
}
