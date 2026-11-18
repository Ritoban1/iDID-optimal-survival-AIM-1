\name{idid_control}
\alias{idid_control}
\title{Configure iDID Nuisance Fitting and Policy Search}
\description{Select SuperLearner libraries and tune cross-fitting, the survival forest, and policy search.}
\usage{
idid_control(profile = c("default", "lightgbm"), sl_library = NULL,
  outer_folds = 5L, inner_folds = 3L, rf_trees = 500L,
  rf_nodesize = 15L, threads = 1L, lgb_rounds = 200L,
  b_points = 21L, prob_floor = .01, delta_floor = .05,
  surv_floor = .02, first_stage_sign = 1L, b_extra_times = numeric(),
  theta_bound = 3, genoud_population = 2000L, genoud_generations = 40L,
  genoud_wait = 6L)
}
\arguments{
\item{profile}{Default uses a mean learner and a GLM with a period-by-instrument term, and defaults to Cox survival.
The \code{"lightgbm"} profile uses a mean learner and two LightGBM learners
and defaults to randomForestSRC.}
\item{sl_library}{A character vector of SuperLearner wrapper names, or a named
list with separate character vectors for \code{binomial} and \code{gaussian}.
The latter library must support continuous weighted outcomes.}
\item{outer_folds,inner_folds}{Number of cross-fitting and SuperLearner folds,
respectively.}
\item{rf_trees,rf_nodesize,threads}{Forest tree count, node size and number of
threads. Threads also governs the bundled LightGBM wrappers.}
\item{lgb_rounds}{Boosting rounds for the bundled LightGBM learners.}
\item{b_points}{Minimum number of equally spaced risk-set regression knots.}
\item{prob_floor,delta_floor,surv_floor}{Probability floors for the period and
instrument models, oriented first stage, and survival curves. The period and
instrument probabilities are also capped at one minus their floor.}
\item{first_stage_sign}{Expected sign of the four-cell first stage, 1 or -1.
Choose it according to the design and instrument coding. The sign is assumed
common across covariate values.}
\item{b_extra_times}{Optional additional risk-set knots beyond the automatically
chosen grid and requested landmarks.}
\item{theta_bound}{Positive bound on the intercept and free coefficients.}
\item{genoud_population,genoud_generations,genoud_wait}{Population size,
generation limit and stopping patience for the policy search.}
}
\value{An object of class \code{idid_control}, for use with
\code{\link{idid_fit_scores}} and \code{\link{idid_learn}}.}
