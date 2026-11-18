\name{idid_fit_scores}
\alias{idid_fit_scores}
\title{Compute Cross-Fitted Scores at User-Specified Landmarks}
\description{Fit the nuisance functions by cross-fitting and compute count and survival scores for IPW1, IPW2, and two-part AIPW.}
\usage{
idid_fit_scores(data, events, times, tau, covariates,
  methods = c("IPW1", "IPW2", "AIPW"),
  columns = c(A = "A", L = "L", Z = "Z", X = "X", delta = "delta", J = "J"),
  censor = NULL, death = NULL, propensity = c("superlearner", "logistic"),
  control = idid_control(), seed = 1L, fold_id = NULL,
  learner_env = parent.frame(), verbose = FALSE)
}
\arguments{
\item{data}{Data frame with one row per subject.}
\item{events}{List of vectors of observed recurrent-event times, in data row
order. Use \code{numeric(0)} for no observed events. A repeated time represents
multiple events at that time.}
\item{times}{Distinct numeric landmarks in \eqn{(0,\tau]}, in the desired
output order. Fractional times are allowed.}
\item{tau}{Administrative end of follow-up, in the same units as \code{times}
and \code{events}.}
\item{covariates}{Names of finite numeric baseline predictors. Encode categorical
predictors as numeric indicators before calling this function.}
\item{methods}{Any subset of IPW1, IPW2 and AIPW.}
\item{columns}{Named character mapping from internal field names to data columns.
\code{A}, \code{L}, and \code{Z} are binary treatment, period and instrument.
\code{X} is observed follow-up, \code{delta} indicates observation through death
or \code{tau}, and \code{J} indicates an observed death.}
\item{censor,death}{Independently select \code{"cox"} or
\code{"randomForestSRC"}. NULL selects the profile default. Death is fitted
only for AIPW.}
\item{propensity}{SuperLearner or logistic period, instrument, and treatment
regressions. Logistic treatment fits are by period-instrument cell.}
\item{control}{Settings from \code{\link{idid_control}}.}
\item{seed}{Nonnegative integer seed.}
\item{fold_id}{Optional labels for a prespecified cross-fit split, one per row.}
\item{learner_env}{Environment in which to look for custom SuperLearner wrappers
and for overrides of built-in wrappers. Custom prediction methods must be
registered or accessible to R's S3 dispatch.}
\item{verbose}{Print fitting progress.}
}
\details{
Event times must be greater than zero and at most \code{X}. IPW2 sums the
recurrence jumps divided by censoring survival just before each jump. IPW1
uses the completed outcome and the censoring survival just before observed
follow-up ends. AIPW uses separate completion and weighted outcome models
for the marginal means and the recurrent-event risk-set term.

Within a fit call, the methods share the outer folds and the required
nuisances that are common across methods and times. If a training fold has no
subjects observed to the administrative horizon, its risk-set regression
ends at its largest observed follow-up and is carried forward beyond that
point. Diagnostics record the number of such folds. Choose landmarks with
adequate follow-up and inspect the amount of clipping and extrapolation.
}
\value{An \code{idid_scores} list with one \code{scores} element per method.
Each method contains \code{N} and \code{Y} subject-by-landmark score matrices
and per-landmark \code{diagnostics}. The object also returns the landmark
vector, covariates, fold assignments, fitted first stage, fitted joint
period-instrument probability, and selected survival methods.}
