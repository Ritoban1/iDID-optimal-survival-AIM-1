\name{idid_learn}
\alias{idid_learn}
\alias{predict.idid_rule}
\title{Learn a Linear Rule Under a Survival Constraint}
\description{Minimize the recurrent-event objective over a linear rule class subject to an estimated survival constraint.}
\usage{
idid_learn(fit, t, reference, method = "AIPW",
  rule_covariates = names(fit$covariates),
  normalize_by = tail(rule_covariates, 1L), seed = 1L,
  control = fit$control)

\method{predict}{idid_rule}(object, newdata, ...)
}
\arguments{
\item{fit}{An object returned by \code{\link{idid_fit_scores}}.}
\item{t,method}{One fitted landmark and method.}
\item{reference}{A function of the covariate data frame, a vector of 0/1
assignments in row order, or a fitted \code{idid_rule}.}
\item{rule_covariates}{Covariate names used for the decision rule.
Nuisance fitting may use additional covariates.}
\item{normalize_by}{Continuous covariate whose coefficient is set to 1 or -1,
with both signs searched.}
\item{seed}{Nonnegative integer seed for search.}
\item{control}{Settings from \code{\link{idid_control}}.}
\item{object}{A fitted \code{idid_rule}.}
\item{newdata}{Data frame containing the covariates used by the rule.}
\item{...}{Reserved for S3 method dispatch.}
}
\details{A genetic search minimizes the estimated count objective over linear
rules whose estimated survival contrast against the reference is nonnegative.
An error is raised if no feasible rule is found. The search requires rgenoud.}
\value{\code{idid_learn} returns the coefficients, estimated training objective
and survival advantage, time and method. Prediction returns an integer vector
of treatment assignments.}
