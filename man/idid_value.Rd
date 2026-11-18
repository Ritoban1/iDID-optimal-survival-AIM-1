\name{idid_value}
\alias{idid_value}
\title{Evaluate Count and Survival Contrasts Against a Reference Rule}
\description{Compute policy value contrasts from evaluation-sample scores.}
\usage{
idid_value(fit, t, rule, reference, method = "AIPW", conf_level = 0.95)
}
\arguments{
\item{fit}{Evaluation scores from \code{\link{idid_fit_scores}}.}
\item{t,method}{One fitted landmark and method.}
\item{rule,reference}{Each is a function of the covariate data frame, a vector
of 0/1 assignments in row order, or a fitted \code{idid_rule}.}
\item{conf_level}{Confidence level for AIPW Wald intervals.}
}
\details{When evaluating a learned rule, construct the evaluation scores with
data independent of rule training. The contrast compares the supplied rule
with the supplied reference. IPW1 and IPW2 return point estimates.}
\value{A data frame with one row per outcome and the estimated contrast,
standard error, and confidence limits. Standard errors and confidence limits
are available for AIPW.}
