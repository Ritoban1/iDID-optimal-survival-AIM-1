# ididPolicy

R source package for instrumented difference-in-differences policy learning
with recurrent outcomes and a competing terminal event. It implements IPW1,
IPW2 and two-part AIPW, corresponding to the completed-history, event-integral
and augmented scores in the study. It provides cross-fitted scores, a
survival-constrained linear treatment rule and value contrasts against a
user-supplied reference rule.

## Install locally

Install `SuperLearner` and `survival` in your R library. Install
`randomForestSRC` for the forest option, `lightgbm` for the paper learner
library and `rgenoud` for rule learning. From the repository root, run:

```sh
R CMD INSTALL .
```

## Run the example

With `rgenoud` installed, run the example from the repository root.

```sh
Rscript examples/example.R
```

The example learns rules with all three methods and evaluates them on an
independent sample. See [validation results](validation/Results.txt) for the
comparison with the study's simulation code and
[reproduction instructions](validation/README.txt) to run it again.

## Prepare data

Use one row per subject. Supply binary treatment `A`, period `L`, instrument
`Z`, observed follow-up `X = min(C, D, tau)`, completion indicator
`delta = I(C >= min(D, tau))` and observed-death indicator
`J = I(D <= C, D <= tau)`. At ordinary administrative completion, supply
`X = tau, delta = 1, J = 0`. Map other column names with `columns` (see
`?idid_fit_scores`). Baseline covariates are finite numeric columns. Create
consistent indicator columns for categorical covariates in any training,
evaluation and prediction data. Scale continuous covariates as appropriate
for the bounded linear policy search.

Supply `events` as a list in data row order. Each vector contains recorded
recurrent event times, including events observed before censoring. Use
`numeric(0)` for subjects with no recurrences. Repeat a time for multiple
events at that time. Use the same time units for follow-up, events, and
landmarks. IPW2 computes the exact counting-process sum of the observed
jumps divided by censoring survival at their left limits. Data containing
only cumulative counts at requested landmarks are insufficient for this
calculation unless the jump times themselves are known.

## Fit and evaluate

```r
library(ididPolicy)
control <- idid_control(
  sl_library = list(binomial = c("SL.mean", "SL.idid_cell_glm"),
                    gaussian = c("SL.mean", "SL.idid_cell_glm"))
)
landmarks <- c(0.8, 2.75, 4)
w <- c("age_scaled", "risk_score", "baseline_marker")
training <- idid_fit_scores(
  train, train_events, times = landmarks, tau = 4, covariates = w,
  censor = "randomForestSRC", death = "cox", control = control, seed = 11L
)

# A user-supplied rule, defined before evaluation.
reference <- function(w) as.integer(w$risk_score > 0)
rule <- idid_learn(training, t = 2.75, reference = reference, method = "AIPW",
                  rule_covariates = c("risk_score", "baseline_marker"),
                  normalize_by = "baseline_marker")
predict(rule, new_subjects)

evaluation <- idid_fit_scores(
  test, test_events, times = landmarks, tau = 4, covariates = w,
  methods = "AIPW", censor = "randomForestSRC", death = "cox",
  control = control, seed = 12L
)
idid_value(evaluation, t = 2.75, rule = rule, reference = reference)
```

`train` and `test` in this example are independent subject samples.
Replace the learned rule with a fixed function to evaluate a prespecified rule.
Choose `method = "IPW1"` or `"IPW2"` in `idid_learn()` to learn with either
IPW score. `idid_value()` returns count and survival contrasts. Its AIPW
interval is computed for the supplied rule on the evaluation sample. The
simulation comparison uses the known population optimum as a separate
benchmark.

Choose any distinct landmarks in `(0, tau]`, including fractional times.
Select `censor` and `death` independently. Default control uses Cox and
`SL.mean` plus a GLM that includes the period-by-instrument interaction
when both variables enter the regression. The `paper` profile uses randomForestSRC and
`SL.mean` plus the bundled depth-2 and depth-4 LightGBM wrappers. Provide
wrapper names as a single vector or separate libraries by outcome family.
For example, add `SL.glmnet` after installing its dependency. Pass
`learner_env` for custom wrappers, with any required S3 prediction methods
available as for ordinary SuperLearner. Set `propensity = "logistic"` and
`censor = "cox"` to reproduce the paper's working parametric IPW1 comparator.

AIPW estimates two-part outcome nuisances using binomial completion and
Gaussian weighted-response fits. The marginal models condition on
`(L, Z, W)` and the recurrent-event risk-set models on `(A, L, Z, W)`.
Completed zero counts are retained. Shared nuisances are reused across the
methods and landmarks within a fit call. The paper profile uses out-of-bag
forest predictions for the training responses. Cox uses the fitted
training-fold Cox model and Breslow baseline hazards.

The estimator relies on the study's causal identification and support
conditions. Inspect `fit$scores$AIPW$diagnostics` for clipping, small risk sets,
and the number of folds that required late-time extrapolation. Set the
first-stage orientation consistently with your instrument coding via
`first_stage_sign`. Consider calibrating the probability floors for the
application. Fit and evaluate on independent samples when conducting inference
for a learned rule.

See the R help pages for the data contract and arguments.
