idid_numeric_frame <- function(x, names_required = NULL) {
  x <- as.data.frame(x)
  if (!is.null(names_required)) {
    if (!all(names_required %in% names(x))) stop("Missing prediction covariates.")
    x <- x[, names_required, drop = FALSE]
  }
  if (!ncol(x) || !all(vapply(x, is.numeric, logical(1))))
    stop("Learner covariates must be nonempty numeric columns.")
  if (any(!is.finite(as.matrix(x)))) stop("Nonfinite learner covariates.")
  x
}

idid_lgb_fit <- function(Y, X, newX, family, depth) {
  cfg <- getOption("idid.ml.cfg")
  if (!requireNamespace("lightgbm", quietly = TRUE))
    stop("Install lightgbm to use the built-in LightGBM learners.")
  X <- idid_numeric_frame(X)
  newX <- idid_numeric_frame(newX, names(X))
  objective <- if (family$family == "binomial") "binary" else "regression"
  params <- list(objective = objective,
    metric = if (objective == "binary") "binary_logloss" else "l2",
    max_depth = as.integer(depth), num_leaves = as.integer(2^depth),
    learning_rate = 0.03, min_data_in_leaf = 20L,
    feature_pre_filter = FALSE, verbosity = -1L,
    num_threads = as.integer(cfg$threads), deterministic = TRUE,
    force_col_wise = TRUE, seed = as.integer(cfg$seed + depth),
    data_random_seed = as.integer(cfg$seed + depth),
    feature_fraction_seed = as.integer(cfg$seed + depth),
    bagging_seed = as.integer(cfg$seed + depth))
  dtrain <- lightgbm::lgb.Dataset(as.matrix(X), label = as.numeric(Y),
    params = list(feature_pre_filter = FALSE))
  object <- lightgbm::lgb.train(params = params, data = dtrain,
    nrounds = as.integer(cfg$lgb_rounds), verbose = -1L,
    serializable = TRUE)
  fit <- list(object = object, feature_names = names(X))
  class(fit) <- "SL.idid_lgb"
  list(pred = as.numeric(predict(object, as.matrix(newX))), fit = fit)
}

SL.idid_lgb2 <- function(Y, X, newX, family, ...) {
  idid_lgb_fit(Y, X, newX, family, depth = 2L)
}

SL.idid_lgb4 <- function(Y, X, newX, family, ...) {
  idid_lgb_fit(Y, X, newX, family, depth = 4L)
}

predict.SL.idid_lgb <- function(object, newdata, ...) {
  newdata <- idid_numeric_frame(newdata, object$feature_names)
  as.numeric(predict(object$object, as.matrix(newdata)))
}

# Include the period-by-instrument term in the built-in GLM. The DID first
# stage relies on an interaction that an additive model can miss.
add_cell_interaction <- function(x) {
  if (all(c("L", "Z") %in% names(x))) x$.idid_LZ <- x$L * x$Z
  x
}

SL.idid_cell_glm <- function(Y, X, newX, family, obsWeights, ...) {
  fit <- SuperLearner::SL.glm(Y, add_cell_interaction(X),
    add_cell_interaction(newX), family, obsWeights, model = FALSE)
  class(fit$fit) <- "SL.idid_cell_glm"
  fit
}

predict.SL.idid_cell_glm <- function(object, newdata, ...) {
  as.numeric(predict(object$object, newdata = add_cell_interaction(newdata),
                     type = "response"))
}

fit_ml <- function(y, x, family = c("gaussian", "binomial"), seed, cfg) {
  family <- match.arg(family)
  x <- idid_numeric_frame(x)
  y <- as.numeric(y)
  if (length(y) != nrow(x) || !length(y) || any(!is.finite(y)))
    stop("Invalid learner response.")
  if (family == "binomial" && any(!y %in% c(0, 1)))
    stop("Binomial learner response must contain only zero and one.")
  constant_reason <- if (length(unique(y)) == 1L) "constant_response" else
    if (length(y) < 10L) "fewer_than_10_rows" else NULL
  if (family == "binomial" && is.null(constant_reason) &&
      min(tabulate(as.integer(y) + 1L, nbins = 2L)) < 2L)
    constant_reason <- "fewer_than_2_in_binary_class"
  if (!is.null(constant_reason)) return(list(type = "constant",
    value = mean(y), feature_names = names(x), family = family,
    diagnostics = list(n = length(y), constant_reason = constant_reason)))
  V <- min(as.integer(cfg$inner_folds), nrow(x))
  if (family == "binomial") V <- min(V,
    min(tabulate(as.integer(y) + 1L, nbins = 2L)))
  if (V < 2L) stop("At least two inner folds are required.")
  cfg$seed <- as.integer(seed)
  old <- options(idid.ml.cfg = cfg)
  on.exit(options(old), add = TRUE)
  set.seed(seed)
  # Built-in learners and screening methods are available without attaching
  # SuperLearner. The supplied environment can override them by name.
  library <- cfg$sl_library[[family]]
  env <- new.env(parent = asNamespace("SuperLearner"))
  wrapper_names <- unique(as.character(unlist(library, use.names = FALSE)))
  for (name in wrapper_names) {
    if (name %in% c("SL.idid_lgb2", "SL.idid_lgb4", "SL.idid_cell_glm")) {
      assign(name, get(name, envir = environment(fit_ml)), envir = env)
    } else if (exists(name, envir = cfg$learner_env, mode = "function",
                      inherits = TRUE)) {
      assign(name, get(name, envir = cfg$learner_env, mode = "function"), envir = env)
    }
  }
  fit <- SuperLearner::SuperLearner(Y = y, X = x,
    family = if (family == "binomial") stats::binomial() else stats::gaussian(),
    SL.library = library, method = "method.NNLS",
    cvControl = list(V = V, stratifyCV = identical(family, "binomial")),
    control = list(saveFitLibrary = TRUE), verbose = FALSE, env = env)
  if (any(fit$errorsInLibrary) || any(fit$errorsInCVLibrary))
    stop("A SuperLearner library fit failed; inspect the recorded warnings.")
  # Prediction needs the fitted library and ensemble weights, not the large
  # training/CV arrays or the enclosing fit environment. This matters when the
  # frozen reference rule is trained on one million observations.
  fit[c("env", "SL.predict", "library.predict", "Z", "validRows",
    "cvFitLibrary", "metaOptimizer")] <- NULL
  list(type = "sl", object = fit, feature_names = names(x), family = family,
    diagnostics = list(n = length(y), inner_folds = V,
      weights = fit$coef, constant_reason = NULL))
}

predict_ml <- function(fit, newdata) {
  newdata <- idid_numeric_frame(newdata, fit$feature_names)
  if (fit$type == "constant") ans <- rep(fit$value, nrow(newdata)) else
    if (fit$type == "sl") ans <- as.numeric(
      predict(fit$object, newdata = newdata, onlySL = TRUE)$pred) else
        stop("Unrecognized learner object.")
  if (length(ans) != nrow(newdata) || any(!is.finite(ans)))
    stop("Nonfinite or malformed learner predictions.")
  # In particular, do not trim treatment propensities at 0.01: some true
  # treatment probabilities are smaller, while the DID first stage stays strong.
  if (fit$family == "binomial") ans <- pmin(1, pmax(0, ans))
  ans
}

# Step-curve lookup. Before the first jump the survival probability is one.
# left = TRUE uses times strictly less than q, including when q is a jump time.
curve_lookup <- function(s, times, q, left = TRUE) {
  s <- as.matrix(s)
  times <- as.numeric(times)
  if (ncol(s) != length(times) || is.unsorted(times, strictly = TRUE))
    stop("Malformed survival curve grid.")
  if (length(q) == 1L) q <- rep(q, nrow(s))
  if (length(q) != nrow(s) || anyNA(q)) stop("Malformed survival query.")
  ix <- findInterval(q, times, left.open = left)
  ans <- rep(1, nrow(s))
  has_jump <- ix > 0L
  ans[has_jump] <- s[cbind(which(has_jump), ix[has_jump])]
  ans
}

curve_eval <- function(s, times, grid, left = FALSE) {
  s <- as.matrix(s)
  times <- as.numeric(times)
  if (ncol(s) != length(times) || is.unsorted(times, strictly = TRUE) ||
      anyNA(grid)) stop("Malformed survival curve grid.")
  ix <- findInterval(grid, times, left.open = left)
  ans <- matrix(1, nrow = nrow(s), ncol = length(grid))
  has_jump <- ix > 0L
  ans[, has_jump] <- s[, ix[has_jump], drop = FALSE]
  ans
}

idid_check_survival_data <- function(train, newdata, cfg) {
  required <- cfg$qcols
  x <- idid_numeric_frame(train, required)
  nx <- idid_numeric_frame(newdata, required)
  if (!all(c("X", "delta", "J") %in% names(train)) ||
      any(!is.finite(train$X)) || any(train$X < 0) ||
      any(!train$delta %in% c(0, 1)) || any(!train$J %in% c(0, 1)))
    stop("Invalid survival training outcomes.")
  list(x = x, nx = nx)
}

fit_survival <- function(train, newdata, event = c("censor", "death"),
                         seed, cfg) {
  event <- match.arg(event)
  checked <- idid_check_survival_data(train, newdata, cfg)
  model <- if (event == "censor") cfg$censor else cfg$death
  if (model == "cox") return(fit_cox_survival(train, newdata, event, cfg))
  if (!requireNamespace("randomForestSRC", quietly = TRUE))
    stop("Install randomForestSRC to use the forest survival option.")
  previous_cores <- options(rf.cores=as.integer(cfg$threads))
  on.exit(options(previous_cores),add=TRUE)
  status <- if (event == "censor") 1 - train$delta else train$J
  if (!sum(status)) return(list(times = numeric(),
    test = matrix(numeric(), nrow(newdata), 0L),
    train_oob_at_X = rep(1, nrow(train)),
    diagnostics = list(type = "rf", event = event, events = 0L,
      test_floored = 0L, train_oob_floored = 0L, grid_points = 0L)))
  data <- data.frame(X = train$X, status = status, checked$x)
  invalid_endpoint <- if(event=='censor') train$X[status==1L]>=cfg$tau else train$X[status==1L]>cfg$tau
  if (any(invalid_endpoint))
    stop("An administrative endpoint was incorrectly marked as an event.")
  # RF-SRC selects this many actual event times. Administrative completion
  # carries status zero and therefore creates no artificial terminal jump.
  grid_n <- max(150L, ceiling(4 * sqrt(nrow(train))))
  Surv <- survival::Surv
  formula <- stats::reformulate(cfg$qcols, response = "Surv(X, status)",
    env = environment())
  set.seed(seed)
  # Growing subsampling follows the RF-SRC fast-forest approach, with a more
  # conservative minimum for this study. OOB predictions remain available.
  sample_n <- max(1L,floor(min(nrow(train)*.632,
    max(cfg$rf_sample_min,nrow(train)^cfg$rf_sample_power))))
  object <- randomForestSRC::rfsrc(formula, data = data,
    ntree = as.integer(cfg$rf_trees), nodesize = as.integer(cfg$rf_nodesize),
    splitrule = "logrank", samptype = "swor", sampsize = sample_n,
    ntime = as.integer(grid_n),
    perf.type = "none", importance = "none",
    seed = -as.integer(abs(seed) %% (.Machine$integer.max - 1L) + 1L))
  pred <- predict(object, newdata = checked$nx)
  times <- as.numeric(object$time.interest)
  test_raw <- as.matrix(pred$survival)
  oob_raw <- as.matrix(object$survival.oob)
  if (length(times) != ncol(test_raw) || length(times) != ncol(oob_raw) ||
      nrow(test_raw) != nrow(newdata) || nrow(oob_raw) != nrow(train) ||
      any(!is.finite(test_raw)) || any(!is.finite(oob_raw)))
    stop("Nonfinite or malformed RF survival predictions; no in-bag fallback is used.")
  oob_at_X <- curve_lookup(oob_raw, times, train$X, left = TRUE)
  test <- pmax(pmin(test_raw, 1), cfg$surv_floor)
  list(times = times, test = test,
    train_oob_at_X = pmax(cfg$surv_floor, pmin(1, oob_at_X)),
    diagnostics = list(type = "rf", event = event, events = sum(status),
      grid_points = length(times), subsample_size = sample_n,
      test_floored = sum(test_raw < cfg$surv_floor),
      train_oob_floored = sum(oob_at_X < cfg$surv_floor),
      min_test_raw = min(test_raw), min_train_oob_raw = min(oob_at_X)))
}

# Cox fit with Breslow jumps, using the same curve interface as the forest.
fit_cox_survival <- function(train, newdata, event, cfg) {
  checked <- idid_check_survival_data(train, newdata, cfg)
  status <- if (event == "censor") 1 - train$delta else train$J
  if (!sum(status)) return(list(times = numeric(),
    test = matrix(numeric(), nrow(newdata), 0L), train_oob_at_X = rep(1, nrow(train)),
    diagnostics = list(type = "cox", event = event, events = 0L,
      test_floored = 0L, train_oob_floored = 0L, grid_points = 0L)))
  data <- data.frame(X = train$X, status = status, checked$x)
  Surv <- survival::Surv
  formula <- stats::reformulate(cfg$qcols, response = "Surv(X, status)",
    env = environment())
  object <- survival::coxph(formula, data = data, ties = "breslow",
    singular.ok = TRUE)
  if (any(is.infinite(stats::coef(object)))) stop("Nonfinite Cox coefficients.")
  bh <- survival::basehaz(object, centered = FALSE)
  jumps <- c(bh$hazard[1L], diff(bh$hazard)) > 0
  times <- bh$time[jumps]
  hazard <- bh$hazard[jumps]
  lp_test <- as.numeric(predict(object, newdata = checked$nx,
    type = "lp", reference = "zero"))
  lp_train <- as.numeric(predict(object, newdata = checked$x,
    type = "lp", reference = "zero"))
  test_raw <- exp(-outer(exp(lp_test), hazard))
  ix <- findInterval(train$X, times, left.open = TRUE)
  at_X <- rep(0, nrow(train))
  has_jump <- ix > 0L
  at_X[has_jump] <- hazard[ix[has_jump]]
  train_raw <- exp(-exp(lp_train) * at_X)
  if (any(!is.finite(test_raw)) || any(!is.finite(train_raw)))
    stop("Nonfinite Cox predictions.")
  list(times = times, test = pmax(pmin(test_raw, 1), cfg$surv_floor),
    # The name matches the RF interface; these are fitted Cox probabilities,
    # not OOB predictions. The enclosing outer split still excludes test rows.
    train_oob_at_X = pmax(cfg$surv_floor, pmin(1, train_raw)),
    diagnostics = list(type = "cox", event = event, events = sum(status),
      grid_points = length(times), test_floored = sum(test_raw < cfg$surv_floor),
      train_oob_floored = sum(train_raw < cfg$surv_floor),
      min_test_raw = min(test_raw), min_train_oob_raw = min(train_raw)))
}
