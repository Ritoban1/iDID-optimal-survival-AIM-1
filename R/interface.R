whole_number <- function(x, min, max) {
  is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x >= min && x <= max && x == floor(x)
}

idid_control <- function(profile = c("default", "lightgbm"), sl_library = NULL,
                         outer_folds = 5L, inner_folds = 3L,
                         rf_trees = 500L, rf_nodesize = 15L, threads = 1L,
                         lgb_rounds = 200L, b_points = 21L,
                         prob_floor = .01, delta_floor = .05, surv_floor = .02,
                         first_stage_sign = 1L, b_extra_times = numeric(),
                         theta_bound = 3, genoud_population = 2000L,
                         genoud_generations = 40L, genoud_wait = 6L) {
  profile <- if (missing(profile)) "default" else
    match.arg(profile, c("default", "lightgbm", "paper"))
  if (profile == "paper") profile <- "lightgbm"
  if (is.null(sl_library)) {
    library <- if (profile == "lightgbm") c("SL.mean", "SL.idid_lgb2", "SL.idid_lgb4")
               else c("SL.mean", "SL.idid_cell_glm")
    sl_library <- list(binomial = library, gaussian = library)
  } else if (is.character(sl_library)) {
    sl_library <- list(binomial = sl_library, gaussian = sl_library)
  }
  valid_library <- function(x) {
    is.character(x) && length(x) > 0L && !anyNA(x) && all(nzchar(x))
  }
  if (!is.list(sl_library) || !all(c("binomial", "gaussian") %in% names(sl_library)) ||
      !all(vapply(sl_library[c("binomial", "gaussian")], valid_library, logical(1)))) {
    stop("sl_library must be a character vector or a list of binomial and gaussian character vectors.")
  }
  if (!whole_number(outer_folds, 2, 100) || !whole_number(inner_folds, 2, 100) ||
      !all(vapply(list(rf_trees, rf_nodesize, threads, lgb_rounds, b_points,
                       genoud_population, genoud_generations, genoud_wait),
                 whole_number, logical(1), min = 1, max = 100000L)) ||
      b_points < 2L || genoud_population < 10L) {
    stop("Invalid fold, learner, knot, or optimizer settings.")
  }
  positive <- function(x, limit) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0 && x < limit
  }
  if (!positive(prob_floor, .5) || !positive(delta_floor, 2) ||
      !positive(surv_floor, 1) || !positive(theta_bound, Inf) ||
      length(first_stage_sign) != 1L || !first_stage_sign %in% c(-1, 1) ||
      !is.numeric(b_extra_times) || any(!is.finite(b_extra_times) | b_extra_times <= 0)) {
    stop("Invalid probability floors, first-stage sign, extra times, or coefficient bound.")
  }
  structure(list(profile = profile, sl_library = sl_library,
    outer_folds = as.integer(outer_folds), inner_folds = as.integer(inner_folds),
    rf_trees = rf_trees, rf_nodesize = rf_nodesize, threads = threads,
    rf_sample_min = 1000L, rf_sample_power = .75,
    lgb_rounds = lgb_rounds, b_points = b_points, prob_floor = prob_floor,
    delta_floor = delta_floor, surv_floor = surv_floor,
    first_stage_sign = first_stage_sign, b_extra_times = b_extra_times,
    theta_bound = theta_bound, genoud_population = genoud_population,
    genoud_generations = genoud_generations, genoud_wait = genoud_wait),
    class = "idid_control")
}

cached <- function(key, cfg, compute) {
  if (!exists(key, envir = cfg$cache, inherits = FALSE)) {
    assign(key, compute(), envir = cfg$cache)
  }
  get(key, envir = cfg$cache, inherits = FALSE)
}

weighted_events <- function(events, curve, t, grid) {
  if (!is.null(grid)) {
    times <- grid[grid <= t]
    if (!length(times)) return(numeric(length(events)))
    increments <- vapply(events, function(x) tabulate(match(x, times), nbins = length(times)),
                         integer(length(times)))
    increments <- matrix(increments, nrow = length(times))
    return(rowSums(t(increments) / curve_eval(curve$test, curve$times, times, left = TRUE)))
  }
  vapply(seq_along(events), function(i) {
    u <- events[[i]]
    u <- u[u <= t]
    ix <- findInterval(u, curve$times, left.open = TRUE)
    k <- c(1, curve$test[i, ])[ix + 1L]
    sum(1 / k)
  }, numeric(1))
}

idid_fit_scores <- function(data, events, times, tau, covariates,
                            methods = c("IPW1", "IPW2", "AIPW"),
                            columns = c(A = "A", L = "L", Z = "Z", X = "X",
                                        delta = "delta", J = "J"),
                            censor = NULL, death = NULL,
                            propensity = c("superlearner", "logistic"),
                            control = idid_control(), seed = 1L,
                            fold_id = NULL, learner_env = parent.frame(),
                            verbose = FALSE) {
  if (!inherits(control, "idid_control")) stop("Use idid_control() to create control.")
  cfg <- control
  default_survival <- if (cfg$profile %in% c("lightgbm", "paper")) "randomForestSRC" else "cox"
  cfg$censor <- match.arg(if (is.null(censor)) default_survival else censor,
                          c("cox", "randomForestSRC"))
  cfg$death <- match.arg(if (is.null(death)) default_survival else death,
                         c("cox", "randomForestSRC"))
  cfg$propensity <- match.arg(propensity)
  methods <- match.arg(methods, c("IPW1", "IPW2", "AIPW"), several.ok = TRUE)
  if (anyDuplicated(methods)) stop("methods must be unique.")
  fields <- c("A", "L", "Z", "X", "delta", "J")
  if (!is.data.frame(data) || anyDuplicated(names(data)) ||
      !is.character(columns) || !setequal(names(columns), fields) ||
      length(columns) != length(fields) || anyNA(columns) || anyDuplicated(columns) ||
      !all(columns %in% names(data))) stop("Supply one data column for each of A, L, Z, X, delta, and J.")
  if (!is.character(covariates) || !length(covariates) || anyNA(covariates) ||
      anyDuplicated(covariates) || !all(covariates %in% setdiff(names(data), columns))) {
    stop("covariates must name distinct numeric baseline columns.")
  }
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau <= 0 ||
      !is.numeric(times) || !length(times) || any(!is.finite(times)) ||
      any(times <= 0 | times > tau) || anyDuplicated(times) ||
      any(cfg$b_extra_times > tau)) stop("Use distinct landmarks in (0, tau].")
  if (!whole_number(seed, 0, 2145000000) || !is.environment(learner_env) ||
      !is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("Invalid seed, learner environment, or verbose setting.")
  }
  n <- nrow(data)
  if (n < 10L) stop("At least ten subjects are required for cross-fitting.")
  w <- idid_numeric_frame(data[covariates])
  main <- setNames(data[unname(columns[fields])], fields)
  for (v in setdiff(fields, "X")) {
    if (!(is.numeric(main[[v]]) || is.logical(main[[v]])) ||
        any(!main[[v]] %in% c(0, 1))) stop(v, " must contain only 0 and 1.")
    main[[v]] <- as.integer(main[[v]])
  }
  if (!is.numeric(main$X) || any(!is.finite(main$X) | main$X <= 0 | main$X > tau)) {
    stop("X must be in (0, tau].")
  }
  if (any(main$J > main$delta) || any(main$X == tau & main$delta == 0L) ||
      any(main$X < tau & main$delta == 1L & main$J == 0L)) {
    stop("Inconsistent X, delta, and J. Completed histories end at death or tau.")
  }
  if (length(unique(2L * main$L + main$Z)) != 4L) {
    stop("All four period-instrument cells must be represented.")
  }
  if (!is.list(events) || is.data.frame(events) || length(events) != n)
    stop("events must be a list of recorded event times, one vector per subject.")
  for (i in seq_len(n)) {
    u <- events[[i]]
    if (length(u) && (!is.numeric(u) || any(!is.finite(u) | u <= 0 | u > main$X[i]))) {
      stop("Invalid event times for subject ", i, ". Use (0, X] in the same units as X.")
    }
  }
  # Avoid collisions with internal survival response and logistic response names.
  reserved <- c(fields, "status", "y", "Y", "obsWeights", ".idid_LZ")
  cfg$wcols <- tail(make.unique(c(reserved, make.names(covariates))), length(covariates))
  cfg$xcols <- c("L", "Z", cfg$wcols)
  cfg$qcols <- c("A", cfg$xcols)
  d <- cbind(main, setNames(w, cfg$wcols))
  if (!is.null(fold_id)) {
    if (length(fold_id) != n || anyNA(fold_id)) stop("fold_id must label every row.")
    cfg$fold_id <- match(fold_id, unique(fold_id))
    cfg$outer_folds <- length(unique(cfg$fold_id))
  }
  if (cfg$outer_folds < 2L || cfg$outer_folds > min(n, 100L)) stop("Invalid number of outer folds.")
  cfg$tau <- tau
  cfg$learner_env <- learner_env
  cfg$verbose <- verbose
  cfg$cache <- new.env(parent = emptyenv())
  counts <- vapply(times, function(t) vapply(events, function(u) sum(u <= t), integer(1)),
                   integer(n))
  if ("IPW2" %in% methods) {
    distinct <- sort(unique(unlist(events, use.names = FALSE)))
    if (length(distinct) <= 128L) cfg$event_grid <- distinct
  }
  sim <- list(data = d, events = events, obs_N = counts, landmark_times = times,
              b_times = sort(unique(c(times, cfg$b_extra_times))))
  # Preserve the caller's random-number state across nuisance fits.
  rng <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(restore_rng(rng))
  fit <- setNames(lapply(methods, function(method) {
    estimates <- lapply(times, function(t) get_scores(sim, t, method, seed, cfg))
    matrices <- lapply(c("N", "Y"), function(y) {
      x <- vapply(estimates, `[[`, numeric(n), y)
      colnames(x) <- format(times, digits = 16, trim = TRUE)
      x
    })
    names(matrices) <- c("N", "Y")
    matrices$diagnostics <- setNames(lapply(estimates, `[[`, "diagnostics"),
                                    colnames(matrices$N))
    matrices
  }), methods)
  base <- get("base", envir = cfg$cache)
  structure(list(scores = fit, times = times, tau = tau, covariates = w,
                 fold_id = base$folds, first_stage = base$delta,
                 joint_probability = base$pi,
                 survival_methods = c(censor = cfg$censor, death = cfg$death),
                 control = control), class = "idid_scores")
}

restore_rng <- function(rng) {
  if (!is.null(rng)) {
    assign(".Random.seed", rng, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
}
