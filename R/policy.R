policy <- function(theta,w) as.integer(theta[1]+theta[2]*w$V1+theta[3]*w$V2>0)

optimizer_grid_candidates <- function(w, scores, reference, bound = 3,
                                      step = .1, tolerance = 1e-10) {
  n <- nrow(w)
  stopifnot(n > 0L, length(scores$N) == n, length(scores$Y) == n,
            length(reference) == n, all(is.finite(scores$N)), all(is.finite(scores$Y)),
            all(is.finite(w$V1)), all(is.finite(w$V2)),
            all(reference %in% c(0, 1)), is.finite(bound), bound > 0,
            is.finite(step), step > 0)
  points <- sort(unique(c(seq(-bound, bound, by = step), -bound, 0, bound)))
  weights <- cbind(N = as.numeric(scores$N), Y = as.numeric(scores$Y))
  totals <- colSums(weights)
  reference_Y <- mean(reference * weights[, "Y"])
  charts <- list()
  started <- proc.time()[3L]
  for (sign in c(1, -1)) {
    best <- NULL
    for (slope in points) {
      index <- slope * w$V1 + sign * w$V2
      rows <- order(index)
      cumulative <- rbind(c(0, 0), apply(weights[rows, , drop = FALSE], 2L, cumsum))
      cut <- findInterval(-points, index[rows])
      sums <- sweep(-cumulative[cut + 1L, , drop = FALSE], 2L, totals, "+") / n
      constraint <- sums[, "Y"] - reference_Y
      feasible <- constraint >= -tolerance
      value <- ifelse(feasible, sums[, "N"], 1e6 + constraint^2)
      at <- which.min(value)
      candidate <- list(theta = c(points[at], slope, sign),
        objective = unname(sums[at, "N"]), constraint = unname(constraint[at]),
        feasible = unname(feasible[at]), penalized_value = unname(value[at]))
      if (is.null(best) || candidate$penalized_value < best$penalized_value) best <- candidate
    }
    best$rules_evaluated <- length(points)^2
    charts[[as.character(sign)]] <- best
  }
  list(charts = charts, step = step, points_per_axis = length(points),
       rules_evaluated = 2L * length(points)^2,
       elapsed_seconds = unname(proc.time()[3L] - started))
}

learn_two_covariate_policy <- function(w,scores,reference,seed,cfg) {
  best <- NULL; chart <- list(); total_calls <- 0L
  # Screen both charts without outcome-truth information. Retain the original
  # search too, so the extra initialization cannot worsen its empirical optimum.
  grid <- optimizer_grid_candidates(w,scores,reference,bound=cfg$theta_bound,step=.1)
  for(sign in c(1,-1)) {
    calls <- 0L
    objective <- function(beta) {
      theta <- c(beta,sign)
      d <- policy(theta,w)
      constraint <- mean((d-reference)*scores$Y)
      value <- mean(d*scores$N)
      calls <<- calls+1L
      if(constraint >= -1e-10 && (is.null(best) || value<best$objective)) {
        best <<- list(theta=theta,objective=value,constraint=constraint,sign=sign)
      }
      if(constraint < -1e-10) 1e6+constraint^2 else value
    }
    generic <- as.matrix(expand.grid(c(-cfg$theta_bound,0,cfg$theta_bound),
                                    c(-cfg$theta_bound,0,cfg$theta_bound)))
    values <- apply(generic,1,objective)
    starts <- unique(rbind(generic[which.min(values),],
                           grid$charts[[as.character(sign)]]$theta[1:2]))
    runs <- list()
    for(attempt in seq_len(nrow(starts))) {
      start <- starts[attempt,]
      objective(start)
      set.seed(seed+(if(sign>0) 1L else 2L)+100L*(attempt-1L))
      fit <- rgenoud::genoud(objective,nvars=2L,max=FALSE,
        Domains=matrix(rep(c(-cfg$theta_bound,cfg$theta_bound),2),ncol=2,byrow=TRUE),
        starting.values=start,pop.size=cfg$genoud_population,
        max.generations=cfg$genoud_generations,wait.generations=cfg$genoud_wait,
        hard.generation.limit=TRUE,boundary.enforcement=2L,
        gradient.check=FALSE,BFGS=FALSE,print.level=0,MemoryMatrix=TRUE,
        solution.tolerance=1e-6)
      objective(fit$par)
      runs[[attempt]] <- list(start=as.numeric(start),theta=c(fit$par,sign),value=fit$value)
    }
    chart[[as.character(sign)]] <- list(runs=runs,calls=calls)
    total_calls <- total_calls+calls
  }
  if(is.null(best)) stop('No feasible policy found in either normalized parameter chart')
  best$charts <- chart; best$calls <- total_calls; best$grid <- grid
  stopifnot(abs(best$theta[3])==1)
  best
}

idid_score_at <- function(fit, t, method) {
  if (!inherits(fit, "idid_scores")) stop("fit must come from idid_fit_scores().")
  if (!is.numeric(t) || length(t) != 1L || is.na(t) ||
      !is.character(method) || length(method) != 1L || is.na(method)) {
    stop("Specify a single fitted time and method.")
  }
  j <- match(t, fit$times)
  if (is.na(j) || !method %in% names(fit$scores)) stop("Time or method was not fitted.")
  list(N = fit$scores[[method]]$N[, j], Y = fit$scores[[method]]$Y[, j])
}

rule_assignments <- function(rule, w) {
  d <- if (is.function(rule)) rule(w) else
    if (inherits(rule, "idid_rule")) predict(rule, w) else rule
  if (!(is.numeric(d) || is.logical(d)) || !is.null(dim(d)) ||
      length(d) != nrow(w) || any(!d %in% c(0, 1))) {
    stop("A rule must assign one 0/1 treatment per row.")
  }
  as.integer(d)
}

# Search both signs of the coefficient used to normalize the linear index.
learn_general_policy <- function(w, scores, reference, seed, cfg) {
  p <- ncol(w)
  bound <- cfg$theta_bound
  design <- cbind(1, as.matrix(w[, seq_len(p - 1L), drop = FALSE]))
  normalizing <- w[[p]]
  best <- NULL
  for (sign in c(1L, -1L)) {
    objective <- function(beta) {
      d <- as.integer(drop(design %*% beta) + sign * normalizing > 0)
      constraint <- mean((d - reference) * scores$Y)
      value <- mean(d * scores$N)
      if (constraint >= -1e-10 && (is.null(best) || value < best$objective)) {
        best <<- list(theta = c(beta, sign), objective = value, constraint = constraint)
      }
      if (constraint < -1e-10) 1e6 + constraint^2 else value
    }
    candidates <- rbind(rep(0, p), bound * diag(p), -bound * diag(p))
    values <- apply(candidates, 1L, objective)
    set.seed(seed + if (sign > 0) 1L else 2L)
    search <- rgenoud::genoud(objective, nvars = p, max = FALSE,
      Domains = matrix(rep(c(-bound, bound), p), ncol = 2L, byrow = TRUE),
      starting.values = candidates[which.min(values), ],
      pop.size = cfg$genoud_population, max.generations = cfg$genoud_generations,
      wait.generations = cfg$genoud_wait, hard.generation.limit = TRUE,
      boundary.enforcement = 2L, gradient.check = FALSE, BFGS = FALSE,
      print.level = 0, MemoryMatrix = TRUE, solution.tolerance = 1e-6)
    objective(search$par)
  }
  if (is.null(best)) stop("No feasible rule was found in either normalization chart.")
  best
}

idid_learn <- function(fit, t, reference, method = "AIPW",
                       rule_covariates = names(fit$covariates),
                       normalize_by = tail(rule_covariates, 1L), seed = 1L,
                       control = fit$control) {
  scores <- idid_score_at(fit, t, method)
  if (!requireNamespace("rgenoud", quietly = TRUE)) stop("Install rgenoud to learn a rule.")
  if (!is.character(rule_covariates) || !length(rule_covariates) ||
      anyDuplicated(rule_covariates) || !all(rule_covariates %in% names(fit$covariates)) ||
      length(normalize_by) != 1L || !normalize_by %in% rule_covariates ||
      !inherits(control, "idid_control") || !whole_number(seed, 0, 2145000000)) {
    stop("Invalid rule covariates, normalization column, control, or seed.")
  }
  ref <- rule_assignments(reference, fit$covariates)
  order <- c(setdiff(rule_covariates, normalize_by), normalize_by)
  w <- fit$covariates[order]
  rng <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(restore_rng(rng))
  result <- if (length(order) == 2L) {
    learn_two_covariate_policy(setNames(w, c("V1", "V2")), scores, ref, seed, control)
  } else {
    learn_general_policy(w, scores, ref, seed, control)
  }
  # Preserve the operation order used by the search on two-covariate rules.
  structure(list(coefficients = setNames(result$theta, c("(Intercept)", order)),
                 covariates = order, normalized = normalize_by,
                 objective = result$objective, survival_advantage = result$constraint,
                 time = t, method = method), class = "idid_rule")
}

predict.idid_rule <- function(object, newdata, ...) {
  x <- idid_numeric_frame(newdata, object$covariates)
  theta <- unname(object$coefficients)
  if (ncol(x) == 2L) return(policy(theta, setNames(x, c("V1", "V2"))))
  p <- ncol(x)
  index <- drop(cbind(1, as.matrix(x[, seq_len(p - 1L), drop = FALSE])) %*% theta[seq_len(p)])
  as.integer(index + theta[p + 1L] * x[[p]] > 0)
}

idid_value <- function(fit, t, rule, reference, method = "AIPW", conf_level = 0.95) {
  scores <- idid_score_at(fit, t, method)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be between zero and one.")
  }
  difference <- rule_assignments(rule, fit$covariates) -
    rule_assignments(reference, fit$covariates)
  contribution <- cbind(N = difference * scores$N, Y = difference * scores$Y)
  estimate <- colMeans(contribution)
  se <- rep(NA_real_, length(estimate))
  if (method == "AIPW") {
    se <- sqrt(colMeans(sweep(contribution, 2L, estimate, "-")^2) / nrow(contribution))
  }
  half_width <- stats::qnorm((1 + conf_level) / 2) * se
  data.frame(outcome = c("recurrent_count", "survival"), estimate = unname(estimate),
             se = unname(se), lower = unname(estimate) - half_width,
             upper = unname(estimate) + half_width)
}
