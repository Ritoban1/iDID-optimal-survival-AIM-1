# Complete proposed DGP for the user's iDID policy method.
# Requires base R only. Does not modify or source any original simulation file.
# The response-group construction and its restrictions are documented in
# valid_idid_dgp_assumptions.txt. It is a change to the original outcome DGP.

simulate_valid_idid <- function(n = 2500L, seed = NULL,
                                return_potentials = FALSE) {
  stopifnot(length(n) == 1L, is.finite(n), n >= 1, n == as.integer(n),
            is.logical(return_potentials), length(return_potentials) == 1L,
            !is.na(return_potentials))
  n <- as.integer(n)
  if (!is.null(seed)) set.seed(seed)
  times <- 1:5
  tau <- 10
  size <- 1000L
  cumulative <- function(x) {
    for (j in 2:ncol(x)) x[, j] <- x[, j] + x[, j - 1L]
    colnames(x) <- paste0('t', times)
    x
  }
  draw_bridge <- function(n) {
    # For phi=1/2: F(u)=2*atan(exp(u/2))/pi.
    p <- runif(n, 2*atan(exp(-1))/pi, 2*atan(exp(1))/pi)
    2*log(tan(pi*p/2))
  }
  V1 <- qnorm(runif(n, pnorm(-2), pnorm(2)))
  V2 <- qnorm(runif(n, pnorm(-2), pnorm(2)))
  U0 <- draw_bridge(n)
  U1 <- draw_bridge(n)
  Z <- rbinom(n, 1L, .5)
  L <- rbinom(n, 1L, .5)
  group <- sample.int(3L, n, replace = TRUE, prob = c(.10, .10, .80))

  # Original treatment logits. One uniform per period defines both Z regimes.
  treatment_error0 <- runif(n)
  treatment_error1 <- runif(n)
  A00 <- as.integer(treatment_error0 < plogis(2 + .2*U0 + 2*V1))
  A01 <- as.integer(treatment_error0 < plogis(2 - 7 + .2*U0 + 2*V1))
  A10 <- as.integer(treatment_error1 < plogis(-1.5 - .15*U1 + 1.5*V2))
  A11 <- as.integer(treatment_error1 < plogis(-1.5 + 5 - .15*U1 + 1.5*V2))
  A0 <- ifelse(Z == 1L, A01, A00)
  A1 <- ifelse(Z == 1L, A11, A10)
  A <- ifelse(L == 1L, A1, A0)

  beta <- .8*tanh(4*(V1-V2-.5))
  death_base <- .2 - .01*V1 + .02*V2
  increment_base <- outer(rep(1, n), diff(c(0, times)^1.1) + .1*times) +
    .05*V1 + .03*V2
  # Reuse these independent errors in EVERY potential outcome regime.
  death_error <- rexp(n)
  count_error <- matrix(runif(n*length(times)), nrow = n)
  make_outcome <- function(l, a, z) {
    l <- rep_len(l,n)
    a <- rep_len(a,n)
    z <- rep_len(z,n)
    ul <- ifelse(l == 1L, U1, U0)
    rate <- death_base + (group == 2L)*(.05*l + .03*ul) +
      (group == 1L)*.01*z + (group == 3L)*a*.2125*beta
    mu <- increment_base + (group == 2L)*(.2*l + .05*ul) +
      (group == 1L)*.02*z + (group == 3L)*a*beta
    stopifnot(all(rate > 0), all(mu > 0), all(mu < size))
    D <- death_error/rate
    increments <- matrix(qbinom(as.vector(count_error), size = size,
                               prob = as.vector(mu)/size), nrow = n)
    N_star <- cumulative(increments*outer(D, times, '>='))
    list(D = D, N_star = N_star, increments = increments,
         rate = rate, increment_mean = mu)
  }
  realized <- make_outcome(L, A, Z)
  D <- realized$D

  # Original Weibull-PH censoring law: shape 1, scale 2.
  eta_C <- -3 + .3*A + .1*L + .1*Z + .2*V1 + .1*V2
  C <- 2*rexp(n)/exp(eta_C)
  X <- pmin(C, D, tau)
  delta <- as.integer(C >= pmin(D, tau))
  J <- as.integer(D <= C & D <= tau)
  obs_N <- cumulative(realized$increments*outer(X, times, '>='))
  obs_Y <- outer(X, times, '>=')*1L
  Y_star <- outer(D, times, '>=')*1L
  colnames(obs_Y) <- colnames(Y_star) <- paste0('t', times)

  answer <- list(
    data = data.frame(V1 = V1, V2 = V2, Z = Z, L = L, A = A,
                      X = X, delta = delta, J = J),
    landmark_times = times, tau = tau,
    obs_N = obs_N, obs_Y = obs_Y,
    N_star = realized$N_star, Y_star = Y_star,
    parameters = list(group_probabilities = c(Z = .10, L = .10, A = .80),
                      increment_size = size, optimal_theta = c(.5, -1, 1))
  )
  if (return_potentials) {
    answer$latent <- data.frame(U0 = U0, U1 = U1,
                               group = c('Z','L','A')[group], D = D, C = C)
    answer$treatment_potentials <- data.frame(A00, A01, A10, A11)
    answer$potentials <- list()
    for (l in 0:1) for (a in 0:1) for (z in 0:1) {
      answer$potentials[[paste0(l,a,z)]] <- make_outcome(l,a,z)
    }
  }
  answer
}

# Exact CATEs for DEATH-STOPPED counts and survival, conditional on W=(V1,V2).
# This function integrates out the response group and all unmeasured variables.
idid_true_cate <- function(V1, V2, t) {
  stopifnot(length(t) == 1L, is.finite(t), t >= 0, t <= 10,
            length(V1) > 0L, length(V2) > 0L,
            length(V1) == length(V2) || length(V1) == 1L || length(V2) == 1L,
            all(is.finite(V1)), all(is.finite(V2)),
            all(abs(V1) <= 2), all(abs(V2) <= 2))
  n <- max(length(V1), length(V2))
  V1 <- rep_len(V1, n)
  V2 <- rep_len(V2, n)
  beta <- .8*tanh(4*(V1-V2-.5))
  rate0 <- .2 - .01*V1 + .02*V2
  rate1 <- rate0 + .2125*beta
  times <- (1:5)[(1:5) <= t]
  effect_N <- numeric(n)
  if (length(times)) {
    mu0 <- outer(rep(1, n), diff(c(0, 1:5)^1.1)[times] + .1*times) +
      .05*V1 + .03*V2
    effect_N <- rowSums(exp(-outer(rate1, times))*(mu0 + beta) -
                         exp(-outer(rate0, times))*mu0)
  }
  data.frame(tau_N = .8*effect_N,
             tau_Y = .8*(exp(-rate1*t) - exp(-rate0*t)))
}

# Population contrasts by deterministic quadrature over the true truncated-N W.
# A NULL reference_theta denotes the fixed reference rule d_ref(W)=0.
# Both supplied theta vectors use d(W)=I(theta0+theta1*V1+theta2*V2 > 0).
idid_true_policy_contrast <- function(theta = c(.5,-1,1),
                                      reference_theta = NULL, t = 1,
                                      rel.tol = 1e-8) {
  stopifnot(length(theta) == 3L, all(is.finite(theta)),
            is.null(reference_theta) ||
              (length(reference_theta) == 3L && all(is.finite(reference_theta))),
            length(t) == 1L, is.finite(t), t >= 1, t <= 10,
            is.finite(rel.tol), rel.tol > 0)
  policies <- list(theta)
  if (!is.null(reference_theta)) policies <- c(policies, list(reference_theta))
  mass <- pnorm(2)-pnorm(-2)
  inside <- function(x) x[is.finite(x) & x > -2 & x < 2]
  outer_breaks <- c(-2,2)
  for (th in policies) if (th[2] != 0) {
    outer_breaks <- c(outer_breaks, inside(-(th[1]+th[3]*c(-2,2))/th[2]))
  }
  outer_breaks <- sort(unique(outer_breaks))
  solve_component <- function(component) {
    at_V1 <- function(w1) {
      inner_breaks <- c(-2,2)
      for (th in policies) if (th[3] != 0) {
        inner_breaks <- c(inner_breaks, inside(-(th[1]+th[2]*w1)/th[3]))
      }
      inner_breaks <- sort(unique(inner_breaks))
      ans <- 0
      for (j in seq_len(length(inner_breaks)-1L)) {
        lower <- inner_breaks[j]
        upper <- inner_breaks[j+1L]
        mid <- (lower+upper)/2
        d <- as.integer(theta[1]+theta[2]*w1+theta[3]*mid > 0)
        ref <- if (is.null(reference_theta)) 0L else {
          as.integer(reference_theta[1]+reference_theta[2]*w1+reference_theta[3]*mid > 0)
        }
        if (d != ref) {
          ans <- ans + (d-ref)*integrate(function(w2) {
            idid_true_cate(w1,w2,t)[[component]]*dnorm(w2)/mass
          },lower,upper,rel.tol=rel.tol,abs.tol=1e-11)$value
        }
      }
      ans*dnorm(w1)/mass
    }
    sum(vapply(seq_len(length(outer_breaks)-1L), function(j) {
      integrate(function(w1) vapply(w1,at_V1,numeric(1)),
                outer_breaks[j],outer_breaks[j+1L],
                rel.tol=rel.tol,abs.tol=1e-10)$value
    },numeric(1)))
  }
  c(value_contrast = solve_component('tau_N'),
    survival_contrast = solve_component('tau_Y'))
}
