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

learn_policy <- function(w,scores,reference,seed,cfg) {
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

policy_accuracy <- function(theta) {
  mass <- pnorm(2)-pnorm(-2)
  F <- function(x) (pnorm(pmin(2,pmax(-2,x)))-pnorm(-2))/mass
  breaks <- c(-2,2,-1.5)
  if(theta[2]!=0) breaks <- c(breaks,-(theta[1]+theta[3]*c(-2,2))/theta[2])
  if(theta[2]+theta[3]!=0) breaks <- c(breaks,(.5*theta[3]-theta[1])/(theta[2]+theta[3]))
  breaks <- sort(unique(pmin(2,pmax(-2,breaks))))
  error <- sum(vapply(seq_len(length(breaks)-1L),function(j) {
    # Coincident boundaries can differ by a few ulps after grid arithmetic.
    # Their probability mass is negligible, and QUADPACK cannot integrate them.
    if(breaks[j+1L]-breaks[j] <= 64*.Machine$double.eps*
       max(1,abs(breaks[j]),abs(breaks[j+1L]))) return(0)
    integrate(function(x) {
      cut <- -(theta[1]+theta[2]*x)/theta[3]; opt <- x-.5
      distance <- abs(F(cut)-F(opt))
      mismatch <- if(theta[3]>0) distance else 1-distance
      mismatch*dnorm(x)/mass
    },breaks[j],breaks[j+1L],rel.tol=1e-9,abs.tol=1e-10)$value
  },numeric(1)))
  pmin(1,pmax(0,1-error))
}
