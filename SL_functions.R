###
SL.lgb <- function(Y, X, newX, family, obsWeights, id,
                   depth = 3, rate = 0.01, rounds = 200) {
  stopifnot(require("lightgbm"))
  
  params <- list(task = "train",
                 max_depth = depth,
                 learning_rate = rate)
  
  if (family$family == 'gaussian') {
    params$objective <- "regression"
  } else if (family$family == 'binomial') {
    params$objective <- "binary"
  }
  
  m <- lgb.train(params,
                 lgb.Dataset(data = as.matrix(X), label = Y, weight = obsWeights),
                 nrounds = rounds, verbose = -1)
  
  pred <- predict(m, as.matrix(newX))
  
  fit = list(object = m)
  class(fit) <- 'SL.lgb'
  out <- list(pred = pred, fit = fit)
  return(out)
}

predict.SL.lgb <- function(object, newdata, ...) {
  stopifnot(require("lightgbm"))
  pred <- predict(object$object, as.matrix(newdata))
  return(pred)
}
# Define HAL function for SuperLearner
SL.hal <- function(Y, X, newX, family, obsWeights, id,
                   num_knots = 10, max_degree = 1) {
  stopifnot(require("hal9001"))
  params=NULL
  if (family$family == 'gaussian') {
    params$objective <- "regression"
  } else if (family$family == 'binomial') {
    params$objective <- "binary"
  }
  
  m <- fit_hal(as.matrix(X), Y,
               num_knots = num_knots,
               max_degree = max_degree)
  
  pred <- predict(m, as.matrix(newX))
  
  fit = list(object = m)
  class(fit) <- 'SL.hal'
  out <- list(pred = pred, fit = fit)
  return(out)
}

predict.SL.hal <- function(object, newdata, ...) {
  stopifnot(require("hal9001"))
  pred <- predict(object$object, as.matrix(newdata))
  return(pred)
}