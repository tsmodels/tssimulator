#' Ensemble Setup
#'
#' @param ... either a list of valid simulation objects or individual objects
#' passed to the function
#' @return A object of class tssim.mixture ready for ensembling,
#' @details The function performs certain checks on the inputs to ensure they
#' conform to the simulation models in the package and are of the same length.
#' @aliases mixture_modelspec
#' @rdname mixture_modelspec
#' @export
#'
#'
mixture_modelspec <- function(...)
{
  m <- list(...)
  # if its already a list
  if (length(m) == 1) {
    if (is.list(m[[1]])) {
      m <- m[[1]]
    }
  }
  object_class <- sapply(1:length(m), function(i) tail(class(m[[i]]),1))
  if (!(all(object_class %in% c("tssim.component")))) stop("\nThe objects must all inherit tssim.component class.")
  sim_length <- sapply(1:length(m), function(i) m[[i]]$n)
  if (!all(sim_length == sim_length[1])) stop("\nsimulation objects have different lengths.")
  class(m) <- "tssim.mixture"
  return(m)
}

fun_roc <- function(x, log = FALSE)
{
  if (log) {
    out <- na.omit(diff(log(x)))
  } else {
    out <- x[2:length(x)]/x[1:(length(x) - 1)] - 1
  }
  return(out)
}

#' Ensembling of Simulations
#'
#' @param object an object of class tssim.mixture.
#' @param weights the weighting (or probability) matrix for aggregating the
#' simulations (see details).
#' @param difference whether to take the rates of changes first before aggregating
#' and reconverting to levels.
#' @param ... additional parameters.
#' @return A vector of the simulated series.
#' @method tsensemble tssim.mixture
#' @details When mixing dynamics for the same series, and when series are
#' not stationary, differences should be used. In that case the rate of change
#' transformation is applied to each simulated series and then weighted by the
#' weights matrix. Since the weights matrix will have one more row than is required (the
#' first row), this can be used to choose how the initial level is generated.
#' For instance, if we want to use the level of the first simulated series, then
#' the first row would have a 1 on the first column and zeros in the rest.
#' For aggregating series, difference should be set to FALSE since we are looking
#' at summation of data (under the assumption of flow variables). In this case,
#' the p matrix is usually static by column (i.e. the same weights).
#' @aliases tsensemble
#' @rdname tsensemble
#' @export
#'
#'
tsensemble.tssim.mixture <- function(object, weights = NULL, difference = TRUE, ...)
{
  n <- length(object)
  if (!is.matrix(weights)) stop("\nweights must be a matrix")
  required_n <- object[[1]]$n
  if (nrow(weights) != required_n) stop("\nnrow of p must be equal to the simulation length of the series")
  if (ncol(weights) != n) stop("\nncol of p must be equal to length(object)")
  y <- do.call(cbind, lapply(1:length(object), function(i) object[[i]]$simulated))
  if (difference) {
    l0 <- sum(weights[1,] * y[1,])
    r <- do.call(cbind, lapply(1:length(object),function(i) fun_roc(y[,i])))
    r <- rowSums(r * weights[-1,])
    new_sim <- cumprod(c(l0, 1 + r))
    if (!is.null(object[[1]]$index)) new_sim <- xts(new_sim, object[[1]]$index)
  } else {
    new_sim <- rowSums(weights * y)
    if (!is.null(object[[1]]$index)) new_sim <- xts(new_sim, object[[1]]$index)
  }
  return(new_sim)
}
