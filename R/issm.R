#' Simulator Initializer
#'
#' @param x a vector of zero mean errors to use in the model.
#' @param index an optional Date or POSIXct vector of same length as x. Used
#' for indexing the simulated values.
#' @param sampling an optional string denoting the sampling frequency for the
#' simulator. If no index is present, will automatically generate one based on
#' the sampling frequency given with start date 2000-01-01. Valid sampling
#' frequencies are days, weeks, months, years, secs, mins, hours and subintervals
#' of those as documented in the \code{\link{seq.POSIXt}} function.
#' @param model the type of model to initialize a class for.
#' @param ... additional parameters to the function (not currently used).
#' @return A object whose class depends on the type of model used.
#' @aliases initialize_simulator
#' @rdname initialize_simulator
#' @export
#'
initialize_simulator <- function(x, index = NULL, sampling = NULL, model = "issm", ...)
{
  # start with issm model
  model <- match.arg(model[1], c("issm"))
  n <- length(x)
  M <- matrix(c(0, as.numeric(x)), ncol = 1)
  colnames(M) <- "Error"
  simulated <- as.numeric(x)
  if (is.null(index)) {
    if (!is.null(sampling)) {
      index <- future_dates(as.Date("2000-01-01"), sampling, length(x))
    }
  } else {
    if (length(index) != length(x)) stop("\nindex must have the same length as x")
  }
  components_table <- data.table(component = "error", type = "irregular", start = 2, include = TRUE, parameter = "sigma", value = sd(x))
  L <- list(components = M, simulated = simulated, n = n, table = components_table, index = index, sampling = sampling)
  class(L) <- c(switch(model[1], "issm" = "tsissm.component"), "tssim.component")
  return(L)
}

#' Polynomial Trend Component
#'
#' @param x an object of class tsissm.component or other supported class.
#' @param order the order of the polynomial (min 1 and max 2).
#' @param alpha the decay coefficient on the error of the level.
#' @param beta the decay coefficient on the error of the slope.
#' @param phi dampening parameter for the slope.
#' @param l0 initial level.
#' @param b0 initial slope.
#' @param ... additional parameters.
#' @return An object of class tsissm.component updated with the polynomial trend
#' component.
#' @method add_polynomial tsissm.component
#' @aliases add_polynomial
#' @rdname add_polynomial
#' @export
#'
#'
add_polynomial.tsissm.component <- function(x, order = 1, alpha = 0.1, beta = 0.01, phi = 1.0, l0 = 100, b0 = 1.0, ...)
{
  # this component allows replacement in order to re-run in the presence
  # of seasonal normalization (which requires an adjustment to the dynamics)
  type <- NULL
  n <- x$n + 1
  extra_args <- list(...)
  if (any(names(extra_args) == "a")) {
    a <- extra_args$a
  } else {
    a <- rep(0, n)
  }
  e <- x$components[,"Error"]
  v <- rep(0.0, n)
  level <- slope <- rep(0.0, n)
  level[1] <- l0
  if (order == 1) {
    slope[1] <- 0.0
    beta <- 0.0
    phi <- 1.0
  } else {
    slope[1] <- b0
  }
  for (i in 2:n) {
    level[i] <- level[i - 1] + phi * slope[i - 1] + alpha * e[i] + a[i]
    slope[i] <- phi * slope[i - 1] + beta * e[i]
    v[i] <- level[i - 1] + phi * slope[i - 1]
  }
  poly_mat <- matrix(level, ncol = 1)
  colnames(poly_mat) <- "Level"
  if (order == 2) {
    poly_mat <- cbind(poly_mat, matrix(slope, ncol = 1))
    colnames(poly_mat)[2] <- "Slope"
  }
  if (NROW(x$table[type == "polynomial"]) > 0) {
    # replace current
    x$components <- x$components[, -which(colnames(x$components) %in% c("Level","Slope")), drop = FALSE]
    x$components <- cbind(x$components, poly_mat)
    x$simulated <- x$simulated  - x$extra$trend[-1] + v[-1]
    components_table <- data.table(component = "level", type = "polynomial", start = 1, include = TRUE, parameter = "alpha", value = alpha)
    if (order == 2) {
      components_table <- rbind(components_table,
                                data.table(component = "slope", type = "polynomial", start = 1, include = TRUE, parameter = "beta", value = beta))
      components_table <- rbind(components_table,
                                data.table(component = "dampening", type = "polynomial", start = 1, include = TRUE, parameter = "phi", value = phi))
    }
    x$table <- x$table[-which(x$table$type %in% "polynomial")]
    x$table <- rbind(x$table, components_table)
  } else {
    x$components <- cbind(x$components, poly_mat)
    x$simulated <- x$simulated + v[-1]
    components_table <- data.table(component = "level", type = "polynomial", start = 1, include = TRUE, parameter = "alpha", value = alpha)
    if (order == 2) {
      components_table <- rbind(components_table,
                                data.table(component = "slope", type = "polynomial", start = 1, include = TRUE, parameter = "beta", value = beta))
      components_table <- rbind(components_table,
                                data.table(component = "dampening", type = "polynomial", start = 1, include = TRUE, parameter = "phi", value = phi))
    }
    x$table <- rbind(x$table, components_table)
  }
  x$extra <- list(seasonal_normalization = a, trend = v)
  return(x)
}

seasonal_matrix <- function(frequency)
{
  if (frequency < 2) stop("seasonal.matrix: frequency must be > 1")
  Fmat <- matrix(0, frequency, frequency)
  Fmat[1, frequency] <- 1
  Fmat[2:frequency, 1:(frequency - 1)] <- diag(frequency - 1)
  return(Fmat)
}

#' Seasonal Trend Component
#'
#' @param x an object of class tsissm.component or other supported class.
#' @param frequency seasonal frequency.
#' @param gamma the decay coefficient on the error of the seasonal component
#' @param s0 a vector of length frequency - 1 for the initial seasonal component.
#' @param init_harmonics number of harmonics to initialize s0 when this is not provided.
#' @param normalized_seasonality whether normalize the seasonal component based on the
#' method of Roberts and McKenzie. This is applied only to a single seasonal frequency.
#' @param init_scale the scaling multiplier for s0 (when this is not provided).
#' @param ... additional parameters.
#' @return An object of class tsissm.component updated with the seasonal
#' component.
#' @method add_seasonal tsissm.component
#' @aliases add_seasonal
#' @rdname add_seasonal
#' @export
#'
#'
add_seasonal.tsissm.component <- function(x, frequency = 12, gamma = 0.01, s0 = NULL,
                                         init_harmonics = frequency/2,
                                         normalized_seasonality = TRUE,
                                         init_scale = 1, ...)
{
  component <- NULL
  if (any(x$table$component == paste0("seasonal_",frequency))) stop("\nseasonal component with this frequency already exists.")
  e <- x$components[,"Error"]
  if (any(x$table$type  == "seasonal")) normalized_seasonality <- FALSE
  if (normalized_seasonality) {
    # we need to rebuild the polynomial
    if (any(x$table$type == "polynomial")) {
        alpha <- x$table[component == "level"]$value
        beta <- 1e-12
        phi <- 1
        l0 <- x$components[1,"Level"]
        b0 <- 0
        order <- 1
        if (any(x$table$component == "slope")) {
          order <- 2
          b0 <- x$components[1,"Slope"]
          phi <- x$table[component == "dampening"]$value
          beta <- x$table[component == "slope"]$value
        }
        x <- add_polynomial.tsissm.component(x, order = order, alpha = alpha, beta = beta, phi = phi, l0 = l0, b0 = b0,
                                            check_bounds = TRUE, a = gamma/frequency * e)
    }
  }
  n <- x$n + 1
  if (!is.null(s0)) {
    if (length(s0) != (frequency - 1)) stop("\ns0 must be of length frequency - 1")
    s0 <- c(s0, -1.0 * sum(s0))
  } else {
    h <- fourier_series(dates = as.Date(1:(frequency*2), origin = as.Date("1970-01-01")), period = frequency, K = init_harmonics)
    s0 <- (init_scale * h %*% rnorm(ncol(h), 0, 1))[1:(frequency - 1),]
    s0 <- c(s0, -1.0 * sum(s0))
  }
  Fmat <- seasonal_matrix(frequency)
  s <- matrix(0, ncol = frequency, nrow = n)
  s[1,] <- s0
  colnames(s) <- paste0("S",0:(frequency - 1))
  v <- rep(0, n)
  for (i in 2:n) {
    if (normalized_seasonality) {
      a <- gamma/frequency * e[i]
    } else{
      a <- 0
    }
    sx <- (s[i - 1, frequency] + gamma * e[i]) - a
    s[i,] <- (Fmat %*% s[i - 1, ] ) - a
    s[i, 1] <- sx
  }
  seasonal <- s[, frequency, drop = FALSE]
  colnames(seasonal) <- paste0("Seasonal_",frequency)
  x$components <- cbind(x$components, seasonal)
  x$simulated <- x$simulated + seasonal[1:(n - 1)]
  components_table <- data.table(component = paste0("seasonal_",frequency), type = "seasonal", start = 1, include =  TRUE, parameter = "gamma", value = gamma)
  x$table <- rbind(x$table, components_table)
  return(x)
}

#' ARMA Component
#'
#' @param x an object of class tsissm.component or other supported class.
#' @param order the ar and ma orders.
#' @param ar a vector of ar coefficients.
#' @param ma a vector of ma coefficients.
#' @param mu the mean parameter (defaults to zero) of the ARMA process.
#' @param ... additional parameters.
#' @return An object of class tsissm.component updated with the ARMA
#' component.
#' @method add_arma tsissm.component
#' @aliases add_arma
#' @rdname add_arma
#' @export
#'
#'
add_arma.tsissm.component <- function(x, order = c(0,0), ar = 0, ma = 0, mu = 0, ...)
{
  if (any(x$table$component == "arma")) stop("\narma component already exists.")
  if (length(order) != 2) stop("\norder must be a vector of length 2")
  n <- x$n
  order[1] <- max(0, order[1])
  order[2] <- max(0, order[2])
  if (sum(order) == 0) {
    if (mu != 0) {
      sim <-  matrix(mu, ncol = 1, nrow = n + 1)
      colnames(sim) <- "ARMA"
      x$components <- cbind(x$components, sim)
      x$simulated <- x$simulated + sim[-1]
      x$table <- rbind(x$table, data.table(component = "mu", type = "arma", start = 2, include = TRUE, parameter = "mu", value = mu))
    }
    return(x)
  }
  if (length(ar) != order[1]) stop(paste0("ar must have ", order[1]," parameters"))
  if (length(ma) != order[2]) stop(paste0("ma must have ", order[2]," parameters"))
  sim <- arima.sim(n = n, model = list(order = c(order[1],0,order[2]), ar = ar, ma = ma), innov = x$components[-1,"Error"]) + mu
  sim <- matrix(c(0, as.numeric(sim)), ncol = 1)
  colnames(sim) <- "ARMA"
  x$components <- cbind(x$components, sim)
  x$simulated <- x$simulated + sim[-1]
  if (order[1] > 0) {
    ar_table <- do.call(rbind, lapply(1:order[1], function(i){
      data.table(component = paste0("ar",i), type = "arma", start = 2, include = TRUE, parameter = "theta", value = ar[i])
    }))
    x$table <- rbind(x$table, ar_table)
  }
  if (order[2] > 0) {
    ma_table <- do.call(rbind, lapply(1:order[2], function(i){
      data.table(component = paste0("ma",i), type = "arma", start = 2, include = TRUE, parameter = "psi", value = ma[i])
    }))
    x$table <- rbind(x$table, ma_table)
  }
  x$table <- rbind(x$table, data.table(component = "mu", type = "arma", start = 2, include = TRUE, parameter = "mu", value = mu))
  return(x)
}

#' Regressor Component
#'
#' @param x an object of class tsissm.component or other supported class.
#' @param xreg a matrix of regressors.
#' @param pars regressors coefficients.
#' @param ... additional parameters.
#' @return An object of class tsissm.component updated with the regressor
#' components.
#' @method add_regressor tsissm.component
#' @aliases add_regressor
#' @rdname add_regressor
#' @export
#'
#'
add_regressor.tsissm.component <- function(x, xreg = NULL, pars = NULL, ...)
{
  if (is.null(xreg)) return(x)
  if (is.null(pars)) stop("\npars cannot be NULL")
  n <- x$n
  if (NROW(xreg) != n) stop(paste0("\nxreg must have ", n," rows"))
  xreg <- as.matrix(xreg)
  m <- NCOL(xreg)
  if (length(pars) != m) stop("\npars must be of length NCOL(xreg)")
  xr <- xreg %*% pars
  x$simulated <- x$simulated + xr
  if (is.null(colnames(xreg))) {
    xreg_names <- paste0("x_",1:m)
    colnames(xreg) <- xreg_names
  } else {
    xreg_names <- colnames(xreg)
  }
  xreg <- rbind(matrix(0, ncol = m, nrow = 1), xreg)
  X <- matrix(c(0, xr), ncol = 1)
  colnames(X) <- "X"
  x$components <- cbind(x$components, X)
  xreg_table <- do.call(rbind, lapply(1:m, function(i){
    data.table(component = xreg_names[i], type = "regressor", start = 2, include = TRUE, parameter = "rho", value = pars[i])
  }))
  x$table <- rbind(x$table, xreg_table)
  if (!is.null(x$extra$xreg)) {
    x$extra$xreg <- cbind(x$extra$xreg, xreg)
  } else {
    x$extra$xreg <- xreg
  }
  return(x)
}

#' Transform
#'
#' @param x an object of class tsissm.component or other supported class.
#' @param method a valid transform.
#' @param lambda the Box-Cox parameter.
#' @param lower the lower bound for the logit transform.
#' @param upper the upper bound for the logit transform.
#' @param ... additional parameters.
#' @return An object of class tsissm.component updated with the transformation.
#' @method add_transform tsissm.component
#' @aliases add_transform
#' @details The inverse transform is applied to the simulated series based on
#' the \code{\link{tstransform}} method.
#' @rdname add_transform
#' @export
#'
#'
add_transform.tsissm.component <- function(x, method = "box-cox", lambda = 1, lower = 0, upper = 1, ...)
{
  if (any(x$table$type == "transform")) stop("\ntransform already included in object.")
  if (method == "box-cox") {
    t_table <- data.table(component = method[1], type = "transform", start = 1, include = TRUE,
                 parameter = "lambda", value = lambda)
  } else if (method == "logit") {
    t_table <- rbind(
      data.table(component = method[1], type = "transform", start = 1, include = TRUE,
                 parameter = "lower", value = lower),
      data.table(component = method[1], type = "transform", start = 1, include = TRUE,
                 parameter = "upper", value = upper))
  }
  x$table <- rbind(x$table, t_table)
  tr <- tstransform(method = method[1], lambda = lambda, lower = lower, upper = upper)
  y <- tr$inverse(x$simulated, lambda = lambda, lower = lower, upper = upper)
  x$simulated <- y
  return(x)
}

#' State Decomposition
#'
#' @param object an object of class tsissm.component or other supported class.
#' @param ... additional parameters.
#' @return A matrix of the simplified state decomposition.
#' @method tsdecompose tsissm.component
#' @details Creates a simplified decomposition of the states and aligns their
#' time indices so that the sum up to the simulated component per period.
#' @aliases tsdecompose
#' @rdname tsdecompose
#' @export
#'
#'
tsdecompose.tsissm.component <- function(object, ...)
{
  components <- object$components
  cnames <- colnames(components)
  n <- nrow(components)
  Irregular <- components[2:n,"Error"]
  new_components <- NULL
  if (any(cnames == "Level")) {
    Trend <- components[1:(n - 1),"Level"]
    if (any(cnames == "Slope")) {
      Trend <- Trend + components[1:(n - 1),"Slope"]
    }
    new_components <- cbind(new_components, Trend)
  }
  if (any(grepl("Seasonal",cnames))) {
    wix <- which(grepl("Seasonal", cnames))
    Seasonal <- 0
    for (i in 1:length(wix)) Seasonal <- Seasonal + components[1:(n - 1),wix[i]]
    Seasonal <- matrix(Seasonal, ncol = 1)
    colnames(Seasonal) <- "Seasonal"
    new_components <- cbind(new_components, Seasonal)
  }
  if (any(cnames == "X")) {
    X <- components[2:n,"X"]
    new_components <- cbind(new_components, X)
  }
  if (any(cnames  == "ARMA")) {
    Irregular <- Irregular + components[2:n,"ARMA"]
  }
  new_components <- cbind(new_components, Irregular)
  if (!is.null(object$index)) {
    new_components <- xts(new_components, object$index)
  }
  return(new_components)
}


#' Plot Simulation Object
#'
#' @param x an object of class tsissm.component or other supported class.
#' @param y not used.
#' @param type the type of output to plot.
#' @param ... additional parameters passed to the \code{\link{plot.zoo}} function.
#' @method plot tsissm.component
#' @aliases plot
#' @rdname plot
#' @export
#'
#'
plot.tsissm.component <- function(x, y = NULL, type = c("simulated","components"), ...)
{
  type <- match.arg(type[1], c("simulated","components"))
  if (type == "simulated") {
    y <- x$simulated
    if (!is.null(x$index)) {
      y <- xts(y, x$index)
    }
  } else {
    y <- tsdecompose(x)
  }
  simulated <- as.zoo(y)
  plot(simulated, plot.type	 = "multiple", xlab = "", ...)
  grid()
}
