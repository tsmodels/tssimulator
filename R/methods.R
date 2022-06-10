#' @rdname add_polynomial
#' @export
#'
#'
add_polynomial <- function(x, ...)
{
  UseMethod("add_polynomial")
}

#' @rdname add_seasonal
#' @export
#'
#'
add_seasonal <- function(x, ...)
{
  UseMethod("add_seasonal")
}

#' @rdname add_arma
#' @export
#'
#'
add_arma <- function(x, ...)
{
  UseMethod("add_arma")
}

#' @rdname add_regressor
#' @export
#'
#'
add_regressor <- function(x, ...)
{
  UseMethod("add_regressor")
}

#' @rdname add_transform
#' @export
#'
#'
add_transform <- function(x, ...)
{
  UseMethod("add_transform")
}


#' @rdname add_anomaly
#' @export
#'
#'
add_anomaly <- function(x, ...)
{
  UseMethod("add_anomaly")
}
