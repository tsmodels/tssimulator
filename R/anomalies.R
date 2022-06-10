#' Anomaly Component
#'
#' @param x an object of class tsissm.component or other supported class.
#' @param time the numeric index of when the anomaly occurs. If NULL, a ranom
#' time will be chosen.
#' @param delta the autoregressive component determining the type of anomaly. A
#' value of zero results in an additive outlier, a value of 1 in a level shift
#' and anything in between a temporary change with a half life of -log(2)/log(delta).
#' @param ratio the anomaly to series ratio at the time it occurs. For instance, a
#' value of 1 means that the anomaly will jump by 100 percent compared to the
#' data series.
#' @param ... additional parameters.
#' @return An object of class tsissm.component updated with the anomaly
#' component.
#' @method add_anomaly tsissm.component
#' @aliases add_anomaly
#' @rdname add_anomaly
#' @export
#'
#'
add_anomaly.tsissm.component <- function(x, time = NULL, delta = 0, ratio = 0.5, ...)
{
  n <- 1:x$n
  if (!is.null(x$extra$anomaly_times)) {
    a_times <- x$extra$anomaly_times
    if (is.null(time)) {
      valid_seq <- setdiff(n, a_times)
      time <- sample(valid_seq, 1)
    } else {
      time <- as.integer(time)
      if (!time %in% n) {
        stop("\ninvalid time index (outside of data index)")
      }
    }
  } else {
    if (is.null(time)) {
      time <- sample(n, 1)
    } else {
      time <- as.integer(time)
      if (!time %in% n) {
        stop("\ninvalid time index (outside of data index)")
      }
    }
  }
  anomaly <- rep(0, x$n)
  anomaly[time] <- 1
  anomaly <- as.numeric(filter(anomaly, filter = delta, method = "recursive"))
  anomaly <- matrix(c(0, anomaly * x$simulated * ratio), ncol = 1)
  if (delta == 0) {
    anomaly_name <- "AO"
    a_name <- "additive_outlier"
  } else if (delta == 1) {
    anomaly_name <- "LS"
    a_name <- "level_shift"
  } else {
    anomaly_name <- "TC"
    a_name <- "temporary_change"
  }
  anomaly_name <- paste0(anomaly_name, time)
  colnames(anomaly) <- anomaly_name
  x$components <- cbind(x$components, anomaly)
  x$table <- rbind(x$table,
                   data.table(component = a_name, type = "anomaly", start = 2, include = TRUE, parameter = "ratio", value = ratio))
  x$extra$anomaly_times <- c(x$extra$anomaly_times, time)
  x$extra$anomaly_delta <- c(x$extra$anomaly_delta, delta)
  x$simulated <- x$simulated + anomaly[-1,1]
  return(x)
}
