#' @keywords internal
#' @import data.table
#' @import methods
#' @importFrom zoo coredata index as.zoo plot.zoo
#' @importFrom xts xts
#' @importFrom tsaux future_dates fourier_series tstransform
#' @importFrom stats arima.sim rnorm na.omit sd filter
#' @importFrom utils tail
#' @importFrom graphics grid lines
#' @import tsmethods
#' @examples
#' library(magrittr)
#' set.seed(105)
#' sim <- initialize_simulator(rnorm(12*20, 0, 1), sampling = "months")
#' sim <- sim %>% add_polynomial(order = 2, alpha = 0.1, beta = 0.02, l0 = 100, b0 = 0.25) %>%
#' add_seasonal(frequency = 12, gamma = 0.01, init_harmonics = 4, init_scale = 3)
#' # plot(sim, lwd = 2, ylab = "")
#' # plot(sim, type = "components")
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
