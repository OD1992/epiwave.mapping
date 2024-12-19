#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param distance
#' @param correlation
#' @param nu
#' @return
#' @author Nick Golding
#' @importFrom fields Matern.cor.to.range
#' @export
matern_lengthscale_from_correl <- function(distance, correlation, nu) {

  # use fields to get the range parameter for their parameterisation
  range <- fields::Matern.cor.to.range(d = distance,
                                       nu = nu,
                                       cor.target = correlation)

  # convert this to the greta.gp definition (same as INLA, notwithstanding their
  # parameterisation). This was worked out from Wang et al. with much pain and
  # suffering
  0.5 * sqrt(8 * nu) * range

}

# # check deinfition vs fields
# arange <- 40000
# nu <- 5/2
# len <- 0.5 * sqrt(8 * nu) * arange
# k <- mat52(lengthscales = len, variance = 1)
# d <- c(0, 10^(1:5))
# calculate(k(d)[1, ])[[1]][1, ]
# fields::stationary.cov(d,
#                        Covariance = "Matern",
#                        range = arange,
#                        smoothness = nu)[1, ]
