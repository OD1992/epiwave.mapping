#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param upper_value
#' @param prob_higher
#' @param scale
#' @return
#' @author Nick Golding
#' @export

# Find the shape and rate parameters of a gamma distribution that meets the
# criteria that `lower_prob` proportion of the mass is lower than
# `lower_value`, and `upper_prob` proportion of the mass is greater than
# `upper_value`.
find_gamma_parameters <- function(lower_value, lower_prob, upper_value, upper_prob) {
  # Since there is no analytical form for the inverse of the gamma function
  # (used in the CDF), we may as well do 1D optimisation

  param_transform <- function(free_params) {
    params <- exp(free_params)
    list(shape = params[1], rate = params[2])
  }

  # goodness of fit metric (squared error re. upper_value)
  gof <- function(free_params) {
    params <- param_transform(free_params)
    # given the parameters, look up the quantile at (ie. proportion of mass
    # below) the given lower bound probability
    lower_value_est <- qgamma(p = lower_prob,
                              shape = params$shape,
                              rate = params$rate)
    # do the same for the upper bound, but flip the probability
    upper_value_est <- qgamma(p = 1 - upper_prob,
                              shape = params$shape,
                              rate = params$rate)
    # return the SSE versus the target values
    errors <- c(lower_value_est, upper_value_est) - c(lower_value, upper_value)
    sum(errors ^ 2)
  }

  # optimise these
  init <- c(0, 0)
  o <- optim(init, gof)

  # convert back to the parameter values, and return
  param_transform(o$par)

}
