#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param bcb_setup
#' @param n_times
#' @param space_kernel
#' @param time_correlation
#' @param sigma
#' @return
#' @author Nick Golding
#' @export
sim_bcb_gp <- function(bcb_setup = bcb_setup,
                       space_kernel,
                       sigma = 1,
                       n_times = 1,
                       time_correlation = 1) {

  # given the bcb setup and number of timesteps to simulate define the latent
  # white noise parameter
  z <- normal(0, 1, dim = c(bcb_setup$n, n_times))

  # evaluate kernel on the bcb space
  stop("TO DO")

  # apply the bcb fft step
  stop("TO DO")

  # maybe apply the AR1 process to this matrix
  if (n_times > 1) {

    stop("TO CHECK")
    f_mat <- ar1ify(f_mat,
                    time_correlation = time_correlation)

  }

  # apply the marginal sigma
  f_mat <- f_mat * sigma

  # attach the white noise as an attribute of this object
  attr(f_mat, "white_noise") <- z

  # return the matrix of spatially- and temporally-correlated results
  f_mat




}
