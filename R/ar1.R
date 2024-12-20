#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param rho
#' @param innovations
#' @return
#' @author Nick Golding
#' @importFrom greta `%*%`
#' @noRd
# create a greta operation for an AR(1) process, given innovations and parameters

# given a scalar greta array for `rho` (between -1 and 1) and a matrix greta
# array of random normal (but not necessarily IID) deviates for 'innovations',
# with columns giving times and rows giving locations, return a matrix of the
# same dimension as `innovations`, with AR(1) process correlations. Do this by
# iterating forward in time.

# Simulate AR1 efficiently in TensorFlow by expanding the iterative equation:
#   X_t = \rho X_{t-1} + epsilon_t
# to
#   X_t = \Sum_{i=0}^{t}(\epsilon_{t - i} \rho ^ i)
# which can be simulated with elementwise matrix operations, and a matrix
# multiplication (50% sparsity, but treated as dense).
ar1 <- function (rho, innovations) {


  # build a matrix of time modifiers
  n_times <- ncol(innovations)
  n_sites <- nrow(innovations)
  t_seq <- seq_len(n_times)
  t_mat <- outer(t_seq, t_seq, FUN = "-")
  t_mat <- pmax(t_mat, 0)

  # which elements to include (don't want upper triangular ones)
  mask <- lower.tri(t_mat, diag = TRUE)

  # matrix of rho ^ n contributions
  rho_mat <- (rho ^ t_mat) * mask

  # multiply by innovations to get multiple correlated timeseries
  # t(greta:::`%*%`(rho_mat, t(innovations)))
  t(rho_mat %*% t(innovations))

}


