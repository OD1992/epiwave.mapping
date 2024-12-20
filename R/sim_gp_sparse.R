#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param inducing_points
#' @param grid_raster
#' @param n_times
#' @param space_kernel
#' @param time_correlation
#' @param sigma
#' @return
#' @author Nick Golding
#' @export
sim_gp_sparse <- function(inducing_points,
                          grid_raster,
                          space_kernel,
                          sigma = 1,
                          n_times = 1,
                          time_correlation = 1,
                          tol = 1e-4) {

  # pull out the cell coordiantes from the raster
  cell_coords <- xyFromCell(grid_raster, cells(grid_raster))

  # given the inducing point setup and number of times to simulate, use greta.gp
  # to create the spatially-correlated innovations across all non-NA cells
  f <- greta.gp::gp(x = cell_coords,
                    kernel = space_kernel,
                    inducing = inducing_points,
                    n = n_times,
                    tol = tol)

  # pull out the white noise variables
  gp_info <- attr(f, "gp_info")
  v <- gp_info$v

  # maybe apply the AR1 process to this matrix
  if (n_times > 1) {

    f <- ar1(
      rho = time_correlation,
      innovations = f
    )

  }

  # apply the marginal variance
  f <- f * sigma

  # attach the white noise variables as an attribute of this object, for
  # initialisation later
  attr(f, "v") <- v

  # return the matrix of spatially- and temporally-correlated results
  f

}
