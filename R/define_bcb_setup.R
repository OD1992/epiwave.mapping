#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param grid_raster a terra spatRaster object for which to create a
#'   block-circulant basis representation. This should have a reasonable
#'   projected CRS, so that cells can be assumed to be in a regular grid on a
#'   two-dimensional plane
#' @param check_lonlat whether to error if the user provides a raster
#'   that appears to have a longitude/latitude (and therefore unprojected) CRS,
#'   as tested by `terra::is.lonlat()`
#' @return
#' @author Nick Golding
#' @export
define_bcb_setup <- function(grid_raster, check_lonlat = TRUE) {

  # check this raster doesn't look unprojected
  looks_lonlat <- terra::is.lonlat(grid_raster)
  if (looks_lonlat & check_lonlat) {
    stop("grid_raster appears to have a longitude/latitude CRS rather than a
         projected CRS. Please apply a reasonable projection, or set
         check_lonlat = FALSE if this is incorrect")
  }

  # stash the original, to preserve
  original_raster <- grid_raster

  # remove NA space around the edges
  grid_raster <- trim(grid_raster)

  # get the dimensions
  res <- res(grid_raster)
  m <- ncol(grid_raster)
  n <- nrow(grid_raster)

  # create an expanded grid (representing a mapping to locations on the surface
  # of a torus) that has dimensions that are a power of 2, and are at least
  # twice as big in each dimension as the original grid. This ensures that
  # simulation of the GP is exact (given the planar regular grid assumption),
  # since the minimum distance between two points is across the plain (within
  # the original grid), rather than around the back of the torus.
  M <- ceiling2(2 * m)
  N <- ceiling2(2 * n)

  # get the distances from each cell to the centre of the extended torus grid
  x_sq <- centred_squared_distance(M, res[1])
  y_sq <- centred_squared_distance(N, res[2])
  dist_squared <- outer(x_sq, y_sq, FUN = "+")
  dist <- sqrt(dist_squared)

  # represent the centre location in the spectral domain
  adjust <- matrix(0, nrow = M, ncol = N)
  adjust[M / 2, N / 2] <- 1
  fft_adjust <- (stats::fft(adjust) * M * N)

  # find the locations of the non-missing cells in the original raster
  # not_missing_idx <- which(!is.na(values(grid_raster)))

  # find their locations in the extended torus grid (using the top left corner)
  extract_coords <- cbind(rowFromCell(grid_raster, cells(grid_raster)),
                          colFromCell(grid_raster, cells(grid_raster)))

  # return this info as a list
  list(dim_all = c(M, N),
       dim_sub = c(m, n),
       dist = dist,
       fft_adjust = fft_adjust,
       extract_coords = extract_coords,
       grid_raster = grid_raster)

}

# find the smallest integer, larger than `n`, that is an integer-power of 2
ceiling2 <- function(n) {
  2 ^ ceiling(log2(n))
}

# given a number of cells `n`, and a cell width `diff`, return a vector giving
# the squared distance from the centrepoint of the vector
centred_squared_distance <- function(n, diff) {
  seq <- seq_len(n) * diff
  mid <- diff * n / 2
  (seq - mid) ^ 2
}
