#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param grid_raster
#' @return
#' @author Nick Golding
#' @export
define_bcb_setup <- function(grid_raster) {

  stop("TO DO")

  # check this raster is projected

  # get the dimensions of the raster

  # double them to get the toroidal dimensions

  # build the BCB objects: basis, number of elements

  # get an index to the cells we care about in the torus

  # return a list of these things
  list(
    n = basis_size,
    basis = basis,
    dim = basis_dimensions,
    original_dim = original_dimensions)

}
