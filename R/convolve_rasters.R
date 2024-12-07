#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param rast
#' @param kernel
#' @param diffs
#' @return
#' @author Nick Golding
#' @export
# Given a SpatRaster `rast`, with layers giving rasters for contiguous times of
# a timeseries, and a function `kernel` taking in an integer value for a time
# difference (in the same discrete units as the raster layers), apply the kernel
# to convolve the raster through time. If kernel was created using
# epiwave.mapping:: transform_convolutional_kernel(), the maximum timeseries
# differences to convolve should be automatically detected. If not, they can be
# passed via the `diffs` argument. It is the user's responsibility to check that
# the kernel and diffs are correct in this case.
convolve_rasters <- function(rast, kernel, diffs = NULL) {

  # maybe use user-provided diff levels
  if (is.null(diffs)) {
    e <- environment(q)
    diffs <- e$table$timeperiod_diff
  }

  # check this makes sense
  if(!identical(sort(diffs), seq_along(diffs) - 1)) {
    stop("diffs must be a continuous increasing sequence of integers, starting at 0")
  }
  n_diffs <- length(diffs)

  # blank raster for padding earlier versions
  blank <- rast[[1]] * 0

  # amke a list of weighted raster stacks, to sum together
  weighted_list <- list()

  # for each diff level
  for (i in seq_len(n_diffs)) {
    # compute the weight
    diff <- diffs[i]
    weight <- q(diff)

    # maybe pad the start of the timeseries
    if (diff > 0) {
      rast_offset <- c(
        rep(blank, diff),
        rast[[-diff]]
      )
    } else {
      rast_offset <- rast
    }

    # apply the weights for this diff
    weighted_list[[i]] <- rast_offset * weight
  }

  # combine the weighted rasters and return
  convolved <- do.call(`+`, weighted_list)

  convolved

}
