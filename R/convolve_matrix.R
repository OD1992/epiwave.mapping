#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param timeperiod_new_infections_prev_locs
#' @param kernel
#' @return
#' @author Nick Golding
#' @export

# Given a `matrix` with rows timeseries over different locations and columns
# giving values in a series of contiguous and increasing timeperiods, and a
# function `kernel` taking in an integer value for a time difference (in the
# same time period units as the matrix columns), apply the kernel to convolve
# the raster through time. If kernel was created using epiwave.mapping::
# transform_convolutional_kernel(), the maximum timeseries differences to
# convolve should be automatically detected. If not, they can be passed via the
# `diffs` argument. It is the user's responsibility to check that the kernel and
# diffs are correct in this case.
convolve_matrix <- function(matrix,
                            kernel,
                            diffs = NULL) {

  # maybe use user-provided diff levels
  if (is.null(diffs)) {
    e <- environment(kernel)
    diffs <- e$table$timeperiod_diff
  }

  # check this makes sense
  if(!identical(sort(diffs), seq_along(diffs) - 1)) {
    stop("diffs must be a continuous increasing sequence of integers, starting at 0")
  }
  n_diffs <- length(diffs)

  # blank matrix column for padding earlier versions
  blank <- matrix[, 1] * 0

  # make a list of weighted matrices, to sum together
  weighted_list <- list()

  # for each diff level
  for (i in seq_len(n_diffs)) {
    # compute the weight
    diff <- diffs[i]
    weight <- kernel(diff)

    # maybe pad the start of the timeseries
    if (diff > 0) {
      diff_seq <- 1:diff
      matrix_offset <- cbind(
        rep(blank, diff),
        matrix[, -(1:diff)]
      )
    } else {
      matrix_offset <- matrix
    }

    # apply the weights for this diff
    weighted_list[[i]] <- matrix_offset * weight
  }

  # combine the weighted rasters and return
  convolved <- do.call(`+`, weighted_list)

  convolved

}
