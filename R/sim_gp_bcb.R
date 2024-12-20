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
sim_gp_bcb <- function(bcb_setup = bcb_setup,
                       space_kernel,
                       sigma = 1,
                       n_times = 1,
                       time_correlation = 1) {

  # given the bcb setup and number of times to simulate, define the latent white
  # noise parameters for each times
  v <- normal(0, 1, dim = c(bcb_setup$n_variables, n_times))

  # evaluate kernel on the bcb space and get simulated values for all non-NA
  # pixels in the template raster
  f <- bcb_gp(v = v,
              kernel = space_kernel,
              bcb_setup = bcb_setup)

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

# given standard normal matrix v, kernel, and bcb setup, colour the spatial
# correlation and then return only the elements of the compute grid that were in
# the original raster
bcb_gp <- function (v, kernel, bcb_setup) {
  f <- bcb_gp_colour(bcb_setup, v, kernel)
  bcb_gp_extract(bcb_setup, f)
}

# do the colouring (applying covariance structure) of the random variables
bcb_gp_colour <- function (bcb_setup, v, kernel) {

  # apply the greta.gp spatial kernel to the minimal vector of distances (one
  # per gridcell on the torus, so ~4x the number of gridcells in a square region
  # of interest) representing the BCB setup, to get the corresponding
  # covariances

  # since the kernel is expecting locations, but the bcb_setup provides
  # distances from the centroid, we need to flatten the grid of distances,
  # compute covariance from 0, and then reshape
  covar_vec <- kernel(bcb_setup$dist_vec, 0)

  spectral_colour(v, covar_vec, bcb_setup)

}

# given two greta arrays, with matching dimensions, for pointwise covariances
# and standard normal random variables, 'colour' them (apply
# correlation/covariance) and return the greta array of coloured random
# variables
spectral_colour <- function (v, covar_vec, bcb_setup) {

  # check inputs are both greta arrays
  if (!inherits(covar_vec, "greta_array") |
      !inherits(v, "greta_array")) {
    stop ("'covar_vec' and 'v' must both be greta arrays",
          call. = FALSE)
  }

  # check they have the same first dimensions
  if (!(dim(covar_vec)[1] == dim(v)[1])) {
    stop ("'covar_vec' and 'v' must have the same first dimensions, ",
          "but 'covar_vec' had first dimension ",
          dim(covar_vec)[1],
          " and 'v' had first dimension ",
          dim(v)[1],
          call. = FALSE)
  }

  # create and return an operation greta array for the result
  op <- greta::.internals$nodes$constructors$op

  op("spectral_colour",
     v, covar_vec,
     operation_args = list(bcb_setup = bcb_setup),
     tf_operation = "tf_spectral_colour")

}

# tensorflow code for the spectral colouring operation
tf_spectral_colour <- function (v, covar_vec, bcb_setup) {

  tf <- tensorflow::tf

  # load the fft adjustment vector, and reshape to align with dimensions
  # of other arrays (needed to avoid a crash at complex division step!)
  fft_adjust_vec <- bcb_setup$fft_adjust_vec
  dim(fft_adjust_vec) <- c(1, length(fft_adjust_vec), 1)

  # find the complex type corresponding to the user's float type
  float_type <- options()$greta_tf_float
  complex_type <- switch(float_type,
                         float32 = tf$complex64,
                         float64 = tf$complex128)

  # cast the covariance, random normals, and fft adjustment,
  # to complex (only type fft accepts)
  covar_vec <- tf$cast(covar_vec, complex_type)
  v <- tf$cast(v, complex_type)
  fft_adjust_vec <- tf$cast(fft_adjust_vec, complex_type)

  # TF's fft operates on the innermost (last) dimension of the tensor, but we
  # want it on the second dimension, so move batch dimensions around to do
  # calculations
  perm <- c(0L, 2L, 1L)
  reverse_perm <- c(0L, 2L, 1L)
  covar_vec <- tf$transpose(covar_vec, perm)
  v <- tf$transpose(v, perm)
  fft_adjust_vec <- tf$transpose(fft_adjust_vec, perm)

  # cast covariance and colourless normals into spectral domain
  covar_vec_fft <- tf$signal$fft(covar_vec) / fft_adjust_vec
  v_fft <- tf$signal$fft(v)

  # get coloured function in spectral domain
  f_fft <- tf$math$sqrt(covar_vec_fft) * v_fft

  # bring coloured function back from the spectral domain
  f <- tf$math$real(tf$signal$ifft(f_fft))

  # rescale by the dimension size
  root_nelem <- greta:::fl(sqrt(dim(v)[1]))
  f <- f / root_nelem

  # cast the float type back
  f <- tf$cast(f, options()$greta_tf_float)

  # permute it again, to put the spatial dimension second again, and return
  f <- tf$transpose(f, reverse_perm)
  f

}

# given the bcb setup information and the full coloured random variables,
# extract only those the user wants
bcb_gp_extract <- function (bcb_setup, f) {

    f[bcb_setup$extract_index_vec, ]

}


