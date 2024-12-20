#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param template_rast
#' @param n_times
#' @return
#' @author Nick Golding
#' @noRd
# simulate a spatio-temporal random field (Matern x AR1) using the fields R
# package, and return as a SpatRaster
sim_epsilon <- function(template_rast,
                        n_times,
                        sigma = 0.6,
                        phi_space = 0.8,
                        theta_time = 0.9) {

  # fields' BCB setup
  cells <- terra::cells(template_rast)
  ext <- terra::ext(template_rast)
  dims <- dim(template_rast)
  y_grid <- seq(ext[3], ext[4], length.out = dims[1])
  x_grid <- seq(ext[1], ext[2], length.out = dims[2])

  # convert phi from the greta definition to fields definition
  phi_space_fields <- phi_space / (0.5 * sqrt(8 * 2.5))
  cov_setup <- fields::matern.image.cov(
    grid = list(x = x_grid,
                y = y_grid),
    setup = TRUE,
    aRange = phi_space_fields,
    smoothness = 2.5 # matern 5/2
  )

  # empty matrix for random white noise values
  Y <- matrix(0, dims[1], dims[2])

  # now loop through times, doing AR1 process
  eps_list <- replicate(n_times,
                        template_rast * 0,
                        simplify = FALSE)
  eps <- do.call(c, eps_list)
  names(eps) <- paste0("time_", seq_len(n_times))

  # loop through the time points innovating the spatial process
  for (t in seq_len(n_times)) {
    # new white noise
    Y[] <- rnorm(length(Y))
    # temporal correlation - can't work out why we need to divide by 30 to get
    # correct marginal variance, but apparently we do
    sim <- (sigma / 30) * fields::matern.image.cov(Y = Y, cov.obj = cov_setup)

    # in the first time point, just use it as a raw value
    if (t == 1) {
      eps[[t]][cells] <- sim[cells]
    } else {
      eps[[t]][cells] <- theta_time * eps[[t - 1]][cells] + sim[cells]
    }

  }

  eps

}

