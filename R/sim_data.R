covariates_rast <- bioclim_kenya[[1:5]]
population_rast <- pop_kenya
n_health_facilities <- 100
years <- 2020:2024
n_prev_surveys <- 30
n_times <- length(years) * 12


# simulate infection incidence data

# cells <- cells(covariates_rast)
# all_coords <- xyFromCell(covariates_rast, cells)
# X <- extract(covariates_rast_scaled, cells)
alpha <- rnorm(1, -8, 1)
beta <- rnorm(nlyr(covariates_rast))
epsilon <- sim_epsilon(template_rast = population_rast,
                       n_times = n_times)

# scale the covariates
covariates_rast_scaled <- scale(covariates_rast)

# create some fixed effects raster part
fixef <- alpha + sum(covariates_rast_scaled * beta)

# create eta
eta <- fixef + epsilon

# create new infections
new_infections <- population_rast * exp(eta)

# simulate infection prevalence on rasters

# given a function of test positivity over time since infection (days), convert
# this into the convolution function q, for months

# for each day post-infection this is the (entirely fake) probability of testing
# positive on our hypothetical test, if tested that day. Hard-coded to be 0 for
# negative times or beyond 30 days
max_infection_detectability <- 30
q_daily <- function(days) {
  up <- plogis((days * 2) - 5)
  down <- 1 - plogis(days / 2 - 10)
  q <- up  * down
  q[days < 0] <- 0
  q[days > 30] <- 0
  q
}



# create a discrete kernel, in days, for the delay from infection to case
# reporting, then convert it to a monthly kernel

# now make it a discrete daily distribution
pi_daily <- function(day_difference) {
  # use a truncated lognormal (truncated 0-14)
  mu <- 1.2
  sigma <- 0.5
  upper_limit <- 14
  day_difference_old <- day_difference
  day_difference <- round(day_difference)
  if (!max(abs(day_difference - day_difference_old) < 1e-4)) {
    warning("day_difference was non-integer, using rounded values")
  }

  upper <- plnorm(day_difference + 1, mu, sigma)
  lower <- plnorm(day_difference, mu, sigma)
  dens <- upper - lower

  # correct the densities
  norm <- plnorm(upper_limit, mu, sigma)
  res <- dens / norm
  res[day_difference > upper_limit] <- 0
  res[day_difference < 0] <- 0
  res
}

plot(pi_continuous, xlim = c(0, 21))

# check I did this right
# integrate(pi_continuous, 0, 14)

# make the kernel for convolving this on our monthly timestep
month_length_days <- 365.25/12

q <- transform_convolution_kernel(kernel_daily = q_daily,
                                  max_diff_days = max_infection_detectability,
                                  timeperiod_days = month_length_days)

# now make this work as a convolution on the discrete timeperiods
pi <- transform_convolution_kernel(kernel_daily = pi_daily,
                                   max_diff_days = 14,
                                   timeperiod_days = month_length_days)
e <- environment(pi)
sum(e$q_table$q)
x <- seq(0, 16, length.out = 1000)
plot(pi_daily(x) ~x, type = "l")
pi_daily(14.1)
space_sims <- replicate(
  n_times,
  {
    matern.image.cov(Y = Y, cov.obj = cov_setup)
  },
  simplify = FALSE
)

# loop through these

image.plot(space_sims[[3]])
# simulate some prevalence surveys
sim_prev_surveys <- function(prevalence_rast,
                             population_rast,
                             n_surveys,
                             n_samples)

prev_survey_locs <- sim_clusters(log1p(population_rast), n = 30)

# simulate some health facilities
health_facilities <- sim_health_facilities(
  population_rast = population_rast,
  n_health_facilities = n_health_facilities)

plot(population_rast)
points(prev_survey_locs)

plot(population_rast)
points(health_facilities$health_facility_loc)
