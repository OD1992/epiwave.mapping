covariates_rast <- bioclim_kenya[[1:5]]
population_rast <- pop_kenya
n_health_facilities <- 100
years <- 2020:2024
n_prev_surveys <- 30
n_times <- length(years) * 12

# length of a month in days
month_length_days <- 365.25/12

# simulate infection incidence data
set.seed(1)
alpha <- rnorm(1, log(1e-4), 1)
beta <- rnorm(nlyr(covariates_rast))
gamma <- rbeta(1, 5, 10)
epsilon <- sim_epsilon(template_rast = population_rast,
                       n_times = n_times)

# scale the covariates
covariates_rast_scaled <- scale(covariates_rast)

# create some fixed effects raster part
fixef <- alpha + sum(covariates_rast_scaled * beta)

# create eta
eta <- fixef + epsilon

# get the infection incidence
infection_incidence_daily <- exp(eta)
# create raster of number of new infections per timeperiod/pixel
new_infections <- month_length_days * population_rast * infection_incidence_daily
names(new_infections) <- names(epsilon)

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

# transform to a monthly timestep
q <- transform_convolution_kernel(kernel_daily = q_daily,
                                  max_diff_days = max_infection_detectability,
                                  timeperiod_days = month_length_days)

# do convolution of rasters, and apply the rest of the maths
new_infections_daily <- new_infections / month_length_days
detectable_infections_daily <- convolve_rasters(rast = new_infections_daily,
                                                kernel = q)
prevalence <- detectable_infections_daily / population_rast

# now do the same for the pixel map of clinical case counts
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

# now make this work as a convolution on the discrete timeperiods
pi <- transform_convolution_kernel(kernel_daily = pi_daily,
                                   max_diff_days = 14,
                                   timeperiod_days = month_length_days)

# convolve the incidence, and multiply by the probability of reporting
reporting_convolved_infections_pixel <- convolve_rasters(rast = new_infections,
                                                        kernel = pi)
clinical_cases_pixel <- reporting_convolved_infections_pixel * gamma

prev_surveys <- sim_prev_surveys(one_prevalence_rast = prevalence[[16]],
                                 population_rast = population_rast,
                                 n_surveys = 30,
                                 n_samples = 1000)



# simulate some health facilities
health_facilities <- sim_health_facilities(
  population_rast = population_rast,
  n_health_facilities = n_health_facilities)

plot(population_rast)
points(prev_survey_locs)

plot(population_rast)
points(health_facilities$health_facility_loc)
