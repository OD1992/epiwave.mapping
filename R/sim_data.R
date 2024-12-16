#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param covariates_rast
#' @param population_rast
#' @param years
#' @param n_health_facilities
#' @param n_prev_surveys
#' @param case_missingness
#' @param prev_survey_time
#' @return
#' @author Nick Golding
#' @export

# Given a multi-layer SpatRaster 'covariates_rast' of temporally-static
# covariate values, a single-layer SpatRaster of 'population_rast' of temporally
# static population distribution, a vector of 'years' to consider, and the
# number of health facilities ('n_health_facilities') and prevalence surveys
# locations ('n_prev_surveys'), the proportion of month/healthfacility
# combinations for which clinical case count data are missing
# 'case_missingness', and the time period in which the prevalence survey
# occurred 'prevalence_survey_time' (in months since the start of the first
# year, by default the midpoint of 'years'), simulate health facilities,
# health-facility level clinical case counts, and infection prevalence survey
# locations and results, and return these along with the true parameters, latent
# quantities of interest and parameters of the surveillance processes.
sim_data <- function(covariates_rast,
                     population_rast,
                     years = 2020:2024,
                     n_health_facilities = 100,
                     n_prev_surveys = 30,
                     case_missingness = 0.15,
                     prev_survey_time = round(length(years) * 12 / 2)) {

  # use a monthly timestep, so create times for these years
  n_times <- length(years) * 12

  # length of a month in days
  month_length_days <- 365.25/12

  # random parameters
  alpha <- rnorm(1, log(1e-4), 1)
  beta <- rnorm(terra::nlyr(covariates_rast))
  r <- rbeta(1, 10, 5)
  sigma <- rbeta(1, 12, 8)
  phi <- rbeta(1, 16, 4)
  theta <- rbeta(1, 18, 2)

  # simulate infection incidence data
  epsilon <- sim_epsilon(template_rast = population_rast,
                         n_times = n_times,
                         sigma = sigma,
                         phi_space = phi,
                         theta_time = theta)

  # scale the covariates
  covariates_rast_scaled <- terra::scale(covariates_rast)

  # create the fixed effects part as a raster
  fixef <- alpha + sum(covariates_rast_scaled * beta)

  # get linear predictor and transform to daily infection incidence rate
  eta <- fixef + epsilon
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

  max_case_delay <- 14

  # now do the same for the pixel map of clinical case counts
  pi_daily <- function(day_difference) {
    upper_limit <- max_case_delay
    # use a truncated lognormal (truncated to 0-upper_limit)
    mu <- 1.2
    sigma <- 0.5
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

  # transform to a monthly timestep
  q <- transform_convolution_kernel(kernel_daily = q_daily,
                                    max_diff_days = max_infection_detectability,
                                    timeperiod_days = month_length_days)

  # create a function g for the prevalence-incidence relationship
  g <- function(prevalence) {
    (1/365) * (1 - exp(-prevalence * 10))
  }

  # infections to reported clinical cases, as a function of infection incidence
  q_star <- sum(q_daily(seq(0, max_infection_detectability)))
  compute_gamma <- function(infection_incidence_daily) {
    r * g(q_star * infection_incidence_daily) / infection_incidence_daily
  }

  # do convolution of rasters, and apply the rest of the maths
  new_infections_daily <- new_infections / month_length_days
  detectable_infections_daily <- convolve_rasters(rast = new_infections_daily,
                                                  kernel = q)
  prevalence <- detectable_infections_daily / population_rast
  names(prevalence) <- names(epsilon)

  # transform delay probability distribution to discrete timeperiods
  pi <- transform_convolution_kernel(kernel_daily = pi_daily,
                                     max_diff_days = max_case_delay,
                                     timeperiod_days = month_length_days)

  # compute the number of new clinical cases by infection date convolve the
  # incidence, and multiply by the probability of reporting
  gamma_rast <- compute_gamma(infection_incidence_daily)
  new_clinical_cases <- gamma_rast * new_infections
  reported_clinical_cases_pixel <- convolve_rasters(rast = new_clinical_cases,
                                                           kernel = pi)
  names(reported_clinical_cases_pixel) <- names(epsilon)

  # back out the true number of clinical cases (since reporting rate r is
  # wrapped into gamma)
  clinical_cases_pixel <- reported_clinical_cases_pixel / r

  # simulate data collection

  # simulate a prevalence survey
  prev_slice <- prevalence[[prev_survey_time]]
  prevalence_surveys <- sim_prev_surveys(one_prevalence_rast = prev_slice,
                                         population_rast = population_rast,
                                         n_surveys = n_prev_surveys,
                                         n_samples = 1000)

  # add on the time period
  prevalence_surveys <- prevalence_surveys %>%
    dplyr::mutate(
      time = prev_survey_time
    )

  # simulate some health facilities at which to observe case counts
  health_facilities <- sim_health_facilities(
    population_rast = population_rast,
    n_health_facilities = n_health_facilities)

  # aggregate expected case counts by these facilities and simulate case data
  clinical_cases <- sim_clinical_cases(
    clinical_cases_rast = reported_clinical_cases_pixel,
    health_facilities = health_facilities,
    missingness = case_missingness)

  # return all these objects as a structured list
  list(
    epi_data = list(
      # prevalence survey data, locations, and time
      prevalence_surveys = prevalence_surveys,
      # clinical case counts
      clinical_cases = clinical_cases
    ),
    surveillance_information = list(
      # health facility locations and catchment probabilities
      health_facilities = health_facilities,
      # daily probability of positive result in prevalence diagnostic
      prev_detectability_daily_fun = q_daily,
      # maximum detectability in days
      prev_detectability_max_days =  max_infection_detectability,
      # probability mass function of case reporting delay distribution in days
      case_delay_distribution_daily_fun = pi_daily,
      # maximum delay in case reporting in days
      case_delay_max_days = max_case_delay,
      # clinical incidence as a function of prevalence
      prev_inc_function = g
    ),
    # objects contianing the truth
    truth = list(
      # model fixed hyperparameters
      parameters = list(
        alpha = alpha,
        beta = beta,
        gamma = gamma,
        sigma = sigma,
        phi = phi,
        theta = theta
      ),
      rasters = list(
        # space-time random effect, as a raster
        epsilon_rast = epsilon,
        # unobserved quantities, as rasters
        # daily incidence rate of new infections
        infection_incidence_rast = infection_incidence_daily,
        # population prevalence of infection
        prevalence_rast = prevalence,
        # number of clinical cases per timeperiod
        clinical_cases_rast = clinical_cases_pixel
      )
    )
  )

}

