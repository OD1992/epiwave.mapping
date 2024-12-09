#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param prevalence_rast
#' @param population_rast
#' @param n_surveys
#' @param n_samples
#' @return
#' @author Nick Golding
#' @noRd
# simulate some prevalence surveys. one_prevalence_rast should be a single layer
# of prevalence, giving the timepoint in which to simulate the surveys,
# population_rast should be a static raster giving the population counts,
# n_surveysnis the number of survey locations, and n_samples the number of
# individuals selected at random during each survey
sim_prev_surveys <- function(one_prevalence_rast,
                             population_rast,
                             n_surveys = 30,
                             n_samples = 1000) {

  # simulate some survey locations, with a bias away from population
  prev_survey_locs <- sim_clusters(log1p(population_rast), n = n_surveys)

  # look up the prevalences
  prevs <- terra::extract(one_prevalence_rast, prev_survey_locs)[, 1]

  # simulate the prevalence survey data
  dplyr::tibble(
    n_sampled = rep(n_samples, n_surveys),
    n_positive = rbinom(n_surveys, size = n_sampled, prob = prevs)
  ) %>%
    dplyr::bind_cols(
      prev_survey_locs
    ) %>%
    as.data.frame()

}

