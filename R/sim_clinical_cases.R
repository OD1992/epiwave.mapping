#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param clinical_cases_rast
#' @param health_facilities
#' @return
#' @author Nick Golding
#' @export

# Given a temporally-varying raster 'clinical_cases_rast' of the number of
# clinical cases per pixel and timeperiod, and an object 'health_facilities'
# (created by sim_health_facilities()) describing the health facilities in the
# region, and a scalar fraction 'missingness' giving the proportion of case
# counts that are missing (set to NA), compute and return timeseries for each
# health facility of the number of clinical cases observed over time
sim_clinical_cases <- function(clinical_cases_rast,
                               health_facilities,
                               missingness = 0) {

  # pull out the expected clinical cases for each cell
  cell_ids <- terra::cells(clinical_cases_rast)
  clinical_cases_pixel <- terra::extract(clinical_cases_rast,
                                         cell_ids)

  # matrix-multiply by the health facility weights to get the expected counts
  # per health facility
  t_clinical_cases_hf <- t(clinical_cases_pixel) %*% health_facilities$health_facility_weights
  clinical_cases_hf <- t(t_clinical_cases_hf)
  dim(clinical_cases_hf)

  # pivot these into long form and add Poisson noise and missingness
  clinical_cases_hf %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(
      health_facility = seq_len(nrow(health_facilities$health_facility_loc)),
      .before = everything()
    ) %>%
    tidyr::pivot_longer(
      cols = starts_with("time_"),
      names_to = "time",
      names_prefix = "time_",
      values_to = "expected_cases"
    ) %>%
    dplyr::mutate(
      cases = rpois(dplyr::n(), expected_cases),
      cases = dplyr::case_when(
        rbinom(dplyr::n(), 1, missingness) == 1 ~ NA,
        .default = cases
      )
      # set some to NAs
    ) %>%
    dplyr::select(
      # drop the truth
      -expected_cases
    ) %>%
    as.data.frame()

}
