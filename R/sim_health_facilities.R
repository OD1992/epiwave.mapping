#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param population_rast
#' @param n_health_facilities
#' @return
#' @author Nick Golding
#' @noRd
# simulate the health facility locations and their catchments
sim_health_facilities <- function(population_rast, n_health_facilities) {

  # simulate health facility locations biased by population
  locs <- sim_clusters(population_rast, n = n_health_facilities)

  # find the distance from each pixel to each health facility
  all_coords <- xyFromCell(population_rast, cells(population_rast))
  d <- distance(all_coords,
                locs,
                lonlat = TRUE)

  # model the relative probability of travel there
  prob <- exp(-d / 1e5)
  # prob[prob > quantile(prob, 0.5)] <- 0
  weights <- sweep(prob, 1, rowSums(prob), FUN = "/")

  # make rasters of these
  catchments_list <- list()
  for (h in seq_len(n_health_facilities)) {
    catchment <- population_rast * 0
    names(catchment) <- paste0("catchment_", h)
    catchment[cells(catchment)] <- weights[, h]
    catchments_list[[h]] <- catchment
  }
  catchments <- do.call(c, catchments_list)

  # return all these things
  list(
    health_facility_loc = locs,
    health_facility_weights = weights,
    catchments_rast = catchments)

}
