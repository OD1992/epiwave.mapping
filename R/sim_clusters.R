#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param nameme1
#' @param n
#' @return
#' @author Nick Golding
#' @noRd
sim_clusters <- function(population_rast = population_rast,
                         n = n,
                         n_initial = 1e4) {

  pts <- terra::spatSample(
    population_rast ^ 2,
    size = 1e4,
    method = "weights",
    replace = TRUE,
    na.rm = TRUE,
    xy = TRUE,
    values = FALSE
  )

  # simulate cluster locations with k-means
  kmn <- kmeans(pts,
                centers = n)
  kmn$centers
}
