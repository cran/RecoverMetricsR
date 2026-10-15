#' @title Assess Recovery-Time Threshold Sensitivity
#'
#' @description
#' Computes time-to-recovery (TTR) for multiple recovery thresholds and
#' quantifies the sensitivity of TTR to the selected threshold.
#'
#' @param raster_stack A \code{terra::SpatRaster} containing yearly spectral
#'   index layers.
#' @param years Numeric vector of years corresponding to the raster layers.
#' @param event_year Numeric disturbance year.
#' @param thresholds Numeric vector of recovery proportions. Values must be
#'   in the interval \code{(0, 1]}.
#' @param pre_window Number of observations before the disturbance used to
#'   calculate the pre-disturbance baseline.
#' @param persistence Number of consecutive post-disturbance observations
#'   that must meet or exceed the recovery threshold.
#' @param mask_layer Optional spatial mask. Can be a \code{terra::SpatRaster},
#'   \code{terra::SpatVector}, or \code{sf} object, passed to
#'   \code{\link{compute_recovery_time}}.
#' @param return_band Logical. If \code{TRUE}, also returns the range of TTR
#'   values obtained across the specified thresholds. Default is
#'   \code{TRUE}.
#'
#' @details
#' Recovery time is computed independently for each threshold using
#' \code{\link{compute_recovery_time}}. The threshold range is calculated
#' as the difference between the maximum and minimum TTR across thresholds:
#' \deqn{
#' TTR_{range} = \max(TTR_p) - \min(TTR_p)
#' }
#' where \code{p} denotes the selected recovery proportions.
#'
#' The resulting \code{TTR_Threshold_Range} represents sensitivity of the
#' estimated recovery time to the choice of recovery threshold. It should
#' not be interpreted as a statistical confidence interval or a formal
#' measure of uncertainty.
#'
#' Threshold-based recovery times have been used in post-fire vegetation
#' recovery studies, including definitions based on a specified proportion
#' of an undisturbed reference condition (Chompuchan and Lin, 2017).
#'
#' @return
#' A \code{terra::SpatRaster} containing one layer for each recovery
#' threshold. Layer names are formatted as \code{"TTR_70"}, \code{"TTR_80"},
#' etc. If \code{return_band = TRUE} and more than one threshold is supplied,
#' an additional layer named \code{"TTR_Threshold_Range"} is returned.
#'
#' @references
#' Chompuchan, C., & Lin, C.-Y. (2017). Assessment of forest recovery
#' at the deciduous forest in a complex terrain using Landsat time
#' series. \emph{Ecological Indicators}, 79, 249--259.
#' \doi{10.1016/j.ecolind.2017.04.038}
#'
#' @examples
#' library(terra)
#'
#' # Small toy example -- runs in well under 5 seconds
#' r <- rast(nrows = 2, ncols = 2, nlyrs = 6)
#' values(r) <- rep(c(0.60, 0.20, 0.35, 0.42, 0.50, 0.56), each = 4)
#' names(r) <- paste0("NBR_", 2020:2025)
#'
#' result <- compute_recovery_time_uncertainty(
#'   raster_stack = r,
#'   years = 2020:2025,
#'   event_year = 2021,
#'   thresholds = c(0.7, 0.8, 0.9),
#'   pre_window = 1,
#'   persistence = 1
#' )
#'
#' result
#'
#' \dontshow{
#' # Internal checks (not shown to users)
#' stopifnot(inherits(result, "SpatRaster"))
#' stopifnot(terra::nlyr(result) == 4L)
#' stopifnot("TTR_Threshold_Range" %in% names(result))
#'
#' # return_band = FALSE suppresses the range layer
#' result_nb <- compute_recovery_time_uncertainty(
#'   raster_stack = r,
#'   years = 2020:2025,
#'   event_year = 2021,
#'   thresholds = c(0.7, 0.8, 0.9),
#'   pre_window = 1,
#'   persistence = 1,
#'   return_band = FALSE
#' )
#' stopifnot(terra::nlyr(result_nb) == 3L)
#' stopifnot(!"TTR_Threshold_Range" %in% names(result_nb))
#' }
#'
#' \donttest{
#' # Without the sensitivity band, and with a spatial mask
#' library(sf)
#'
#' p <- st_sfc(
#'   st_polygon(list(rbind(
#'     c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)
#'   ))),
#'   crs = crs(r)
#' )
#' p <- st_sf(id = 1, geometry = p)
#'
#' result_masked <- compute_recovery_time_uncertainty(
#'   raster_stack = r,
#'   years = 2020:2025,
#'   event_year = 2021,
#'   thresholds = c(0.6, 0.7, 0.8, 0.9),
#'   pre_window = 1,
#'   persistence = 2,
#'   mask_layer = p,
#'   return_band = FALSE
#' )
#'
#' result_masked
#' }
#'
#' @importFrom terra nlyr rast app
#' @export
compute_recovery_time_uncertainty <- function(
  raster_stack,
  years,
  event_year,
  thresholds = c(0.7, 0.8, 0.9),
  pre_window = 2,
  persistence = 1,
  mask_layer = NULL,
  return_band = TRUE
) {

  # ---- Input validation ----

  if (!inherits(raster_stack, "SpatRaster")) {
    stop("'raster_stack' must be a terra::SpatRaster")
  }

  if (!is.numeric(years) || anyNA(years)) {
    stop("'years' must be a numeric vector without NA values")
  }

  if (length(years) != terra::nlyr(raster_stack)) {
    stop("Length of 'years' must match the number of raster layers")
  }

  if (anyDuplicated(years)) {
    stop("'years' must contain unique values")
  }

  if (!is.numeric(event_year) ||
      length(event_year) != 1L ||
      is.na(event_year)) {
    stop("'event_year' must be a single numeric value")
  }

  if (!event_year %in% years) {
    stop("'event_year' must be one of the values in 'years'")
  }

  if (!is.numeric(thresholds) ||
      length(thresholds) < 1L ||
      anyNA(thresholds)) {
    stop("'thresholds' must be a non-empty numeric vector")
  }

  if (any(thresholds <= 0 | thresholds > 1)) {
    stop("'thresholds' must contain values in the interval (0, 1]")
  }

  if (anyDuplicated(thresholds)) {
    stop("'thresholds' must contain unique values")
  }

  if (!is.numeric(pre_window) ||
      length(pre_window) != 1L ||
      is.na(pre_window) ||
      pre_window < 1) {
    stop("'pre_window' must be a positive integer")
  }

  if (pre_window != as.integer(pre_window)) {
    stop("'pre_window' must be an integer")
  }

  if (!is.numeric(persistence) ||
      length(persistence) != 1L ||
      is.na(persistence) ||
      persistence < 1) {
    stop("'persistence' must be a positive integer")
  }

  if (persistence != as.integer(persistence)) {
    stop("'persistence' must be an integer")
  }

  if (!is.logical(return_band) ||
      length(return_band) != 1L ||
      is.na(return_band)) {
    stop("'return_band' must be TRUE or FALSE")
  }

  if (isTRUE(return_band) && length(thresholds) == 1L) {
    warning(
      "'return_band = TRUE' requires at least two thresholds; ",
      "the sensitivity band will not be computed.",
      call. = FALSE
    )
  }

  # ---- Compute TTR for each threshold ----

  ttr_list <- lapply(thresholds, function(p) {

    compute_recovery_time(
      raster_stack = raster_stack,
      years = years,
      event_year = event_year,
      proportion = p,
      pre_window = as.integer(pre_window),
      persistence = as.integer(persistence),
      mask_layer = mask_layer,
      return_year = FALSE
    )
  })

  # ---- Layer names for threshold-specific TTR ----

  pct <- thresholds * 100

  labels <- ifelse(
    pct == floor(pct),
    formatC(pct, format = "d"),
    formatC(pct, format = "fg", digits = 3)
  )

  names(ttr_list) <- paste0("TTR_", labels)

  # ---- Combine threshold-specific TTR layers ----

  ttr_stack <- terra::rast(ttr_list)

  # ---- Threshold sensitivity range ----

  if (isTRUE(return_band) && length(thresholds) > 1L) {

    ttr_band <- terra::app(
      ttr_stack,
      fun = function(x) {

        if (all(is.na(x))) {
          return(NA_real_)
        }

        max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
      }
    )

    names(ttr_band) <- "TTR_Threshold_Range"

    ttr_stack <- c(ttr_stack, ttr_band)
  }

  ttr_stack
}
