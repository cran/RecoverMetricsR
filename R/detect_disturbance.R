#' @title Detect Disturbance Year from Spectral Time Series
#'
#' @description
#' Detects a disturbance year from a raster spectral time series using
#' either a year-on-year difference rule or a LandTrendr-style local
#' OLS slope approach. A user-provided disturbance year can also be
#' supplied, in which case automatic detection is skipped.
#'
#' @param raster_stack A \code{terra::SpatRaster} containing temporal
#'   spectral index layers.
#' @param years Numeric vector giving the year corresponding to each
#'   raster layer. Its length must equal the number of raster layers.
#' @param disturbance_year Optional numeric value specifying a known
#'   disturbance year. If supplied, automatic detection is skipped.
#' @param mask_layer Optional spatial mask. Can be a
#'   \code{terra::SpatRaster}, \code{terra::SpatVector}, or \code{sf}
#'   polygon identifying the disturbed area. The mask is used only when
#'   \code{disturbance_year} is \code{NULL}. Raster masks are resampled
#'   to \code{raster_stack} if their CRS, resolution, or extent differ.
#'   Vector masks are reprojected to the raster CRS if needed.
#' @param detect Character string specifying the automatic detection
#'   method. Must be either \code{"difference"} or
#'   \code{"landtrendr"}. The argument is ignored when
#'   \code{disturbance_year} is supplied.
#'
#' @details
#' When \code{disturbance_year} is \code{NULL}, the function calculates
#' a spatial mean time series from \code{raster_stack}, optionally
#' restricted to \code{mask_layer}.
#'
#' Two automatic detection methods are available:
#'
#' \itemize{
#'   \item \code{"difference"}: identifies the year following the
#'     largest negative year-on-year change in the spatial mean
#'     spectral index.
#'
#'   \item \code{"landtrendr"}: identifies the year following the
#'     most negative local OLS slope calculated over a short moving
#'     window. The window contains up to four consecutive observations
#'     (three-year span), following a LandTrendr-style local trend
#'     approach. This option is not an implementation of the
#'     LandTrendr segmentation algorithm.
#' }
#'
#' If \code{mask_layer} is supplied and masking produces only
#' \code{NA} values, a warning is issued and the original unmasked
#' raster is used for detection.
#'
#' The detected year corresponds to the observation immediately
#' following the largest negative change or the most negative local
#' slope.
#'
#' The LandTrendr-style slope is computed against positions in the
#' moving window, not calendar years. For uniformly spaced annual data
#' (the intended use), the two are equivalent.
#'
#' The user-provided \code{disturbance_year} is returned as-is. Callers
#' are responsible for ensuring that it matches a value in \code{years}
#' if subsequent functions depend on that.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{event_year}{Detected or user-provided disturbance year.}
#'   \item{mask}{The supplied disturbance mask, if any.}
#' }
#'
#' @references
#' Kennedy, R. E., Yang, Z., & Cohen, W. B. (2010). Detecting trends
#' in forest disturbance and recovery using yearly Landsat time series:
#' 1. LandTrendr -- Temporal segmentation algorithms. \emph{Remote Sensing
#' of Environment}, 114(12), 2897--2910.
#' \doi{10.1016/j.rse.2010.07.008}
#'
#' @examples
#' library(terra)
#'
#' # Small synthetic time series with a clear decrease in 2021
#' r <- rast(
#'   nrows = 2,
#'   ncols = 2,
#'   nlyrs = 6
#' )
#'
#' values(r) <- rep(
#'   c(0.60, 0.20, 0.30, 0.38, 0.45, 0.50),
#'   each = 4
#' )
#'
#' names(r) <- paste0(
#'   "Index_",
#'   2020:2025
#' )
#'
#' # Automatic detection using year-on-year differences
#' det_diff <- detect_disturbance(
#'   raster_stack = r,
#'   years = 2020:2025,
#'   detect = "difference"
#' )
#'
#' det_diff
#'
#' # User-provided disturbance year
#' det_user <- detect_disturbance(
#'   raster_stack = r,
#'   years = 2020:2025,
#'   disturbance_year = 2021
#' )
#'
#' det_user
#'
#' \dontshow{
#' # Internal checks (not shown to users)
#' stopifnot(is.list(det_diff))
#' stopifnot(det_diff$event_year == 2021)
#' stopifnot(is.list(det_user))
#' stopifnot(det_user$event_year == 2021)
#' }
#'
#' \donttest{
#' # LandTrendr-style local slope detection
#' det_lt <- detect_disturbance(
#'   raster_stack = r,
#'   years = 2020:2025,
#'   detect = "landtrendr"
#' )
#'
#' det_lt
#'
#' # Detection restricted to an sf polygon
#' library(sf)
#'
#' p <- st_sfc(
#'   st_polygon(list(rbind(
#'     c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)
#'   ))),
#'   crs = crs(r)
#' )
#'
#' p <- st_sf(id = 1, geometry = p)
#'
#' det_masked <- detect_disturbance(
#'   raster_stack = r,
#'   years = 2020:2025,
#'   mask_layer = p,
#'   detect = "difference"
#' )
#'
#' det_masked
#' }
#'
#' @importFrom terra global mask nlyr project same.crs vect resample res ext
#' @importFrom stats lm coef
#' @export
detect_disturbance <- function(
  raster_stack,
  years,
  disturbance_year = NULL,
  mask_layer = NULL,
  detect = c(
    "difference",
    "landtrendr"
  )
) {

  # ------------------------------------------------------------------
  # Input validation
  # ------------------------------------------------------------------

  detect <- match.arg(detect)

  if (!inherits(raster_stack, "SpatRaster")) {
    stop("raster_stack must be a terra::SpatRaster.")
  }

  if (!is.numeric(years)) {
    stop("years must be numeric.")
  }

  if (length(years) != terra::nlyr(raster_stack)) {
    stop("Length of 'years' must match the number of raster layers.")
  }

  if (anyNA(years)) {
    stop("'years' must not contain NA values.")
  }

  if (anyDuplicated(years)) {
    stop("'years' must not contain duplicated values.")
  }

  # ------------------------------------------------------------------
  # User-provided disturbance year
  #
  # Return early. Convert an sf mask to SpatVector first so that the
  # returned `mask` element has a consistent type across both branches.
  # ------------------------------------------------------------------

  if (!is.null(disturbance_year)) {

    if (!is.numeric(disturbance_year) ||
        length(disturbance_year) != 1 ||
        is.na(disturbance_year)) {
      stop("disturbance_year must be a single non-missing numeric value.")
    }

    if (inherits(mask_layer, "sf")) {
      mask_layer <- terra::vect(mask_layer)
    }

    return(
      list(
        event_year = disturbance_year,
        mask = mask_layer
      )
    )
  }

  # ------------------------------------------------------------------
  # Optional disturbance mask
  # ------------------------------------------------------------------

  if (!is.null(mask_layer)) {

    if (inherits(mask_layer, "sf")) {
      mask_layer <- terra::vect(mask_layer)
    }

    if (!inherits(mask_layer, c("SpatRaster", "SpatVector"))) {
      stop(
        "mask_layer must be a terra::SpatRaster, ",
        "terra::SpatVector, or sf object."
      )
    }

    if (inherits(mask_layer, "SpatRaster")) {

      needs_resample <- !terra::same.crs(raster_stack, mask_layer) ||
        !isTRUE(all.equal(
          terra::res(raster_stack),
          terra::res(mask_layer)
        )) ||
        !isTRUE(all.equal(
          as.vector(terra::ext(raster_stack)),
          as.vector(terra::ext(mask_layer))
        ))

      if (needs_resample) {
        mask_layer <- terra::resample(
          mask_layer,
          raster_stack,
          method = "near"
        )
      }

    } else {

      # Reproject vector masks to the raster CRS.
      if (!terra::same.crs(raster_stack, mask_layer)) {
        mask_layer <- terra::project(mask_layer, raster_stack)
      }
    }

    raster_masked <- terra::mask(raster_stack, mask_layer)

    # Check whether the mask produced usable data.
    global_mean_masked <- terra::global(
      raster_masked,
      "mean",
      na.rm = TRUE
    )[, 1]

    if (all(is.na(global_mean_masked))) {

      warning(
        "Masked raster contains only NA values. ",
        "Using the full raster for detection.",
        call. = FALSE
      )

    } else {

      raster_stack <- raster_masked
    }
  }

  # ------------------------------------------------------------------
  # Spatial mean time series
  # ------------------------------------------------------------------

  ts <- terra::global(
    raster_stack,
    "mean",
    na.rm = TRUE
  )[, 1]

  if (length(ts) < 2) {
    stop("Time series too short to detect disturbance.")
  }

  if (all(is.na(ts))) {
    stop("Time series contains only NA values.")
  }

  # ------------------------------------------------------------------
  # Difference-based detection
  # ------------------------------------------------------------------

  if (detect == "difference") {

    diffs <- diff(ts)

    if (all(is.na(diffs))) {
      stop(
        "Unable to detect disturbance: all year-on-year ",
        "differences are NA."
      )
    }

    idx <- which.min(diffs)

    event_year <- years[idx + 1]

  } else {

    # --------------------------------------------------------------
    # LandTrendr-style local OLS slope detection
    # --------------------------------------------------------------

    window <- min(3L, length(ts) - 1L)

    n_windows <- length(ts) - window

    if (n_windows < 1L) {
      stop("Insufficient observations for local slope detection.")
    }

    slopes <- vapply(
      seq_len(n_windows),
      function(i) {

        y <- ts[i:(i + window)]

        ok <- is.finite(y)

        if (sum(ok) < 2L) {
          return(NA_real_)
        }

        x <- seq_along(y)

        fit <- stats::lm(y[ok] ~ x[ok])

        slope <- stats::coef(fit)[2]

        if (!is.finite(slope)) {
          return(NA_real_)
        }

        slope
      },
      numeric(1)
    )

    if (all(is.na(slopes))) {
      stop("Unable to detect disturbance: all local slopes are NA.")
    }

    idx <- which.min(slopes)

    event_year <- years[idx + 1]
  }

  # ------------------------------------------------------------------
  # Return
  # ------------------------------------------------------------------

  list(
    event_year = event_year,
    mask = mask_layer
  )
}
