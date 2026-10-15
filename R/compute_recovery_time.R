#' @title Compute Pixel-wise Recovery Time
#'
#' @description
#' Computes the number of years required for each pixel to reach a given
#' proportion of its pre-disturbance spectral value and remain above it
#' for a specified number of consecutive years.
#'
#' @param raster_stack A \code{terra::SpatRaster} of yearly spectral
#'   index layers.
#' @param years Numeric vector of years corresponding to raster layers.
#'   Its length must equal the number of raster layers.
#' @param event_year Numeric. Disturbance year.
#' @param proportion Numeric. Proportion of the pre-disturbance value
#'   used as the recovery threshold. Default is \code{0.8}.
#' @param pre_window Integer. Number of pre-disturbance observations used
#'   to compute the baseline.
#' @param persistence Integer. Number of consecutive years the pixel must
#'   remain above the threshold to be considered recovered. Default is 1.
#' @param mask_layer Optional spatial mask. Can be a
#'   \code{terra::SpatRaster}, \code{terra::SpatVector}, or \code{sf}
#'   object.
#' @param return_year Logical. If \code{TRUE}, return the calendar year
#'   of recovery instead of the elapsed time. Default is \code{FALSE}.
#'
#' @details
#' The pre-disturbance baseline is computed as the mean of the
#' \code{pre_window} observations before \code{event_year}. The recovery
#' threshold is \code{proportion * baseline}. For each pixel, the function
#' finds the first run of at least \code{persistence} consecutive
#' post-disturbance observations at or above the threshold.
#'
#' Missing values are treated as "not recovered" for the purposes of run
#' detection. The reported recovery time is the number of years between
#' the disturbance and the first year of the qualifying run.
#'
#' \code{pre_window} and \code{persistence} count raster layers
#' (observations), not calendar years.
#'
#' When the input raster includes unburned areas, those pixels will
#' typically be classified as recovered at the first post-disturbance
#' year, because their spectral values never dropped below the threshold.
#' This produces a spike at TTR = 1 in the unmasked distribution. To
#' obtain meaningful time-to-recovery distributions, restrict the
#' analysis to the disturbance area using \code{mask_layer} or by
#' cropping \code{raster_stack} before calling this function.
#'
#' @return A \code{terra::SpatRaster} with a single layer. If
#'   \code{return_year = FALSE}, the layer holds elapsed years since the
#'   disturbance. If \code{return_year = TRUE}, the layer holds the
#'   calendar year of recovery. Cells that never meet the recovery
#'   criterion are \code{NA}.
#'
#' @examples
#' library(terra)
#'
#' r <- rast(nrows = 4, ncols = 4, nlyrs = 6)
#' values(r) <- rep(c(0.60, 0.20, 0.30, 0.38, 0.45, 0.50), each = 16)
#' names(r) <- paste0("Index_", 2020:2025)
#'
#' # Elapsed years since disturbance
#' ttr <- compute_recovery_time(
#'   raster_stack = r,
#'   years = 2020:2025,
#'   event_year = 2021,
#'   proportion = 0.8,
#'   pre_window = 1,
#'   persistence = 1
#' )
#'
#' ttr
#'
#' \dontshow{
#' # Internal checks (not shown to users)
#' stopifnot(inherits(ttr, "SpatRaster"))
#' stopifnot(terra::nlyr(ttr) == 1L)
#'
#' # Persistence > 1 must enforce the run-length constraint
#' ttr2 <- compute_recovery_time(
#'   raster_stack = r,
#'   years = 2020:2025,
#'   event_year = 2021,
#'   pre_window = 1,
#'   persistence = 2
#' )
#' vals2 <- terra::values(ttr2, mat = FALSE)
#' stopifnot(all(is.na(vals2) | vals2 >= 2))
#' }
#'
#' \donttest{
#' # Return the calendar year of recovery instead
#' rec_year <- compute_recovery_time(
#'   raster_stack = r,
#'   years = 2020:2025,
#'   event_year = 2021,
#'   proportion = 0.8,
#'   pre_window = 1,
#'   persistence = 2,
#'   return_year = TRUE
#' )
#'
#' rec_year
#' }
#'
#' @importFrom terra nlyr app mask vect same.crs resample res ext project
#' @export
compute_recovery_time <- function(
  raster_stack,
  years,
  event_year,
  proportion = 0.8,
  pre_window = 2,
  persistence = 1,
  mask_layer = NULL,
  return_year = FALSE
) {

  # ---- Input validation ----

  if (!inherits(raster_stack, "SpatRaster"))
    stop("raster_stack must be a terra::SpatRaster")

  if (!is.numeric(years))
    stop("'years' must be numeric")

  if (length(years) != terra::nlyr(raster_stack))
    stop("Length of 'years' must match raster layers")

  if (!is.numeric(proportion) ||
      length(proportion) != 1 ||
      is.na(proportion) ||
      proportion <= 0 ||
      proportion > 1) {
    stop("proportion must be a single value in (0, 1]")
  }

  if (!is.numeric(pre_window) ||
      length(pre_window) != 1 ||
      is.na(pre_window) ||
      pre_window < 1 ||
      pre_window != floor(pre_window)) {
    stop("pre_window must be a positive integer")
  }

  if (!is.numeric(persistence) ||
      length(persistence) != 1 ||
      is.na(persistence) ||
      persistence < 1 ||
      persistence != floor(persistence)) {
    stop("persistence must be a positive integer")
  }

  if (!is.logical(return_year) ||
      length(return_year) != 1 ||
      is.na(return_year)) {
    stop("return_year must be TRUE or FALSE")
  }

  pre_window  <- as.integer(pre_window)
  persistence <- as.integer(persistence)

  # ---- Optional mask ----

  if (!is.null(mask_layer)) {

    if (inherits(mask_layer, "sf"))
      mask_layer <- terra::vect(mask_layer)

    if (!inherits(mask_layer, c("SpatVector", "SpatRaster")))
      stop("mask_layer must be sf, terra::SpatVector, terra::SpatRaster, or NULL")

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

      if (needs_resample)
        mask_layer <- terra::resample(mask_layer, raster_stack, method = "near")

    } else {
      if (!terra::same.crs(raster_stack, mask_layer))
        mask_layer <- terra::project(mask_layer, raster_stack)
    }

    raster_stack <- terra::mask(raster_stack, mask_layer)
  }

  # ---- Disturbance year ----

  idx_event <- which(years == event_year)

  if (length(idx_event) == 0)
    stop("event_year not found in years")

  if (length(idx_event) > 1)
    stop("event_year occurs more than once in years")

  # ---- Pre-disturbance baseline ----

  pre_ids <- (idx_event - pre_window):(idx_event - 1)

  if (min(pre_ids) < 1)
    stop("Not enough pre-disturbance observations for pre_window")

  pre_stack <- raster_stack[[pre_ids]]

  pre_mean <- if (terra::nlyr(pre_stack) == 1L) {
    pre_stack[[1]]
  } else {
    terra::app(pre_stack, mean, na.rm = TRUE)
  }

  threshold <- pre_mean * proportion

  # ---- Post-disturbance observations ----

  post_ids <- (idx_event + 1):terra::nlyr(raster_stack)

  if (length(post_ids) < persistence)
    stop("Not enough post-disturbance observations for persistence")

  post_stack <- raster_stack[[post_ids]]
  post_years <- years[post_ids]

  full_stack <- c(threshold, post_stack)

  # ---- Per-pixel recovery time ----

  rec_time <- terra::app(full_stack, fun = function(v) {

    thr <- v[1]
    x   <- v[-1]

    if (is.na(thr) || all(is.na(x)))
      return(NA_real_)

    above <- x >= thr
    above[is.na(above)] <- FALSE

    if (!any(above))
      return(NA_real_)

    r   <- rle(above)
    idx <- which(r$values & r$lengths >= persistence)[1]

    if (is.na(idx))
      return(NA_real_)

    year_idx <- sum(r$lengths[seq_len(idx - 1)]) + 1L

    as.numeric(post_years[year_idx] - event_year)
  })

  if (return_year) {

    vals <- rec_time[]
    rec_year <- terra::rast(rec_time)
    rec_year[] <- as.integer(vals + event_year)
    names(rec_year) <- "Recovery_Year"
    return(rec_year)
  }

  names(rec_time) <- paste0("TTR_", proportion * 100, "_P", persistence)
  rec_time
}
