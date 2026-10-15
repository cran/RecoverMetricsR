#' Compute Pixel-wise Spectral Recovery Metrics
#'
#' Computes pixel-wise spectral recovery metrics from a disturbance-year
#' raster time series. Metrics include disturbance magnitude (DM), OLS
#' recovery slope, early and late recovery rates, recovery-rate change,
#' Relative Recovery Indicator (RRI), Ratio of Eighty Percent (R80P),
#' and Year-on-Year Average (YrYr).
#'
#' The RRI, R80P, and YrYr metrics are based on methods described by
#' Frazier et al. (2018), with adaptations for a flexible number of
#' post-disturbance observations. The recovery-rate change metric is a
#' package-specific extension comparing early and late recovery rates.
#'
#' @param raster_stack A \code{terra::SpatRaster} containing the spectral
#'   index time series. Layers must correspond to the years supplied in
#'   \code{years}.
#' @param years Numeric vector giving the year corresponding to each raster
#'   layer. Must be sorted in strictly increasing order.
#' @param event_year Numeric value specifying the disturbance year.
#' @param pre_window Positive integer specifying the number of observations
#'   before the disturbance used to calculate the pre-disturbance baseline.
#' @param post_window Positive integer specifying the number of observations
#'   after the disturbance used for recovery analysis.
#' @param mask_layer Optional spatial mask. Can be a
#'   \code{terra::SpatVector}, \code{terra::SpatRaster}, or an \code{sf}
#'   object. Raster masks are resampled to \code{raster_stack} if their
#'   CRS, resolution, or extent differ.
#' @param metrics Character vector specifying the metrics to return.
#'   Available metrics are \code{"DM"}, \code{"Slope"},
#'   \code{"YrYr_early"}, \code{"YrYr_late"},
#'   \code{"Recovery_Rate_Change"}, \code{"RRI"}, \code{"R80P"},
#'   and \code{"YrYr"}.
#' @param clamp Logical. If \code{TRUE}, values outside the ranges specified
#'   in \code{clamp_ranges} are replaced by \code{NA}.
#' @param clamp_ranges Named list defining the lower and upper bounds for
#'   individual metrics. Names must correspond to metric names and each
#'   element must be a length-2 numeric vector with \code{lower < upper}.
#'   Metrics not listed are not clamped. Defaults are tuned for NBR-like
#'   indices.
#' @param eps Numeric. Small positive threshold below which denominators in
#'   ratio metrics are treated as zero, producing \code{NA}. Default
#'   \code{1e-6}.
#'
#' @details
#' The pre-disturbance baseline is calculated as the mean of the specified
#' pre-disturbance layers. The disturbance value corresponds to the layer
#' matching \code{event_year}.
#'
#' The disturbance magnitude is calculated as:
#'
#' \deqn{DM = Index_{post,final} - Index_{event}}
#'
#' The OLS recovery slope is estimated from all available post-disturbance
#' observations using the actual year offsets from the disturbance year.
#'
#' \code{YrYr} represents the average annual recovery from the disturbance
#' year to the final post-disturbance observation:
#'
#' \deqn{
#' YrYr = \frac{Index_{post,final} - Index_{event}}
#' {t_{final} - t_{event}}
#' }
#'
#' where the denominator represents the actual elapsed number of years.
#'
#' The early and late recovery rates are calculated from the mean spectral
#' values during the early and late portions of the post-disturbance
#' trajectory, respectively, relative to their mean year offsets.
#'
#' The recovery-rate change is defined as:
#'
#' \deqn{
#' Recovery\_Rate\_Change = YrYr_{early} - YrYr_{late}
#' }
#'
#' Positive values indicate that the recovery rate decreases from the early
#' to the late period, whereas negative values indicate an increase in the
#' recovery rate.
#'
#' The Relative Recovery Indicator is calculated as:
#'
#' \deqn{
#' RRI =
#' \frac{Index_{post,max} - Index_{event}}
#' {Index_{pre} - Index_{event}}
#' }
#'
#' The Ratio of Eighty Percent is calculated as:
#'
#' \deqn{
#' R80P =
#' \frac{Index_{post,max}}
#' {0.8 \times Index_{pre}}
#' }
#'
#' These formulations are based on methods described by Frazier et al.
#' (2018), with adaptations to the flexible post-disturbance window used
#' here.
#'
#' The \code{post_window} argument specifies the number of raster layers
#' used after the disturbance. If the supplied \code{years} contain gaps,
#' the actual elapsed time is used when calculating \code{YrYr} and the
#' OLS slope, but layer selection still follows the supplied layer order.
#'
#' The default bounds in \code{clamp_ranges} are tuned for NBR-like
#' indices that fall roughly within \code{[-1, 1]}. Users working with
#' other spectral indices should supply their own bounds. Clamping
#' replaces values outside the range with \code{NA}, which is useful for
#' removing unphysical fits but produces additional \code{NA} values that
#' would otherwise be finite. Set \code{clamp = FALSE} to inspect raw
#' values.
#'
#' The \code{RRI} and \code{R80P} metrics are ratio-based and become
#' undefined when their denominators are close to zero. Specifically,
#' \code{RRI} requires \code{Index_pre - Index_event != 0}, and
#' \code{R80P} requires \code{Index_pre != 0}. Pixels where the
#' denominator falls below \code{eps} are returned as \code{NA}. This
#' occurs most often in pixels where the pre-disturbance and
#' disturbance-year values are similar (e.g., unburned pixels inside a
#' fire perimeter), and results in \code{RRI} and \code{R80P} containing
#' more \code{NA} values than the difference-based metrics.
#'
#' When the input raster includes large areas of unburned vegetation, the
#' spatial distribution of the metrics will be dominated by pixels with
#' little spectral change. For fire-specific analyses, restrict the
#' analysis to the disturbance area using \code{mask_layer} or by
#' cropping \code{raster_stack} before calling this function. This
#' produces more interpretable distributions and reduces the fraction of
#' pixels affected by clamping.
#'
#' @references
#' Frazier, R. J., Coops, N. C., Wulder, M. A., Kennedy, R. E.,
#' & Hogg, E. H. (2018). Analyzing spatial and temporal variability
#' in short-term rates of post-fire vegetation return from Landsat
#' time series. \emph{Remote Sensing of Environment}, 205, 32--45.
#' \doi{10.1016/j.rse.2017.11.007}
#'
#' Storey, E. A., Stow, D. A., O'Leary, J. F., & Riggan, P. J. (2016).
#' Evaluating post-fire recovery of chaparral vegetation using Landsat
#' time series. \emph{Remote Sensing of Environment}, 183, 53--64.
#' \doi{10.1016/j.rse.2016.05.018}
#'
#' @return A \code{terra::SpatRaster} containing the requested metrics.
#'
#' @examples
#' library(terra)
#'
#' # Small toy example -- runs in well under 5 seconds
#' r <- rast(nrows = 2, ncols = 2, nlyrs = 6)
#' values(r) <- rep(c(0.60, 0.20, 0.30, 0.38, 0.45, 0.50), each = 4)
#' names(r) <- paste0("Index_", 2020:2025)
#'
#' # Compute all available metrics
#' result <- compute_recovery_metrics(
#'   raster_stack = r,
#'   years = 2020:2025,
#'   event_year = 2021,
#'   pre_window = 1,
#'   post_window = 4
#' )
#' result
#'
#' # Select a subset of metrics
#' selected <- compute_recovery_metrics(
#'   raster_stack = r,
#'   years = 2020:2025,
#'   event_year = 2021,
#'   pre_window = 1,
#'   post_window = 4,
#'   metrics = c("RRI", "R80P", "YrYr", "Recovery_Rate_Change")
#' )
#' selected
#'
#' \dontshow{
#' # Internal checks (not shown to users)
#' stopifnot(terra::nlyr(result) == 8L)
#' stopifnot(terra::nlyr(selected) == 4L)
#' stopifnot(
#'   all(names(result) ==
#'     c("DM", "Slope", "YrYr_early", "YrYr_late",
#'       "Recovery_Rate_Change", "RRI", "R80P", "YrYr"))
#' )
#' stopifnot(
#'   all(names(selected) ==
#'     c("RRI", "R80P", "YrYr", "Recovery_Rate_Change"))
#' )
#'
#' # Clamping removes out-of-range values
#' rri_vals <- terra::values(result[["RRI"]], mat = FALSE)
#' stopifnot(all(is.na(rri_vals) | (rri_vals >= -2 & rri_vals <= 3)))
#'
#' # Unclamped call still returns the expected number of layers
#' result_raw <- compute_recovery_metrics(
#'   raster_stack = r,
#'   years = 2020:2025,
#'   event_year = 2021,
#'   pre_window = 1,
#'   post_window = 4,
#'   clamp = FALSE
#' )
#' stopifnot(terra::nlyr(result_raw) == 8L)
#' }
#'
#' \donttest{
#' # Example with a spatial mask (SpatVector)
#' v <- vect(
#'   "POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))",
#'   crs = crs(r)
#' )
#'
#' masked <- compute_recovery_metrics(
#'   raster_stack = r,
#'   years = 2020:2025,
#'   event_year = 2021,
#'   pre_window = 1,
#'   post_window = 4,
#'   mask_layer = v,
#'   metrics = c("RRI", "YrYr")
#' )
#' masked
#' }
#'
#' @importFrom terra app ifel mask nlyr project rast resample same.crs vect
#' @export
compute_recovery_metrics <- function(
  raster_stack,
  years,
  event_year,
  pre_window = 2,
  post_window = 5,
  mask_layer = NULL,
  metrics = c(
    "DM",
    "Slope",
    "YrYr_early",
    "YrYr_late",
    "Recovery_Rate_Change",
    "RRI",
    "R80P",
    "YrYr"
  ),
  clamp = TRUE,
  clamp_ranges = list(
    DM = c(-1, 1),
    Slope = c(-0.2, 0.2),
    YrYr_early = c(-0.5, 0.5),
    YrYr_late = c(-0.5, 0.5),
    Recovery_Rate_Change = c(-1, 1),
    RRI = c(-2, 3),
    R80P = c(-1, 4),
    YrYr = c(-0.5, 0.5)
  ),
  eps = 1e-6
) {

  # ------------------------------------------------------------------
  # Input validation
  # ------------------------------------------------------------------

  if (!inherits(raster_stack, "SpatRaster")) {
    stop("raster_stack must be a terra::SpatRaster.")
  }

  if (!is.numeric(years)) {
    stop("years must be numeric.")
  }

  if (length(years) != terra::nlyr(raster_stack)) {
    stop("Length of years must equal the number of raster layers.")
  }

  if (anyDuplicated(years)) {
    stop("years must not contain duplicated values.")
  }

  if (is.unsorted(years, strictly = TRUE)) {
    stop("years must be sorted in strictly increasing order.")
  }

  if (!is.numeric(event_year) ||
      length(event_year) != 1 ||
      is.na(event_year)) {
    stop("event_year must be a single non-missing numeric value.")
  }

  if (!is.numeric(pre_window) ||
      length(pre_window) != 1 ||
      is.na(pre_window) ||
      pre_window < 1 ||
      pre_window != floor(pre_window)) {
    stop("pre_window must be a positive integer.")
  }

  if (!is.numeric(post_window) ||
      length(post_window) != 1 ||
      is.na(post_window) ||
      post_window < 1 ||
      post_window != floor(post_window)) {
    stop("post_window must be a positive integer.")
  }

  if (!is.numeric(eps) ||
      length(eps) != 1 ||
      is.na(eps) ||
      eps <= 0) {
    stop("eps must be a single positive numeric value.")
  }

  if (!is.logical(clamp) ||
      length(clamp) != 1 ||
      is.na(clamp)) {
    stop("clamp must be TRUE or FALSE.")
  }

  pre_window <- as.integer(pre_window)
  post_window <- as.integer(post_window)

  available_metrics <- c(
    "DM",
    "Slope",
    "YrYr_early",
    "YrYr_late",
    "Recovery_Rate_Change",
    "RRI",
    "R80P",
    "YrYr"
  )

  if (!is.character(metrics)) {
    stop("metrics must be a character vector.")
  }

  if (length(metrics) < 1L) {
    stop("metrics must contain at least one metric name.")
  }

  if (anyNA(metrics) || any(metrics == "")) {
    stop("metrics must contain valid non-empty names.")
  }

  if (anyDuplicated(metrics)) {
    stop("metrics must not contain duplicated names.")
  }

  invalid_metrics <- setdiff(metrics, available_metrics)

  if (length(invalid_metrics) > 0) {
    stop(
      "Unknown metric(s): ",
      paste(invalid_metrics, collapse = ", ")
    )
  }

  # ------------------------------------------------------------------
  # Validate clamp ranges
  # ------------------------------------------------------------------

  if (!is.list(clamp_ranges)) {
    stop("clamp_ranges must be a named list.")
  }

  if (length(clamp_ranges) > 0) {

    if (is.null(names(clamp_ranges)) ||
        any(names(clamp_ranges) == "")) {
      stop("clamp_ranges must be a fully named list.")
    }

    unknown_ranges <- setdiff(
      names(clamp_ranges),
      available_metrics
    )

    if (length(unknown_ranges) > 0) {
      stop(
        "clamp_ranges contains unknown metric(s): ",
        paste(unknown_ranges, collapse = ", ")
      )
    }

    for (nm in names(clamp_ranges)) {

      r <- clamp_ranges[[nm]]

      if (!is.numeric(r) ||
          length(r) != 2 ||
          anyNA(r) ||
          r[1] >= r[2]) {
        stop(
          "clamp_ranges$", nm,
          " must be a length-2 numeric vector with lower < upper."
        )
      }
    }
  }

  # ------------------------------------------------------------------
  # Optional mask
  # ------------------------------------------------------------------

  if (!is.null(mask_layer)) {

    if (inherits(mask_layer, "sf")) {
      mask_layer <- terra::vect(mask_layer)
    }

    if (!inherits(mask_layer, c("SpatVector", "SpatRaster"))) {
      stop(
        "mask_layer must be a terra::SpatVector, terra::SpatRaster, ",
        "or sf object."
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

      if (!terra::same.crs(raster_stack, mask_layer)) {
        mask_layer <- terra::project(
          mask_layer,
          raster_stack
        )
      }
    }

    raster_stack <- terra::mask(
      raster_stack,
      mask_layer
    )
  }

  # ------------------------------------------------------------------
  # Helper functions
  # ------------------------------------------------------------------

  # Raster-safe division. The denominator is masked before division so
  # that Inf / NaN never appear in intermediate results.
  safe_div <- function(a, b, eps = 1e-6) {

    if (!inherits(a, "SpatRaster")) {
      stop("'a' must be a terra::SpatRaster.")
    }

    invalid <- is.na(b) | abs(b) < eps

    if (!inherits(invalid, "SpatRaster")) {
      # Scalar denominator: no need for terra::ifel
      if (isTRUE(invalid)) {
        return(a * NA_real_)
      }
      return(a / b)
    }

    b_safe <- terra::ifel(invalid, NA_real_, b)

    a / b_safe
  }

  # NA-safe mean.
  safe_mean <- function(x) {

    ok <- is.finite(x)

    if (!any(ok)) {
      return(NA_real_)
    }

    mean(x[ok])
  }

  # NA-safe maximum.
  safe_max <- function(x) {

    ok <- is.finite(x)

    if (!any(ok)) {
      return(NA_real_)
    }

    max(x[ok])
  }

  # Replace values outside [lower, upper] with NA.
  clamp_metric <- function(x, range) {

    if (is.null(range)) {
      return(x)
    }

    terra::ifel(
      x < range[1] | x > range[2],
      NA_real_,
      x
    )
  }

  # ------------------------------------------------------------------
  # Disturbance year
  # ------------------------------------------------------------------

  idx_event <- which(years == event_year)

  if (length(idx_event) != 1) {
    stop(
      "event_year must correspond to exactly one value in years."
    )
  }

  # ------------------------------------------------------------------
  # Pre-disturbance baseline
  # ------------------------------------------------------------------

  pre_ids <- (idx_event - pre_window):(idx_event - 1)

  if (min(pre_ids) < 1) {
    stop(
      "Not enough pre-disturbance observations for pre_window."
    )
  }

  pre_stack <- raster_stack[[pre_ids]]

  Index_pre <- if (terra::nlyr(pre_stack) == 1L) {
    pre_stack[[1]]
  } else {
    terra::app(pre_stack, fun = safe_mean)
  }

  Index_event <- raster_stack[[idx_event]]

  # ------------------------------------------------------------------
  # Post-disturbance observations
  # ------------------------------------------------------------------

  last_post_id <- min(
    idx_event + post_window,
    terra::nlyr(raster_stack)
  )

  post_ids <- (idx_event + 1):last_post_id

  if (length(post_ids) < 2) {
    stop(
      "At least two post-disturbance observations are required."
    )
  }

  post_stack <- raster_stack[[post_ids]]

  n_post <- terra::nlyr(post_stack)

  years_post <- years[post_ids] - event_year

  # ------------------------------------------------------------------
  # Early and late periods
  # ------------------------------------------------------------------

  split_id <- floor(n_post / 2)

  early_ids <- seq_len(split_id)
  late_ids  <- (split_id + 1):n_post

  early_offset <- mean(years_post[early_ids])
  late_offset  <- mean(years_post[late_ids])

  # ------------------------------------------------------------------
  # Post-disturbance summaries
  # ------------------------------------------------------------------

  Index_post_max <- terra::app(
    post_stack,
    fun = safe_max
  )

  Index_post_final <- post_stack[[n_post]]

  Index_post_early <- if (length(early_ids) == 1L) {
    post_stack[[early_ids]]
  } else {
    terra::app(
      post_stack[[early_ids]],
      fun = safe_mean
    )
  }

  Index_post_late <- if (length(late_ids) == 1L) {
    post_stack[[late_ids]]
  } else {
    terra::app(
      post_stack[[late_ids]],
      fun = safe_mean
    )
  }

  # ------------------------------------------------------------------
  # Disturbance magnitude
  # ------------------------------------------------------------------

  DM <- Index_post_final - Index_event

  # ------------------------------------------------------------------
  # Relative Recovery Indicator (RRI)
  # ------------------------------------------------------------------

  RRI <- safe_div(
    Index_post_max - Index_event,
    Index_pre - Index_event,
    eps = eps
  )

  # ------------------------------------------------------------------
  # Ratio of Eighty Percent (R80P)
  # ------------------------------------------------------------------

  R80P <- safe_div(
    Index_post_max,
    0.8 * Index_pre,
    eps = eps
  )

  # ------------------------------------------------------------------
  # Year-on-Year Average (YrYr)
  #
  # Uses actual elapsed years rather than the number of layers.
  # ------------------------------------------------------------------

  elapsed_years <- years_post[n_post]

  YrYr <- safe_div(
    Index_post_final - Index_event,
    elapsed_years,
    eps = eps
  )

  # ------------------------------------------------------------------
  # Early recovery rate
  # ------------------------------------------------------------------

  YrYr_early <- safe_div(
    Index_post_early - Index_event,
    early_offset,
    eps = eps
  )

  # ------------------------------------------------------------------
  # Late recovery rate
  # ------------------------------------------------------------------

  YrYr_late <- safe_div(
    Index_post_late - Index_event,
    late_offset,
    eps = eps
  )

  # ------------------------------------------------------------------
  # Recovery-rate change
  #
  # Positive = decelerating recovery.
  # Negative = accelerating recovery.
  #
  # This is a package-specific metric and should not be confused
  # with the Difference Recovery Metric (DRM) of Storey et al.
  # ------------------------------------------------------------------

  Recovery_Rate_Change <- YrYr_early - YrYr_late

  # ------------------------------------------------------------------
  # OLS recovery slope
  #
  # Analytic, NA-safe calculation using actual year offsets.
  # ------------------------------------------------------------------

  Slope <- terra::app(
    post_stack,
    fun = function(x) {

      ok <- is.finite(x)

      if (sum(ok) < 2) {
        return(NA_real_)
      }

      xv <- x[ok]
      yv <- years_post[ok]

      xm <- mean(xv)
      ym <- mean(yv)

      denominator <- sum((yv - ym)^2)

      if (!is.finite(denominator) || denominator == 0) {
        return(NA_real_)
      }

      numerator <- sum((yv - ym) * (xv - xm))

      numerator / denominator
    }
  )

  # ------------------------------------------------------------------
  # Store metrics
  # ------------------------------------------------------------------

  all_metrics <- list(
    DM                   = DM,
    Slope                = Slope,
    YrYr_early           = YrYr_early,
    YrYr_late            = YrYr_late,
    Recovery_Rate_Change = Recovery_Rate_Change,
    RRI                  = RRI,
    R80P                 = R80P,
    YrYr                 = YrYr
  )

  # ------------------------------------------------------------------
  # Optional clamping
  # ------------------------------------------------------------------

  if (isTRUE(clamp)) {

    for (metric_name in names(all_metrics)) {

      if (!is.null(clamp_ranges[[metric_name]])) {

        all_metrics[[metric_name]] <- clamp_metric(
          all_metrics[[metric_name]],
          clamp_ranges[[metric_name]]
        )
      }
    }
  }

  # ------------------------------------------------------------------
  # Return selected metrics
  # ------------------------------------------------------------------

  output <- terra::rast(all_metrics[metrics])

  names(output) <- metrics

  output
}
