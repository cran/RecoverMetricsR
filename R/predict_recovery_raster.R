#' @title Predict Spectral Recovery Raster
#'
#' @description
#' Predicts post-disturbance spectral recovery pixel by pixel using a
#' fitted recovery model. For each pixel, a parametric model is fitted to
#' the post-disturbance trajectory, then extrapolated over a user-defined
#' horizon.
#'
#' @param raster_stack A \code{terra::SpatRaster} of yearly spectral
#'   index layers.
#' @param years Numeric vector of calendar years corresponding to raster
#'   layers. Its length must equal the number of raster layers.
#' @param event_year Numeric disturbance year.
#' @param model Character. One of \code{"exponential"},
#'   \code{"logistic"}, or \code{"gompertz"}.
#' @param horizon Integer. Number of years to predict after the last
#'   observation.
#' @param mask_layer Optional spatial mask. Can be a
#'   \code{terra::SpatRaster}, \code{terra::SpatVector}, or an \code{sf}
#'   object defining pixels to include. Raster masks are resampled to
#'   \code{raster_stack} if their CRS, resolution, or extent differ.
#'   Vector masks are reprojected to the raster CRS if needed.
#' @param max_abs Numeric. Maximum allowed absolute predicted value.
#'   Predictions exceeding this threshold are treated as failed fits and
#'   returned as \code{NA}. Default is 2.
#'
#' @details
#' The fitted trajectory consists of the layers strictly after the
#' disturbance year. The time axis is defined in years since the
#' disturbance. A pixel is skipped (predicted as \code{NA}) if fewer than
#' four finite observations are available, if the model fails to converge,
#' or if the predicted values fall outside the plausible range specified
#' by \code{max_abs}.
#'
#' Multi-layer output from \code{terra::app()} requires
#' \pkg{terra} >= 1.5-0. The package DESCRIPTION declares
#' \code{terra >= 1.8}, so this requirement is always satisfied.
#'
#' @return A \code{terra::SpatRaster} with one layer per prediction year,
#'   named \code{Year_<calendar year>}.
#'
#' @examples
#' library(terra)
#'
#' r <- rast(nrows = 4, ncols = 4, nlyrs = 6)
#' values(r) <- rep(c(0.60, 0.20, 0.30, 0.38, 0.45, 0.50), each = 16)
#' names(r) <- paste0("Index_", 2020:2025)
#'
#' pred <- predict_recovery_raster(
#'   raster_stack = r,
#'   years = 2020:2025,
#'   event_year = 2021,
#'   model = "exponential",
#'   horizon = 3
#' )
#'
#' pred
#'
#' \dontshow{
#' # Internal checks (not shown to users)
#' n_post <- sum(2020:2025 > 2021)
#' stopifnot(inherits(pred, "SpatRaster"))
#' stopifnot(terra::nlyr(pred) == n_post + 3L)
#' stopifnot(all(grepl("^Year_", names(pred))))
#' }
#'
#' \donttest{
#' # Logistic model, prediction horizon of five years
#' pred_log <- predict_recovery_raster(
#'   raster_stack = r,
#'   years = 2020:2025,
#'   event_year = 2021,
#'   model = "logistic",
#'   horizon = 5
#' )
#'
#' pred_log
#' }
#'
#' @importFrom terra nlyr mask vect app rast same.crs resample res ext project
#' @export
predict_recovery_raster <- function(
  raster_stack,
  years,
  event_year,
  model = c("exponential", "logistic", "gompertz"),
  horizon = 5,
  mask_layer = NULL,
  max_abs = 2
) {

  # ---- Input validation ----

  model <- match.arg(model)

  if (!inherits(raster_stack, "SpatRaster"))
    stop("raster_stack must be a terra::SpatRaster")

  if (!is.numeric(years))
    stop("'years' must be numeric (calendar years)")

  if (length(years) != terra::nlyr(raster_stack))
    stop("Length of 'years' must match raster layers")

  if (!is.numeric(horizon) ||
      length(horizon) != 1 ||
      is.na(horizon) ||
      horizon < 1 ||
      horizon != floor(horizon)) {
    stop("'horizon' must be a positive integer.")
  }

  horizon <- as.integer(horizon)

  if (!is.numeric(max_abs) ||
      length(max_abs) != 1 ||
      is.na(max_abs) ||
      max_abs <= 0) {
    stop("'max_abs' must be a single positive numeric value.")
  }

  # ---- Optional mask ----

  if (!is.null(mask_layer)) {

    mask_r <- if (inherits(mask_layer, "SpatRaster")) {
      mask_layer
    } else if (inherits(mask_layer, "SpatVector")) {
      mask_layer
    } else if (inherits(mask_layer, "sf")) {
      terra::vect(mask_layer)
    } else {
      stop(
        "'mask_layer' must be a terra::SpatRaster, terra::SpatVector, ",
        "or sf object"
      )
    }

    if (inherits(mask_r, "SpatRaster")) {

      needs_resample <- !terra::same.crs(raster_stack, mask_r) ||
        !isTRUE(all.equal(
          terra::res(raster_stack),
          terra::res(mask_r)
        )) ||
        !isTRUE(all.equal(
          as.vector(terra::ext(raster_stack)),
          as.vector(terra::ext(mask_r))
        ))

      if (needs_resample) {
        mask_r <- terra::resample(mask_r, raster_stack, method = "near")
      }

    } else {

      if (!terra::same.crs(raster_stack, mask_r)) {
        mask_r <- terra::project(mask_r, raster_stack)
      }
    }

    raster_stack <- terra::mask(raster_stack, mask_r)
  }

  # ---- Post-disturbance window ----

  idx_event <- which(years == event_year)

  if (length(idx_event) == 0)
    stop("event_year not found in years")

  if (length(idx_event) > 1)
    stop("event_year occurs more than once in years")

  # Exclude the disturbance year so the trajectory is strictly
  # post-disturbance, consistent with compute_recovery_metrics().
  post_ids <- (idx_event + 1):terra::nlyr(raster_stack)

  if (length(post_ids) < 4)
    stop("At least four post-disturbance observations are required.")

  r_post <- raster_stack[[post_ids]]

  t_obs  <- years[post_ids] - event_year
  t_ok   <- is.finite(t_obs)

  t_pred <- max(t_obs) + seq_len(horizon)
  t_all  <- c(t_obs, t_pred)

  # ---- Per-pixel prediction function ----

  pred_fun <- function(x) {

    ok <- is.finite(x) & t_ok

    if (sum(ok) < 4)
      return(rep(NA_real_, length(t_all)))

    y <- x[ok]
    t <- t_obs[ok]

    fit <- fit_recovery_model_core(t, y, model = model)

    if (is.null(fit))
      return(rep(NA_real_, length(t_all)))

    pars <- fit$coefficients

    y_hat <- switch(
      model,
      exponential = pars["A"] * (1 - exp(-pars["k"] * t_all)),
      logistic    = pars["A"] /
        (1 + exp(-pars["k"] * (t_all - pars["x0"]))),
      gompertz    = pars["A"] *
        exp(-exp(-pars["k"] * (t_all - pars["x0"])))
    )

    y_hat <- as.numeric(y_hat)

    if (any(!is.finite(y_hat)) || any(abs(y_hat) > max_abs))
      return(rep(NA_real_, length(t_all)))

    y_hat
  }

  # ---- Apply per pixel ----
  #
  # terra::app() produces one output layer per returned element when the
  # function returns a fixed-length vector. Requires terra >= 1.5-0,
  # guaranteed by the package DESCRIPTION.
  pred_stack <- terra::app(r_post, pred_fun)

  names(pred_stack) <- paste0("Year_", event_year + t_all)

  pred_stack
}
