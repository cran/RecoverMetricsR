#' @title Fit Spectral Recovery Model with Optional Model Fit
#'
#' @description
#' Plots mean spectral recovery trajectories following disturbance.
#' Optionally fits parametric recovery models and overlays fitted curves.
#'
#' @param raster_stack A \code{terra::SpatRaster} of yearly spectral index
#'   layers.
#' @param years Numeric vector of years. Its length must equal the number
#'   of raster layers.
#' @param event_year Numeric disturbance year.
#' @param mask_layer Optional spatial mask. Can be a
#'   \code{terra::SpatRaster}, \code{terra::SpatVector}, or an \code{sf}
#'   object defining regions. If an \code{sf} or \code{SpatVector} object
#'   with multiple features is supplied, a separate trajectory is computed
#'   for each feature.
#' @param index_name Character. Y-axis label. Default
#'   \code{"Spectral index"}.
#' @param ci Either \code{"sd"} to show plus/minus one standard deviation
#'   around the mean, or \code{FALSE} to suppress the error bars.
#' @param model \code{NULL} to skip model fitting, a single model name
#'   (\code{"exponential"}, \code{"logistic"}, or \code{"gompertz"}), or
#'   \code{"compare"} to fit all three and rank them by AIC.
#'
#' @details
#' Model fitting pools all mask regions into a single mean trajectory.
#' To fit separate models per region, call this function once per
#' region.
#'
#' When \code{model} is a single model name, the fitted curve is overlaid
#' on the plot and an equation label is added. When
#' \code{model = "compare"}, all three models are fitted, but no curve
#' is overlaid; the metrics table includes AIC and BIC for ranking.
#' Models that fail to converge are skipped with a warning.
#'
#' The fitted trajectory consists of the layers strictly after the
#' disturbance year, consistent with the other functions in the package.
#'
#' When \code{mask_layer = NULL}, the standard deviation shown by the
#' error bars reflects the spatial spread of all pixels in the raster,
#' including unburned areas. To plot the spread within a specific region,
#' supply \code{mask_layer}.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{plot}{A \code{ggplot2} object.}
#'   \item{metrics}{A data frame of fit statistics, or an empty data
#'     frame when \code{model = NULL}.}
#' }
#'
#' @examples
#' library(terra)
#'
#' r <- rast(nrows = 4, ncols = 4, nlyrs = 6)
#' values(r) <- rep(c(0.60, 0.20, 0.30, 0.38, 0.45, 0.50), each = 16)
#' names(r) <- paste0("Index_", 2020:2025)
#'
#' # Trajectory only, no model fit
#' out <- fit_recovery_model(
#'   raster_stack = r,
#'   years = 2020:2025,
#'   event_year = 2021
#' )
#'
#' out$plot
#'
#' \dontshow{
#' # Internal checks (not shown to users)
#' stopifnot(is.list(out))
#' stopifnot(inherits(out$plot, "ggplot"))
#' stopifnot(is.data.frame(out$metrics))
#' }
#'
#' \donttest{
#' # Fit a single exponential recovery model
#' fit_exp <- fit_recovery_model(
#'   raster_stack = r,
#'   years = 2020:2025,
#'   event_year = 2021,
#'   model = "exponential"
#' )
#'
#' fit_exp$plot
#' fit_exp$metrics
#'
#' # Compare exponential, logistic, and Gompertz
#' fit_cmp <- fit_recovery_model(
#'   raster_stack = r,
#'   years = 2020:2025,
#'   event_year = 2021,
#'   model = "compare"
#' )
#'
#' fit_cmp$metrics
#' }
#'
#' @importFrom stats sd quantile AIC BIC predict coef
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_vline
#'   geom_errorbar geom_label labs theme_classic
#' @importFrom grid unit
#' @importFrom terra nlyr global mask vect rasterize
#'   same.crs resample res ext project
#' @export
fit_recovery_model <- function(
  raster_stack,
  years,
  event_year,
  mask_layer = NULL,
  index_name = "Spectral index",
  ci = "sd",
  model = NULL
) {

  # ---- Input validation ----

  if (!inherits(raster_stack, "SpatRaster"))
    stop("raster_stack must be a terra::SpatRaster")

  if (length(years) != terra::nlyr(raster_stack))
    stop("Length of 'years' must match raster layers")

  if (!identical(ci, "sd") && !identical(ci, FALSE))
    stop("'ci' must be either \"sd\" or FALSE.")

  # ---- Build mask_layer list ----
  #
  # For vector masks (sf or SpatVector) with N features, one mask raster
  # is built per feature. Raster masks are used directly.

  mask_layer_list <- list(All = NULL)

  if (!is.null(mask_layer)) {

    # -- Raster mask -------------------------------------------------
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

      mask_layer_list <- list(Mask = mask_layer)

      # -- Vector mask (sf or SpatVector) ------------------------------
    } else {

      # Convert sf to SpatVector if needed
      if (inherits(mask_layer, "sf"))
        mask_layer <- terra::vect(mask_layer)

      if (!inherits(mask_layer, "SpatVector"))
        stop(
          "'mask_layer' must be a terra::SpatRaster, terra::SpatVector, ",
          "or sf object"
        )

      if (nrow(mask_layer) == 0)
        stop("'mask_layer' has zero features.")

      # Reproject to the raster CRS if needed
      if (!terra::same.crs(raster_stack, mask_layer)) {
        mask_layer <- terra::project(mask_layer, raster_stack)
      }

      mask_layer_list <- lapply(
        seq_len(nrow(mask_layer)),
        function(i) {
          terra::rasterize(mask_layer[i, ], raster_stack[[1]])
        }
      )

      names(mask_layer_list) <- paste0(
        "Region_",
        seq_along(mask_layer_list)
      )
    }
  }

  # ---- Build tidy trajectory data ----

  df <- do.call(
    rbind,
    lapply(names(mask_layer_list), function(nm) {

      r_use <- raster_stack

      if (!is.null(mask_layer_list[[nm]]))
        r_use <- terra::mask(r_use, mask_layer_list[[nm]])

      data.frame(
        Year = as.numeric(years),
        Mean = terra::global(r_use, "mean", na.rm = TRUE)[, 1],
        SD   = terra::global(r_use, "sd",   na.rm = TRUE)[, 1],
        Mask = nm,
        stringsAsFactors = FALSE
      )
    })
  )

  # ---- Base plot ----

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(.data$Year, .data$Mean, color = .data$Mask)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_vline(xintercept = event_year, linetype = "dashed") +
    ggplot2::labs(y = index_name, x = "Year") +
    ggplot2::theme_classic()

  if (identical(ci, "sd")) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$Mean - .data$SD,
                   ymax = .data$Mean + .data$SD),
      width = 0.2, alpha = 0.6
    )
  }

  if (is.null(model))
    return(list(plot = p, metrics = data.frame()))

  # ---- Model fitting ----
  #
  # Strictly post-disturbance years, consistent with the rest of the
  # package (compute_recovery_metrics, compute_recovery_time,
  # predict_recovery_raster).

  df_post <- df[df$Year > event_year, , drop = FALSE]
  x <- df_post$Year - event_year
  y <- df_post$Mean

  if (sum(is.finite(y)) < 4)
    stop("Not enough finite post-disturbance observations to fit a model")

  models <- if (identical(model, "compare"))
    c("exponential", "logistic", "gompertz") else model

  # ---- Fit all requested models, skipping failures gracefully ----

  metrics <- do.call(
    rbind,
    lapply(models, function(m) {

      fit <- fit_recovery_model_core(x, y, m)

      if (is.null(fit)) {
        warning("Model '", m, "' failed to converge; skipping.",
                call. = FALSE)
        return(NULL)
      }

      pred <- tryCatch(
        stats::predict(fit$fit, newdata = data.frame(t = x)),
        error = function(e) NULL
      )

      if (is.null(pred)) {
        warning("Model '", m, "' produced no predictions; skipping.",
                call. = FALSE)
        return(NULL)
      }

      resid <- y - pred

      if (all(!is.finite(resid))) {
        ME   <- NA_real_
        RMSE <- NA_real_
        R2   <- NA_real_
      } else {
        ME   <- mean(resid, na.rm = TRUE)
        RMSE <- sqrt(mean(resid^2, na.rm = TRUE))
        R2   <- 1 - sum(resid^2, na.rm = TRUE) /
          sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
      }

      data.frame(
        Model = m,
        ME    = ME,
        RMSE  = RMSE,
        R2    = R2,
        AIC   = tryCatch(stats::AIC(fit$fit), error = function(e) NA_real_),
        BIC   = tryCatch(stats::BIC(fit$fit), error = function(e) NA_real_),
        stringsAsFactors = FALSE
      )
    })
  )

  if (is.null(metrics) || nrow(metrics) == 0)
    stop("No model converged on the provided trajectory.")

  metrics$DeltaAIC <- metrics$AIC - min(metrics$AIC, na.rm = TRUE)

  # ---- Overlay single-model fit (robust to convergence failure) ----

  if (length(models) == 1) {

    fit <- fit_recovery_model_core(x, y, models)

    if (is.null(fit)) {

      warning("Model '", models,
              "' failed to converge; returning plot without fitted curve.",
              call. = FALSE)

    } else {

      fitted <- tryCatch(
        stats::predict(fit$fit, newdata = data.frame(t = x)),
        error = function(e) NULL
      )

      if (!is.null(fitted)) {

        df_post$Fitted <- fitted

        y_lab <- stats::quantile(df_post$Mean, 0.1, na.rm = TRUE)

        if (!is.finite(y_lab))
          y_lab <- min(df_post$Mean, na.rm = TRUE)

        p <- p +
          ggplot2::geom_line(
            data = df_post,
            ggplot2::aes(.data$Year, .data$Fitted),
            color = "black", linewidth = 1.5
          ) +
          ggplot2::geom_label(
            data = data.frame(
              x = mean(df_post$Year, na.rm = TRUE),
              y = y_lab,
              label = .format_equation(models, stats::coef(fit$fit)),
              stringsAsFactors = FALSE
            ),
            ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
            hjust = 0.5, vjust = 1, size = 5,
            label.size = 0.5,
            label.padding = grid::unit(0.25, "lines"),
            fill = "white", color = "black"
          )
      }
    }
  }

  list(plot = p, metrics = metrics)
}
