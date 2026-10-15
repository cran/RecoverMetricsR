#' @title Classify Recovery Typology
#'
#' @description
#' Classifies post-disturbance recovery into four typological classes
#' according to mean time-to-recovery (TTR) and sensitivity of TTR to the
#' recovery threshold.
#'
#' @param ttr_stack A \code{terra::SpatRaster} produced by
#'   \code{\link{compute_recovery_time_uncertainty}}. It must contain at
#'   least two threshold-specific TTR layers and a
#'   \code{TTR_Threshold_Range} layer.
#' @param fast_TTR Numeric. Maximum TTR, in years, defining fast recovery.
#' @param unc_threshold Numeric. Maximum threshold range, in years, defining
#'   low sensitivity of the recovery-time estimate to the selected threshold.
#' @param plot_map Logical. If \code{TRUE}, plots the recovery typology map.
#' @param plot_stats Logical. If \code{TRUE}, plots the recovery typology
#'   statistics.
#' @param area Logical. If \code{TRUE}, computes the area of each recovery
#'   type in hectares.
#'
#' @details
#' Mean recovery time is calculated across all threshold-specific TTR
#' layers. The threshold range is obtained from the
#' \code{TTR_Threshold_Range} layer and represents the sensitivity of the
#' recovery-time estimate to the choice of recovery threshold.
#'
#' The typology combines two dimensions:
#' \itemize{
#'   \item recovery speed, based on whether mean TTR is less than or equal
#'     to \code{fast_TTR};
#'   \item threshold sensitivity, based on whether the TTR range is less
#'     than or equal to \code{unc_threshold}.
#' }
#'
#' Four classes are therefore produced:
#' \itemize{
#'   \item 1: Fast and stable recovery;
#'   \item 2: Fast but threshold-sensitive recovery;
#'   \item 3: Slow and stable recovery;
#'   \item 4: Slow and threshold-sensitive recovery.
#' }
#'
#' The terms "stable" and "threshold-sensitive" refer specifically to the
#' sensitivity of TTR to the selected recovery threshold. They do not imply
#' ecological stability or statistical uncertainty.
#'
#' The area calculation assumes the raster is in a projected coordinate
#' system with meters as the linear unit. For rasters in geographic
#' coordinates (degrees), project to a suitable CRS with
#' \code{terra::project()} first.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{TTR_Mean}{Raster containing mean TTR across recovery thresholds.}
#'   \item{TTR_Combined}{Raster containing mean TTR and threshold-range
#'     layers.}
#'   \item{Recovery_Type}{Integer raster containing the four recovery
#'     typology classes.}
#'   \item{Stats}{Data frame containing pixel counts, percentages, and,
#'     optionally, area in hectares for each recovery type.}
#'   \item{Stats_Plot}{A \code{ggplot2} object showing recovery typology
#'     percentages, or \code{NULL} if \code{plot_stats = FALSE}.}
#' }
#'
#' @examples
#' library(terra)
#'
#' # Small toy example -- runs in well under 5 seconds
#' r <- rast(nrows = 2, ncols = 2, nlyrs = 4)
#' values(r) <- c(
#'   2, 3, 5, 7,
#'   3, 4, 6, 8,
#'   2, 3, 5, 7,
#'   3, 4, 6, 8
#' )
#' names(r) <- c(
#'   "TTR_70", "TTR_80", "TTR_90",
#'   "TTR_Threshold_Range"
#' )
#'
#' result <- classify_recovery_typology(
#'   ttr_stack = r,
#'   fast_TTR = 3,
#'   unc_threshold = 2,
#'   plot_map = FALSE,
#'   plot_stats = FALSE,
#'   area = TRUE
#' )
#'
#' result$Stats
#'
#' \dontshow{
#' # Internal checks (not shown to users)
#' stopifnot(is.list(result))
#' stopifnot(inherits(result$Recovery_Type, "SpatRaster"))
#' stopifnot(is.data.frame(result$Stats))
#' stopifnot(nrow(result$Stats) == 4L)
#' stopifnot(all(c("Count", "Percent", "Area_ha") %in%
#'                 names(result$Stats)))
#' stopifnot(all(as.character(result$Stats$Recovery_Type) %in%
#'   c("Fast & stable recovery",
#'     "Fast but threshold-sensitive recovery",
#'     "Slow & stable recovery",
#'     "Slow & threshold-sensitive recovery")))
#' }
#'
#' \donttest{
#' # Statistics plot
#' result_plot <- classify_recovery_typology(
#'   ttr_stack = r,
#'   fast_TTR = 3,
#'   unc_threshold = 2,
#'   plot_map = FALSE,
#'   plot_stats = TRUE,
#'   area = FALSE
#' )
#'
#' result_plot$Stats_Plot
#' }
#'
#' @importFrom terra nlyr app as.int values res
#' @importFrom ggplot2 ggplot aes geom_col geom_text scale_fill_manual labs
#'   theme_classic theme element_text
#' @importFrom dplyr filter group_by summarise mutate n
#' @export
classify_recovery_typology <- function(
  ttr_stack,
  fast_TTR = 3,
  unc_threshold = 2,
  plot_map = TRUE,
  plot_stats = TRUE,
  area = TRUE
) {

  # ---- Input validation ----

  if (!inherits(ttr_stack, "SpatRaster")) {
    stop("'ttr_stack' must be a terra::SpatRaster")
  }

  if (!is.numeric(fast_TTR) ||
      length(fast_TTR) != 1L ||
      is.na(fast_TTR) ||
      fast_TTR < 0) {
    stop("'fast_TTR' must be a single non-negative numeric value")
  }

  if (!is.numeric(unc_threshold) ||
      length(unc_threshold) != 1L ||
      is.na(unc_threshold) ||
      unc_threshold < 0) {
    stop("'unc_threshold' must be a single non-negative numeric value")
  }

  if (!is.logical(plot_map) || length(plot_map) != 1L || is.na(plot_map)) {
    stop("'plot_map' must be TRUE or FALSE")
  }

  if (!is.logical(plot_stats) ||
      length(plot_stats) != 1L ||
      is.na(plot_stats)) {
    stop("'plot_stats' must be TRUE or FALSE")
  }

  if (!is.logical(area) || length(area) != 1L || is.na(area)) {
    stop("'area' must be TRUE or FALSE")
  }

  if (!"TTR_Threshold_Range" %in% names(ttr_stack)) {
    stop(
      "ttr_stack must contain a 'TTR_Threshold_Range' layer. ",
      "Call compute_recovery_time_uncertainty(return_band = TRUE)."
    )
  }

  ttr_layers <- grep(
    "^TTR_[0-9]+(?:\\.[0-9]+)?$",
    names(ttr_stack),
    value = TRUE
  )

  if (length(ttr_layers) < 2L) {
    stop("At least two threshold-specific TTR layers are required.")
  }

  # ---- Mean TTR ----

  TTR_mean <- terra::app(
    ttr_stack[[ttr_layers]],
    fun = function(x) {
      if (all(is.na(x))) {
        return(NA_real_)
      }
      mean(x, na.rm = TRUE)
    }
  )

  names(TTR_mean) <- "TTR_Mean"

  # ---- Combine mean TTR and threshold range ----

  TTR_combined <- c(
    TTR_mean,
    ttr_stack[["TTR_Threshold_Range"]]
  )

  names(TTR_combined) <- c(
    "TTR_Mean",
    "TTR_Threshold_Range"
  )

  # ---- Classify recovery typology ----

  Recovery_type <- terra::app(
    TTR_combined,
    fun = function(x) {

      TTRm <- x[1]
      range_TTR <- x[2]

      if (is.na(TTRm) || is.na(range_TTR)) {
        return(NA_integer_)
      }

      if (TTRm <= fast_TTR && range_TTR <= unc_threshold) return(1L)
      if (TTRm <= fast_TTR && range_TTR > unc_threshold)  return(2L)
      if (TTRm > fast_TTR  && range_TTR <= unc_threshold) return(3L)

      return(4L)
    }
  )

  Recovery_type <- terra::as.int(Recovery_type)
  names(Recovery_type) <- "Recovery_Type"

  # ---- Class labels ----

  class_labels <- c(
    "1" = "Fast & stable recovery",
    "2" = "Fast but threshold-sensitive recovery",
    "3" = "Slow & stable recovery",
    "4" = "Slow & threshold-sensitive recovery"
  )

  levels(Recovery_type) <- data.frame(
    value = 1:4,
    label = unname(class_labels)
  )

  # ---- Recovery typology map ----

  if (plot_map) {
    terra::plot(
      Recovery_type,
      col = c(
        "#238b45",
        "#74c476",
        "#fdae6b",
        "#d7301f"
      ),
      main = "Post-disturbance Recovery Typology"
    )
  }

  # ---- Statistics ----

  vals <- terra::values(Recovery_type, mat = FALSE)

  df <- data.frame(
    Recovery_Type = vals
  )

  df <- dplyr::filter(df, !is.na(Recovery_Type))

  stats <- df |>
    dplyr::group_by(Recovery_Type) |>
    dplyr::summarise(
      Count = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      Percent = 100 * Count / sum(Count)
    )

  # ---- Add missing classes if necessary ----

  all_classes <- data.frame(
    Recovery_Type = 1:4
  )

  stats <- merge(
    all_classes,
    stats,
    by = "Recovery_Type",
    all.x = TRUE,
    sort = TRUE
  )

  stats$Count[is.na(stats$Count)] <- 0L
  stats$Percent[is.na(stats$Percent)] <- 0

  # ---- Convert class codes to labeled factor ----

  level_order <- unname(class_labels[c("1", "2", "3", "4")])

  stats$Recovery_Type <- factor(
    class_labels[as.character(stats$Recovery_Type)],
    levels = level_order
  )

  # ---- Area calculation ----

  if (area) {

    px_area <- prod(terra::res(Recovery_type))

    stats <- dplyr::mutate(
      stats,
      Area_ha = Count * px_area / 1e4
    )
  }

  # ---- Statistics plot ----

  stats_plot <- NULL

  if (plot_stats) {

    stats_plot <- ggplot2::ggplot(
      stats,
      ggplot2::aes(
        x = Recovery_Type,
        y = Percent,
        fill = Recovery_Type
      )
    ) +
      ggplot2::geom_col(color = "black") +
      ggplot2::geom_text(
        ggplot2::aes(
          label = paste0(round(Percent, 1), "%")
        ),
        vjust = -0.5,
        size = 3.5
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          "Fast & stable recovery"                  = "#238b45",
          "Fast but threshold-sensitive recovery"   = "#74c476",
          "Slow & stable recovery"                  = "#fdae6b",
          "Slow & threshold-sensitive recovery"     = "#d7301f"
        )
      ) +
      ggplot2::labs(
        title = "Recovery Typology Statistics",
        x = "Recovery Type",
        y = "Percentage (%)"
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1
        ),
        legend.position = "none",
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          size = 14,
          face = "bold"
        )
      )
  }

  # ---- Return ----

  list(
    TTR_Mean = TTR_mean,
    TTR_Combined = TTR_combined,
    Recovery_Type = Recovery_type,
    Stats = stats,
    Stats_Plot = stats_plot
  )
}
