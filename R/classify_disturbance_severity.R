#' Classify Spectral Disturbance and Recovery Across Post-Disturbance Years
#'
#' @description
#' Classifies the magnitude and direction of spectral change between a
#' pre-disturbance reference layer and multiple post-disturbance layers.
#' Positive values of the pre-minus-post difference indicate a decrease
#' in the spectral index, whereas negative values indicate an increase.
#'
#' @param raster_stack A \code{terra::SpatRaster} containing a spectral
#'   index time series.
#' @param pre_year Layer index, layer name, or year corresponding to the
#'   pre-disturbance reference condition.
#' @param event_year Layer index, layer name, or year corresponding to the
#'   disturbance year. The argument is validated for workflow consistency
#'   but the event layer itself is not used in the post-disturbance
#'   classification.
#' @param post_years Vector of layer indices, layer names, or years
#'   corresponding to post-disturbance observations. The event year is
#'   normally excluded so that all layers describe spectral change
#'   relative to the pre-disturbance reference.
#' @param mask_layer Optional spatial mask. Can be a
#'   \code{terra::SpatRaster}, \code{terra::SpatVector}, or \code{sf}
#'   object used to restrict the analysis spatially.
#' @param plot Logical. If \code{TRUE}, produces a plot of spectral-change
#'   class distributions when \code{combine = TRUE}.
#' @param combine Logical. If \code{TRUE}, combines statistics from all
#'   post-disturbance years into a single data frame and plot.
#'
#' @details
#' For each post-disturbance layer, spectral change is calculated as:
#'
#' \deqn{
#' S_{pre} - S_{post}
#' }
#'
#' where \code{S_pre} is the pre-disturbance spectral index value and
#' \code{S_post} is the value observed in the corresponding
#' post-disturbance year.
#'
#' The resulting difference is assigned to seven classes ranging from
#' high spectral disturbance to high spectral recovery. The classes
#' describe the magnitude and direction of spectral change relative to
#' the pre-disturbance reference.
#'
#' A near-zero difference indicates a spectral condition close to the
#' pre-fire reference. Positive values indicate that the post-disturbance
#' spectral index remains below the pre-fire reference, whereas negative
#' values indicate that the post-disturbance spectral index exceeds the
#' pre-fire reference.
#'
#' The classification thresholds are:
#' \itemize{
#'   \item High spectral disturbance: difference > 0.66
#'   \item High-moderate spectral disturbance: 0.44 < difference <= 0.66
#'   \item Low-moderate spectral disturbance: 0.27 < difference <= 0.44
#'   \item Low spectral disturbance: 0.10 < difference <= 0.27
#'   \item Near pre-fire spectral condition: -0.10 < difference <= 0.10
#'   \item Low spectral recovery: -0.25 < difference <= -0.10
#'   \item High spectral recovery: difference <= -0.25
#' }
#'
#' The \code{Near pre-fire spectral condition} class represents a
#' spectral state close to the pre-fire reference, within approximately
#' +/- 0.10 spectral-index units. It does not imply complete ecological
#' or structural recovery.
#'
#' The \code{Low spectral recovery} and \code{High spectral recovery}
#' classes indicate that the post-disturbance spectral index exceeds
#' the pre-fire reference. They describe spectral recovery relative to
#' the pre-fire reference and should not be interpreted as direct
#' measurements of ecological regeneration or structural recovery.
#'
#' These classes describe spectral change and should not be interpreted
#' as direct measurements of ecological, structural, or combustion
#' severity without independent validation.
#'
#' Layer selection accepts three forms: a numeric index in
#' \code{1..nlyr(raster_stack)}, an exact layer name, or a year value
#' that appears within the layer names. Year matching uses a
#' digit-boundary search, so a year like \code{2000} will match
#' \code{"NBR_2000"} but not \code{"NBR_20000"}.
#'
#' By convention, the event year is excluded from \code{post_years} so
#' that all layers describe post-fire spectral change relative to the
#' pre-disturbance reference. Users who wish to map immediate fire
#' severity relative to the pre-disturbance reference can include the
#' event-year layer in \code{post_years} explicitly. In that case, the
#' resulting classes describe the immediate spectral disturbance
#' associated with the event.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{severity_stack}{A \code{terra::SpatRaster} with one layer per
#'     value in \code{post_years}, named after the corresponding year
#'     and containing integer spectral-change classes.}
#'   \item{stats_by_year}{A named list of data frames containing pixel
#'     counts and percentages for each class and post-disturbance year.
#'     Names match the layer names in \code{severity_stack}.}
#'   \item{combined_stats}{A data frame combining spectral-change class
#'     statistics across post-disturbance years.}
#'   \item{plot}{A \code{ggplot2} object showing the distribution of
#'     spectral-change classes, or \code{NULL} if plotting is disabled.}
#' }
#'
#' @examples
#' library(terra)
#'
#' # Small toy example -- runs in well under 5 seconds
#' r <- rast(
#'   nrows = 2,
#'   ncols = 2,
#'   nlyrs = 5
#' )
#'
#' values(r) <- rep(
#'   c(0.70, 0.30, 0.50, 0.55, 0.60),
#'   each = 4
#' )
#'
#' names(r) <- paste0(
#'   "NBR_",
#'   2020:2024
#' )
#'
#' # Event year excluded from post_years
#' result <- classify_disturbance_severity(
#'   raster_stack = r,
#'   pre_year = "NBR_2020",
#'   event_year = "NBR_2021",
#'   post_years = c(
#'     "NBR_2022",
#'     "NBR_2023",
#'     "NBR_2024"
#'   ),
#'   plot = FALSE
#' )
#'
#' result$combined_stats
#'
#' # Years are also accepted
#' result2 <- classify_disturbance_severity(
#'   raster_stack = r,
#'   pre_year = 2020,
#'   event_year = 2021,
#'   post_years = 2022:2024,
#'   plot = FALSE
#' )
#'
#' identical(
#'   result$combined_stats$Count,
#'   result2$combined_stats$Count
#' )
#'
#' \dontshow{
#' stopifnot(is.list(result))
#' stopifnot(inherits(
#'   result$severity_stack,
#'   "SpatRaster"
#' ))
#'
#' stopifnot(
#'   terra::nlyr(result$severity_stack) == 3L
#' )
#'
#' stopifnot(
#'   identical(
#'     names(result$severity_stack),
#'     c(
#'       "NBR_2022",
#'       "NBR_2023",
#'       "NBR_2024"
#'     )
#'   )
#' )
#'
#' stopifnot(
#'   identical(
#'     names(result2$severity_stack),
#'     c(
#'       "2022",
#'       "2023",
#'       "2024"
#'     )
#'   )
#' )
#'
#' stopifnot(
#'   is.data.frame(result$combined_stats)
#' )
#'
#' stopifnot(
#'   all(
#'     result$combined_stats$Class %in%
#'       c(
#'         "High spectral disturbance",
#'         "High-moderate spectral disturbance",
#'         "Low-moderate spectral disturbance",
#'         "Low spectral disturbance",
#'         "Near pre-fire spectral condition",
#'         "Low spectral recovery",
#'         "High spectral recovery"
#'       )
#'   )
#' )
#'
#' stopifnot(
#'   identical(
#'     result$combined_stats$Count,
#'     result2$combined_stats$Count
#'   )
#' )
#' }
#'
#' \donttest{
#' # Combined statistics plot
#' result_plot <- classify_disturbance_severity(
#'   raster_stack = r,
#'   pre_year = "NBR_2020",
#'   event_year = "NBR_2021",
#'   post_years = c(
#'     "NBR_2022",
#'     "NBR_2023",
#'     "NBR_2024"
#'   ),
#'   plot = TRUE,
#'   combine = TRUE
#' )
#'
#' result_plot$plot
#' }
#'
#' @importFrom terra values rast mask nlyr vect ifel levels as.int coltab
#'   same.crs resample res ext project
#' @importFrom ggplot2 ggplot aes geom_bar labs theme_classic theme
#'   element_text scale_fill_brewer
#' @importFrom dplyr group_by summarise n bind_rows
#' @export
classify_disturbance_severity <- function(
  raster_stack,
  pre_year,
  event_year,
  post_years,
  mask_layer = NULL,
  plot = TRUE,
  combine = TRUE
) {

  # ---- Input validation --------------------------------------------

  if (!inherits(raster_stack, "SpatRaster")) {
    stop("'raster_stack' must be a terra::SpatRaster")
  }

  if (length(post_years) < 1L) {
    stop("'post_years' must contain at least one layer")
  }

  if (!is.logical(plot) ||
      length(plot) != 1L ||
      is.na(plot)) {
    stop("'plot' must be TRUE or FALSE")
  }

  if (!is.logical(combine) ||
      length(combine) != 1L ||
      is.na(combine)) {
    stop("'combine' must be TRUE or FALSE")
  }

  # ---- Helper: resolve a layer specifier to a raster index --------

  get_layer_idx <- function(r, l) {

    # Case 1: exact layer name
    if (is.character(l)) {

      idx <- match(l, names(r))

      if (is.na(idx)) {
        stop("Layer name not found: ", l)
      }

      return(idx)
    }

    # Case 2: numeric -- index or year
    if (is.numeric(l) &&
        length(l) == 1L &&
        !is.na(l)) {

      # 2a: valid layer index
      if (l == as.integer(l) &&
          l >= 1L &&
          l <= terra::nlyr(r)) {

        return(as.integer(l))
      }

      # 2b: treat as a year and search layer names
      idx <- grep(
        pattern = paste0(
          "(^|[^0-9])",
          l,
          "([^0-9]|$)"
        ),
        names(r)
      )

      if (length(idx) == 1L) {
        return(idx)
      }

      if (length(idx) > 1L) {
        stop(
          "Ambiguous year '", l,
          "': matches multiple layer names (",
          paste(names(r)[idx], collapse = ", "),
          "). Please supply a layer name instead."
        )
      }

      stop(
        "Cannot resolve '", l,
        "' to a raster layer. ",
        "It is not a valid index (1..",
        terra::nlyr(r),
        ") and does not match any layer name."
      )
    }

    stop(
      "Layer must be a valid layer index, layer name, or year."
    )
  }

  get_layer <- function(r, l) {
    r[[get_layer_idx(r, l)]]
  }

  # ---- Resolve reference and event layers -------------------------

  pre <- get_layer(
    raster_stack,
    pre_year
  )

  # Validated for workflow consistency.
  # The event layer is not used unless explicitly included in
  # post_years.
  event <- get_layer(
    raster_stack,
    event_year
  )

  # Avoid an unused-variable NOTE in R CMD check.
  invisible(event)

  # ---- Optional mask -----------------------------------------------

  mask_r <- NULL

  if (!is.null(mask_layer)) {

    mask_r <- if (inherits(mask_layer, "SpatRaster")) {

      mask_layer

    } else if (inherits(mask_layer, "SpatVector")) {

      mask_layer

    } else if (inherits(mask_layer, "sf")) {

      terra::vect(mask_layer)

    } else {

      stop(
        "'mask_layer' must be a terra::SpatRaster, ",
        "terra::SpatVector, or sf object"
      )
    }

    if (inherits(mask_r, "SpatRaster")) {

      needs_resample <-
        !terra::same.crs(
          raster_stack,
          mask_r
        ) ||
        !isTRUE(
          all.equal(
            terra::res(raster_stack),
            terra::res(mask_r)
          )
        ) ||
        !isTRUE(
          all.equal(
            as.vector(terra::ext(raster_stack)),
            as.vector(terra::ext(mask_r))
          )
        )

      if (needs_resample) {

        mask_r <- terra::resample(
          mask_r,
          raster_stack,
          method = "near"
        )
      }

    } else {

      if (!terra::same.crs(
        raster_stack,
        mask_r
      )) {

        mask_r <- terra::project(
          mask_r,
          raster_stack
        )
      }
    }

    pre <- terra::mask(
      pre,
      mask_r
    )
  }

  # ---- Spectral-change class definitions --------------------------

  class_defs <- data.frame(
    code = 1:7,

    label = c(
      "High spectral disturbance",
      "High-moderate spectral disturbance",
      "Low-moderate spectral disturbance",
      "Low spectral disturbance",
      "Near pre-fire spectral condition",
      "Low spectral recovery",
      "High spectral recovery"
    ),

    lower = c(
      0.66,
      0.44,
      0.27,
      0.10,
      -0.10,
      -0.25,
      -Inf
    ),

    upper = c(
      Inf,
      0.66,
      0.44,
      0.27,
      0.10,
      -0.10,
      -0.25
    ),

    stringsAsFactors = FALSE
  )

  severity_layers <- list()
  stats_by_year <- list()
  combined_stats <- data.frame()

  # ---- Process each post-disturbance year -------------------------

  for (yr in post_years) {

    idx <- get_layer_idx(
      raster_stack,
      yr
    )

    # Output layer name follows the supplied post_years value.
    lyr_name <- as.character(yr)

    post <- raster_stack[[idx]]

    if (!is.null(mask_r)) {

      post <- terra::mask(
        post,
        mask_r
      )
    }

    # Pre-disturbance minus post-disturbance spectral index
    sev <- pre - post

    # ---- Classify spectral change -------------------------------

    class_r <- terra::ifel(
      is.na(sev),
      NA_integer_,
      0L
    )

    for (i in seq_len(
      nrow(class_defs)
    )) {

      cond <- (
        sev > class_defs$lower[i] &
          sev <= class_defs$upper[i]
      )

      class_r <- terra::ifel(
        cond,
        class_defs$code[i],
        class_r
      )
    }

    class_r <- terra::as.int(
      class_r
    )

    names(class_r) <- lyr_name

    severity_layers[[lyr_name]] <- class_r

    # ---- Statistics ----------------------------------------------

    vals <- terra::values(
      class_r,
      mat = FALSE
    )

    df <- data.frame(
      Class = factor(
        vals,
        levels = class_defs$code,
        labels = class_defs$label
      )
    )

    df <- df[
      !is.na(df$Class),
      ,
      drop = FALSE
    ]

    if (nrow(df) > 0L) {

      stats <- df |>
        dplyr::group_by(
          Class,
          .drop = FALSE
        ) |>
        dplyr::summarise(
          Count = dplyr::n(),
          .groups = "drop"
        ) |>
        dplyr::mutate(
          Percent = 100 * Count / sum(Count)
        )

    } else {

      stats <- data.frame(
        Class = factor(
          class_defs$label,
          levels = class_defs$label
        ),
        Count = rep(
          0L,
          length(class_defs$label)
        ),
        Percent = rep(
          0,
          length(class_defs$label)
        )
      )
    }

    stats$Year <- lyr_name

    stats_by_year[[lyr_name]] <- stats

    combined_stats <- dplyr::bind_rows(
      combined_stats,
      stats
    )
  }

  # ---- Build multi-layer raster -----------------------------------

  severity_stack <- terra::rast(
    severity_layers
  )

  # ---- Attach categorical levels ---------------------------------

  lvl_df <- data.frame(
    value = class_defs$code,
    label = class_defs$label
  )

  levels(severity_stack) <- rep(
    list(lvl_df),
    terra::nlyr(severity_stack)
  )

  # ---- Attach colors ----------------------------------------------

  ct <- data.frame(
    value = class_defs$code,
    color = c(
      "#8c2d04",
      "#d94801",
      "#f16913",
      "#fdae6b",
      "#bdbdbd",
      "#74c476",
      "#238b45"
    )
  )

  terra::coltab(
    severity_stack
  ) <- rep(
    list(ct),
    terra::nlyr(severity_stack)
  )

  names(severity_stack) <- vapply(
    post_years,
    as.character,
    character(1)
  )

  # ---- Plot -------------------------------------------------------

  p <- NULL

  if (plot &&
      combine &&
      nrow(combined_stats) > 0L) {

    combined_stats$Year <- factor(
      combined_stats$Year,
      levels = unique(
        combined_stats$Year
      )
    )

    p <- ggplot2::ggplot(
      combined_stats,
      ggplot2::aes(
        x = Year,
        y = Percent,
        fill = Class
      )
    ) +
      ggplot2::geom_bar(
        stat = "identity"
      ) +
      ggplot2::labs(
        x = "Post-disturbance Year",
        y = "Percentage (%)",
        title = "Evolution of Spectral Disturbance and Recovery Classes"
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          size = 12
        ),
        axis.text.y = ggplot2::element_text(
          size = 12
        ),
        plot.title = ggplot2::element_text(
          size = 14,
          face = "bold"
        )
      ) +
      ggplot2::scale_fill_brewer(
        palette = "Set2"
      )
  }

  # ---- Return -----------------------------------------------------

  list(
    severity_stack = severity_stack,
    stats_by_year = stats_by_year,
    combined_stats = combined_stats,
    plot = p
  )
}
