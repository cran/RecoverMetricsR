#' @title Preprocess Raster Data for Spectral Recovery Metrics
#'
#' @description
#' Aggregates monthly raster stacks to yearly values and optionally applies
#' a spatial mask. The resulting yearly stack is intended as input to
#' \code{\link{compute_recovery_metrics}}.
#'
#' @param raster_stack A \code{terra::SpatRaster}. Yearly, monthly, or
#'   pre-aggregated stack.
#' @param data_type Character. Either \code{"yearly"} or \code{"monthly"}.
#' @param months_to_aggregate Numeric vector of months to aggregate. Used
#'   only when \code{data_type == "monthly"}.
#' @param years Numeric vector of years corresponding to raster layers.
#' @param mask_layer Optional spatial mask. Can be a
#'   \code{terra::SpatRaster}, \code{terra::SpatVector}, or an \code{sf}
#'   object. Vector masks are reprojected to the raster CRS if needed, and
#'   raster masks are resampled to the raster grid if their CRS,
#'   resolution, or extent differ.
#' @param fun Character. Aggregation function, either \code{"mean"} or
#'   \code{"median"}. Default is \code{"mean"}.
#'
#' @details
#' When \code{data_type = "monthly"}, layer names must contain a year and
#' month in the form \code{"YYYY_MM"} or \code{"YYYY-MM"} (for example,
#' \code{"2020_01"}). Layers matching each requested year and the months
#' in \code{months_to_aggregate} are aggregated with \code{fun} to produce
#' a single yearly raster. Years with no matching layers are silently
#' skipped; an error is raised if no year matches at all.
#'
#' The mask is applied after aggregation. Raster masks are aligned to
#' \code{raster_stack} by resampling if their grid differs.
#'
#' @return A \code{terra::SpatRaster} containing one layer per year,
#'   ready for metric computation.
#'
#' @examples
#' library(terra)
#'
#' # Yearly input -- no aggregation required
#' r <- rast(nrows = 2, ncols = 2, nlyrs = 3)
#' values(r) <- rep(c(0.5, 0.6, 0.7), each = 4)
#' names(r) <- c("2020", "2021", "2022")
#'
#' yearly <- preprocess_raster(
#'   raster_stack = r,
#'   data_type = "yearly",
#'   years = 2020:2022
#' )
#' yearly
#'
#' # Mask with a SpatVector (no sf dependency)
#' v <- vect(
#'   "POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))",
#'   crs = crs(r)
#' )
#' masked <- preprocess_raster(
#'   raster_stack = r,
#'   data_type = "yearly",
#'   years = 2020:2022,
#'   mask_layer = v
#' )
#' masked
#'
#' # Monthly input -- aggregate selected months to yearly means
#' rm <- rast(nrows = 2, ncols = 2, nlyrs = 6)
#' values(rm) <- rep(c(0.5, 0.6, 0.7, 0.4, 0.5, 0.6), each = 4)
#' names(rm) <- c(
#'   "2020_06", "2020_07", "2020_08",
#'   "2021_06", "2021_07", "2021_08"
#' )
#'
#' agg <- preprocess_raster(
#'   raster_stack = rm,
#'   data_type = "monthly",
#'   months_to_aggregate = 6:8,
#'   years = 2020:2021,
#'   fun = "mean"
#' )
#' agg
#'
#' \dontshow{
#' # Internal checks (not shown to users)
#' stopifnot(inherits(yearly, "SpatRaster"))
#' stopifnot(terra::nlyr(yearly) == 3L)
#' stopifnot(inherits(masked, "SpatRaster"))
#' stopifnot(terra::nlyr(masked) == 3L)
#' stopifnot(inherits(agg, "SpatRaster"))
#' stopifnot(terra::nlyr(agg) == 2L)
#' }
#'
#' \donttest{
#' # Example with a spatial mask (sf polygon)
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
#' masked_sf <- preprocess_raster(
#'   raster_stack = r,
#'   data_type = "yearly",
#'   years = 2020:2022,
#'   mask_layer = p
#' )
#' masked_sf
#' }
#'
#' @importFrom terra app rast mask vect nlyr same.crs res resample project ext
#' @export
preprocess_raster <- function(raster_stack,
                              data_type = c("yearly", "monthly"),
                              months_to_aggregate = NULL,
                              years,
                              mask_layer = NULL,
                              fun = "mean") {

  data_type <- match.arg(data_type)

  if (!fun %in% c("mean", "median"))
    stop("fun must be 'mean' or 'median'")

  if (!inherits(raster_stack, "SpatRaster"))
    stop("raster_stack must be a terra::SpatRaster")

  if (data_type == "monthly") {

    if (is.null(months_to_aggregate))
      stop("For monthly data, months_to_aggregate must be provided")

    lyr_names <- names(raster_stack)

    m <- regexec("(\\d{4})[-_](\\d{1,2})", lyr_names)
    parsed <- regmatches(lyr_names, m)

    layer_year <- as.numeric(vapply(
      parsed,
      function(p) if (length(p) >= 2) p[2] else NA_character_,
      character(1)
    ))

    layer_month <- as.numeric(vapply(
      parsed,
      function(p) if (length(p) >= 3) p[3] else NA_character_,
      character(1)
    ))

    if (anyNA(layer_year) || anyNA(layer_month))
      stop(
        "Layer names must contain a 'YYYY_MM' or 'YYYY-MM' pattern ",
        "(e.g., '2020_01')."
      )

    fun_obj <- match.fun(fun)
    yearly_rasters <- list()

    for (y in years) {

      sel <- which(layer_year == y & layer_month %in% months_to_aggregate)

      if (length(sel) == 0) next

      selected_layers <- raster_stack[[sel]]

      yearly_rasters[[as.character(y)]] <-
        terra::app(selected_layers, fun_obj, na.rm = TRUE)
    }

    if (length(yearly_rasters) == 0)
      stop("No layers matched the provided years/months")

    raster_stack <- terra::rast(yearly_rasters)
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

  raster_stack
}
