## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>"
)


## -----------------------------------------------------------------------------
library(RecoverMetricsR)
library(terra)

# --- Small synthetic dataset -------------------------------------
set.seed(42)

r_stack_processed <- rast(nrows = 20, ncols = 20, nlyrs = 8)
names(r_stack_processed) <- paste0("NBR_", 1998:2005)
years <- 1998:2005

# Simulate a burn in 2000
vals <- matrix(
  runif(20 * 20 * 8, 0.4, 0.7),
  ncol = 20 * 20
)
vals[, ] <- vals
for (i in 1:8) {
  values(r_stack_processed[[i]]) <- vals[i, ]
}
values(r_stack_processed[["NBR_2000"]]) <- values(r_stack_processed[["NBR_2000"]]) * 0.3

Fire_scars <- vect(
  "POLYGON ((2 2, 18 2, 18 18, 2 18, 2 2))",
  crs = crs(r_stack_processed)
)

# --- 1. Compute pixel-wise recovery metrics ----------------------
metrics_stack <- compute_recovery_metrics(
  raster_stack = r_stack_processed,
  years        = years,
  event_year   = 2000,
  pre_window   = 2,
  post_window  = 4,
  clamp        = TRUE,
  metrics      = c("DM", "Slope", "YrYr_early",
                   "YrYr_late", "Recovery_Rate_Change",
                   "RRI", "R80P", "YrYr"),
  mask_layer   = Fire_scars
)

print(metrics_stack)

# --- 2. Compute TTR and uncertainty ------------------------------
TTR_uncertainty <- compute_recovery_time_uncertainty(
  raster_stack = r_stack_processed,
  years        = years,
  event_year   = 2000,
  thresholds   = c(0.7, 0.8, 0.9),
  persistence  = 1,
  mask_layer   = Fire_scars,
  return_band  = TRUE
)

print(TTR_uncertainty)

# --- 3. Classify recovery typology -------------------------------
Recovery_results <- classify_recovery_typology(
  ttr_stack  = TTR_uncertainty,
  fast_TTR   = 3,
  unc_threshold = 2,
  plot_map   = FALSE,
  plot_stats = FALSE,
  area       = TRUE
)

print(Recovery_results$Stats)

