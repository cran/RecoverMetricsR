########################################################################################
# ===========================================================
# global.R
# Register global variables to satisfy R CMD check
# ===========================================================

utils::globalVariables(
  c(
    # Recovery trajectory plotting (fit_recovery_model)
    "Year", "Mean", "SD", "Mask",
    # Model overlay
    "Fitted", "x", "y", "label",
    # Recovery typology statistics
    "Recovery_Type", "Count", "Percent", "Area_ha",
    # Disturbance severity
    "Class"
  )
)
