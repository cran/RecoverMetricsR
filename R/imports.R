# imports.R
# Ensures R CMD check recognizes imported functions

#' @importFrom stats sd quantile na.omit AIC BIC median nls nls.control
#'   predict complete.cases vcov qnorm coef lm
#' @importFrom grid unit
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_vline geom_errorbar
#'   geom_label geom_col geom_text geom_bar labs theme_classic theme_minimal
#'   theme element_text scale_fill_manual scale_fill_brewer
#' @importFrom dplyr filter group_by summarise mutate n bind_rows
#' @importFrom terra rast app mask vect nlyr global values ifel as.int res
#' @importFrom utils globalVariables
NULL
