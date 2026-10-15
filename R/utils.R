#' @title Format a Fitted Recovery Equation as a String
#' @description Internal helper for \code{fit_recovery_model()}.
#'   Not exported.
#' @param model Character. One of "exponential", "logistic", "gompertz".
#' @param coef Named numeric vector of model coefficients
#'   (must contain \code{A}, \code{k}, and, for logistic/Gompertz,
#'   \code{x0}).
#' @return A character string representing the fitted equation.
#' @keywords internal
#' @noRd
.format_equation <- function(model, coef) {
  if (model == "exponential") {
    sprintf("y = %.2f (1 - exp(-%.2f t))", coef["A"], coef["k"])
  } else if (model == "logistic") {
    sprintf("y = %.2f / (1 + exp(-%.2f (t - %.2f)))",
            coef["A"], coef["k"], coef["x0"])
  } else if (model == "gompertz") {
    sprintf("y = %.2f exp(-exp(-%.2f (t - %.2f)))",
            coef["A"], coef["k"], coef["x0"])
  } else {
    stop("Unknown model: ", model)
  }
}
