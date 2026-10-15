#' @title Fit Core Recovery Model
#'
#' @description
#' Fits a parametric recovery model (exponential, logistic, or Gompertz)
#' to a mean post-disturbance spectral trajectory using bounded nonlinear
#' least squares. Returns \code{NULL} if the model fails to converge or
#' if the fitted parameters are invalid.
#'
#' @param t Numeric vector of time since disturbance.
#' @param y Numeric vector of mean spectral index values corresponding to
#'   \code{t}.
#' @param model Character. One of \code{"exponential"}, \code{"logistic"},
#'   or \code{"gompertz"}.
#'
#' @details
#' The exponential model represents recovery as an asymptotic process:
#' \deqn{y(t) = A(1 - exp(-k t))}
#' where \code{A} is the asymptotic recovery level and \code{k} is the
#' recovery-rate parameter. Exponential formulations have been used to
#' characterize post-fire vegetation recovery from satellite time series
#' (Gouveia et al., 2010) and post-fire NBR recovery
#' (Chompuchan and Lin, 2017).
#'
#' The logistic model is defined as:
#' \deqn{y(t) = A / (1 + exp(-k(t - x0)))}
#' where \code{A} is the asymptotic level, \code{k} is the rate parameter,
#' and \code{x0} is the inflection point.
#'
#' The Gompertz model is defined as:
#' \deqn{y(t) = A exp(-exp(-k(t - x0)))}
#' where \code{A} is the asymptotic level, \code{k} is the rate parameter,
#' and \code{x0} is the inflection point. Gompertz functions have been
#' used to characterize vegetation biomass recovery following fire.
#'
#' The logistic and Gompertz formulations are provided as alternative
#' nonlinear representations of recovery trajectories. Model selection
#' should be based on the characteristics and fit of the observed
#' recovery trajectory.
#'
#' Bounds on fitted parameters are \code{A} in \code{[0, 2]},
#' \code{k} in \code{[0.001, 5]}, and \code{x0} in \code{[0, 30]}.
#' These are appropriate for normalized difference spectral indices
#' such as NBR or NDVI in \code{[-1, 1]}. If the input \code{y} falls
#' outside these bounds, the fit will fail and \code{NULL} is returned.
#'
#' The function returns \code{NULL} (rather than raising an error) when
#' the model cannot be fitted or when fitted parameters fall outside
#' plausible ranges, so that it can be used safely inside loops.
#' Callers should check \code{is.null(fit)} before accessing list
#' elements.
#'
#' The \code{$fit} element of the returned list is the raw \code{nls}
#' object, which supports \code{predict()} and \code{summary()}.
#' Convenience extracts (\code{$coefficients}, \code{$vcov},
#' \code{$sigma}) are also provided.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{model}{The fitted model name.}
#'   \item{coefficients}{Estimated model coefficients.}
#'   \item{sigma}{Residual standard error estimated from the fitted model.}
#'   \item{vcov}{Variance-covariance matrix of the estimated coefficients,
#'     when available.}
#'   \item{data}{Data used for model fitting after removal of non-finite
#'     observations.}
#'   \item{fit}{The fitted \code{nls} model object.}
#' }
#' If the model cannot be fitted or the fitted parameters are invalid,
#' \code{NULL} is returned.
#'
#' @references
#' Chompuchan, C., & Lin, C.-Y. (2017). Assessment of forest recovery
#' at the deciduous forest in a complex terrain using Landsat time
#' series. \emph{Ecological Indicators}, 79, 249--259.
#' \doi{10.1016/j.ecolind.2017.04.038}
#'
#' Gouveia, C., DaCamara, C. C., & Trigo, R. M. (2010). Post-fire
#' vegetation recovery in Portugal based on spot/vegetation data.
#' \emph{Natural Hazards and Earth System Sciences}, 10(4), 673--684.
#' \doi{10.5194/nhess-10-673-2010}
#'
#' @examples
#' # Simulated exponential recovery trajectory
#' t <- 0:5
#' y <- 0.70 * (1 - exp(-0.6 * t))
#'
#' fit <- fit_recovery_model_core(
#'   t = t,
#'   y = y,
#'   model = "exponential"
#' )
#'
#' fit$coefficients
#'
#' \dontshow{
#' # Internal checks (not shown to users)
#' stopifnot(!is.null(fit))
#' stopifnot(is.list(fit))
#' stopifnot("A" %in% names(fit$coefficients))
#' stopifnot("k" %in% names(fit$coefficients))
#' stopifnot(is.finite(fit$coefficients["A"]))
#' stopifnot(is.finite(fit$coefficients["k"]))
#'
#' # Degenerate inputs return NULL
#' stopifnot(is.null(fit_recovery_model_core(numeric(0), numeric(0))))
#' stopifnot(is.null(fit_recovery_model_core(1:3, c(0.5, 0.6, 0.7))))
#' stopifnot(is.null(fit_recovery_model_core(
#'   c(1, 1, 1, 1), c(0.5, 0.6, 0.7, 0.8)
#' )))
#' }
#'
#' \donttest{
#' # Logistic and Gompertz alternatives
#' t2 <- 0:6
#' y2 <- 0.70 / (1 + exp(-0.8 * (t2 - 2)))
#'
#' fit_log <- fit_recovery_model_core(
#'   t = t2,
#'   y = y2,
#'   model = "logistic"
#' )
#'
#' fit_log$coefficients
#'
#' y3 <- 0.70 * exp(-exp(-0.8 * (t2 - 2)))
#'
#' fit_gom <- fit_recovery_model_core(
#'   t = t2,
#'   y = y3,
#'   model = "gompertz"
#' )
#'
#' fit_gom$coefficients
#' }
#'
#' @importFrom stats nls nls.control coef vcov median deviance df.residual
#' @export
fit_recovery_model_core <- function(
  t,
  y,
  model = c("exponential", "logistic", "gompertz")) {

  model <- match.arg(model)

  if (!is.numeric(t) || !is.numeric(y))
    stop("'t' and 'y' must be numeric")

  if (length(t) != length(y))
    stop("'t' and 'y' must have the same length")

  if (length(t) < 4)
    return(NULL)

  df <- data.frame(
    t = as.numeric(t),
    y = as.numeric(y)
  )

  df <- df[
    is.finite(df$t) & is.finite(df$y),
    ,
    drop = FALSE
  ]

  if (nrow(df) < 4)
    return(NULL)

  if (length(unique(df$t)) < 2)
    return(NULL)

  # ---- Model specification ----

  if (model == "exponential") {

    form <- y ~ A * (1 - exp(-k * t))

    start <- list(
      A = min(max(max(df$y), 0.1), 2 - 1e-6),
      k = 0.3
    )

    lower <- c(
      A = 0,
      k = 0.001
    )

    upper <- c(
      A = 2,
      k = 5
    )

  } else if (model == "logistic") {

    form <- y ~ A / (1 + exp(-k * (t - x0)))

    start <- list(
      A = min(max(max(df$y), 0.1), 2 - 1e-6),
      k = 0.3,
      x0 = stats::median(df$t)
    )

    lower <- c(
      A = 0,
      k = 0.001,
      x0 = 0
    )

    upper <- c(
      A = 2,
      k = 5,
      x0 = 30
    )

  } else {

    form <- y ~ A * exp(-exp(-k * (t - x0)))

    start <- list(
      A = min(max(max(df$y), 0.1), 2 - 1e-6),
      k = 0.3,
      x0 = stats::median(df$t)
    )

    lower <- c(
      A = 0,
      k = 0.001,
      x0 = 0
    )

    upper <- c(
      A = 2,
      k = 5,
      x0 = 30
    )
  }

  # ---- Bounded nonlinear least-squares fitting ----

  fit <- tryCatch(
    suppressWarnings(
      stats::nls(
        form,
        data = df,
        start = start,
        algorithm = "port",
        lower = lower,
        upper = upper,
        control = stats::nls.control(
          warnOnly = TRUE,
          maxiter = 200,
          minFactor = 1e-4
        )
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fit))
    return(NULL)

  # ---- Extract coefficients ----

  coefs <- tryCatch(
    stats::coef(fit),
    error = function(e) NULL
  )

  if (is.null(coefs))
    return(NULL)

  if (any(!is.finite(coefs)))
    return(NULL)

  # ---- Sanity checks on fitted parameters ----

  if ("k" %in% names(coefs) && coefs["k"] <= 0)
    return(NULL)

  if ("A" %in% names(coefs) &&
      (coefs["A"] < 0 || coefs["A"] > 2))
    return(NULL)

  if ("x0" %in% names(coefs) &&
      (coefs["x0"] < 0 || coefs["x0"] > 30))
    return(NULL)

  # ---- Variance-covariance matrix ----

  vc <- tryCatch(
    stats::vcov(fit),
    error = function(e) NULL
  )

  # ---- Residual standard error ----

  sigma_val <- tryCatch(
    {
      dev <- stats::deviance(fit)
      dfr <- stats::df.residual(fit)

      if (is.finite(dev) &&
          is.finite(dfr) &&
          dfr > 0) {
        sqrt(dev / dfr)
      } else {
        NA_real_
      }
    },
    error = function(e) NA_real_
  )

  # ---- Return fitted model ----

  list(
    model = model,
    coefficients = coefs,
    sigma = sigma_val,
    vcov = vc,
    data = df,
    fit = fit
  )
}
