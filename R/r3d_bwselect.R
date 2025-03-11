#' Bandwidth Selection for Distributional RDD
#'
#' Computes bandwidth(s) for local polynomial estimation in a regression discontinuity
#' setting with distributional (Fréchet-valued) outcomes. Implements a three-step pilot
#' procedure to find either MSE-optimal (per-quantile) bandwidths or a single IMSE-optimal
#' bandwidth, depending on \code{method}.
#'
#' @param X Numeric vector of the running variable.
#' @param Y_list A list of numeric vectors; each entry is the sample of outcomes
#'   from one unit's distribution.
#' @param T (Optional) Numeric or logical vector of treatment statuses for fuzzy design.
#' @param q_grid Numeric vector of quantiles at which local polynomial fits are performed.
#' @param method Either \code{"simple"} (per-quantile MSE-optimal) or
#'   \code{"frechet"} (single IMSE-optimal bandwidth).
#' @param s Integer specifying the order of local polynomial in the pilot stage (often 1).
#' @param p Integer specifying the final local polynomial order (often 2).
#' @param kernel Kernel function for local weighting. Defaults to triangular kernel.
#' @param cutoff Numeric scalar threshold. Data are recentered so \code{X - cutoff} has cutoff at 0.
#' @param fuzzy Logical indicating fuzzy design. Default is FALSE.
#' @param ... Additional arguments for future expansions.
#'
#' @details
#' Implements a three-step procedure:
#' \enumerate{
#'   \item Estimates \( f_X(0) \) using Silverman’s rule and computes pilot bandwidths via global polynomials.
#'   \item Runs pilot local polynomial regressions to estimate bias and variance.
#'   \item Computes MSE-optimal (per-quantile) or IMSE-optimal (single) bandwidths.
#' }
#' In fuzzy RDD, separate bandwidths are computed for the numerator (outcome) and denominator (treatment).
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{method}}{Method used: \code{"simple"} or \code{"frechet"}.}
#'   \item{\code{q_grid}}{Input \code{q_grid}.}
#'   \item{\code{h_star_num}}{Bandwidth(s) for numerator (outcome).}
#'   \item{\code{h_star_den}}{Bandwidth for denominator (treatment, if fuzzy).}
#'   \item{\code{pilot_h_num}}{Pilot bandwidth(s) for numerator.}
#'   \item{\code{pilot_h_den}}{Pilot bandwidth for denominator (if fuzzy).}
#'   \item{\code{s}, \code{p}}{Polynomial orders.}
#'   \item{\code{B_plus}, \code{B_minus}}{Bias estimates for numerator.}
#'   \item{\code{V_plus}, \code{V_minus}}{Variance estimates for numerator.}
#'   \item{\code{f_X_hat}}{Estimated density of \( X \) at cutoff.}
#' }
#'
#' @references
#' Van Dijcke, D. (2025). \emph{Regression Discontinuity Design with Distributional Outcomes (R3D)}.
#' Working paper.
#'
#' @export
r3d_bwselect <- function(X, Y_list, T = NULL,
                         q_grid = seq(0.1, 0.9, 0.1),
                         method = c("simple", "frechet"),
                         s = 1, 
                         p = 2,
                         kernel = function(u) pmax(0, 1 - abs(u)),
                         cutoff = 0,
                         fuzzy = FALSE,
                         ...)
{
  method <- match.arg(method)
  n <- length(X)
  nQ <- length(q_grid)
  
  if (fuzzy && is.null(T)) stop("Fuzzy design requires treatment variable T.")
  
  # Re-center running variable
  Xc <- X - cutoff
  if (all(Xc <= 0) || all(Xc >= 0)) stop("Invalid X data: no data on both sides of cutoff.")
  
  # Step 1: Preliminary Bandwidth Calculation
  # Estimate f_X(0) using Silverman’s rule
  sigma_X <- stats::sd(Xc)
  c_n <- 1.06 * sigma_X * n^(-1/5)
  f_X_hat <- mean(kernel(Xc / c_n)) / c_n
  
  # Compute empirical quantile matrix for Y
  Qmat <- .compute_empirical_qmat(Y_list, q_grid)
  
  # Global polynomial fits for pilot bandwidths
  fit_global_poly <- function(Y, X, order) {
    fit <- lm(Y ~ poly(X, degree = order, raw = TRUE))
    coefs <- coef(fit)
    derivs <- coefs[2:(order + 1)] * factorial(0:(order - 1))
    resid_var <- var(resid(fit))
    list(derivs = derivs, resid_var = resid_var)
  }
  
  # Numerator (outcome)
  B_plus_num <- numeric(nQ)
  B_minus_num <- numeric(nQ)
  V_plus_num <- numeric(nQ)
  V_minus_num <- numeric(nQ)
  
  for (qi in seq_len(nQ)) {
    idx_plus <- Xc >= 0
    idx_minus <- Xc < 0
    if (sum(idx_plus) > (s + 1)) {
      fit_plus <- fit_global_poly(Qmat[idx_plus, qi], Xc[idx_plus], s + 1)
      B_plus_num[qi] <- fit_plus$derivs[s + 1] / factorial(s + 1)
      V_plus_num[qi] <- fit_plus$resid_var
    }
    if (sum(idx_minus) > (s + 1)) {
      fit_minus <- fit_global_poly(Qmat[idx_minus, qi], Xc[idx_minus], s + 1)
      B_minus_num[qi] <- fit_minus$derivs[s + 1] / factorial(s + 1)
      V_minus_num[qi] <- fit_minus$resid_var
    }
  }
  
  # Denominator (treatment) in fuzzy case
  if (fuzzy) {
    fit_T_plus <- fit_global_poly(T[Xc >= 0], Xc[Xc >= 0], s + 1)
    fit_T_minus <- fit_global_poly(T[Xc < 0], Xc[Xc < 0], s + 1)
    B_plus_den <- fit_T_plus$derivs[s + 1] / factorial(s + 1)
    B_minus_den <- fit_T_minus$derivs[s + 1] / factorial(s + 1)
    V_plus_den <- fit_T_plus$resid_var
    V_minus_den <- fit_T_minus$resid_var
  }
  
  # Compute pilot bandwidths
  if (method == "simple") {
    pilot_h_num <- numeric(nQ)
    for (qi in seq_len(nQ)) {
      bias_diff <- B_plus_num[qi] - B_minus_num[qi]
      var_sum <- V_plus_num[qi] + V_minus_num[qi]
      if (abs(bias_diff) > 1e-14) {
        ratio <- (var_sum) / (2 * (s + 1) * bias_diff^2)
        pilot_h_num[qi] <- ratio^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
      } else {
        pilot_h_num[qi] <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
      }
    }
  } else {
    dq <- diff(range(q_grid)) / (length(q_grid) - 1)
    A_s_num <- sum((B_plus_num - B_minus_num)^2) * dq
    B_s_num <- sum(V_plus_num + V_minus_num) * dq
    ratio <- B_s_num / (2 * (s + 1) * max(A_s_num, 1e-14))
    pilot_h_num <- ratio^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
    if (is.na(pilot_h_num) || pilot_h_num <= 0) pilot_h_num <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
  }
  
  if (fuzzy) {
    bias_diff_den <- B_plus_den - B_minus_den
    var_sum_den <- V_plus_den + V_minus_den
    if (abs(bias_diff_den) > 1e-14) {
      ratio_den <- var_sum_den / (2 * (s + 1) * bias_diff_den^2)
      pilot_h_den <- ratio_den^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
    } else {
      pilot_h_den <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
    }
  }
  
  # Step 2: First-Stage Local Polynomial Fits
  B_plus <- numeric(nQ)
  B_minus <- numeric(nQ)
  V_plus <- numeric(nQ)
  V_minus <- numeric(nQ)
  
  for (qi in seq_len(nQ)) {
    h_pilot <- if (method == "simple") pilot_h_num[qi] else pilot_h_num
    idx_plus <- (Xc >= 0) & (abs(Xc) <= h_pilot)
    if (sum(idx_plus) > (s + 1)) {
      Xp <- Xc[idx_plus]
      Yp <- Qmat[idx_plus, qi]
      Wp <- kernel(Xp / h_pilot)
      X_mat_p <- outer(Xp / h_pilot, 0:s, "^")
      W_diag_p <- diag(Wp)
      XWX_p <- t(X_mat_p) %*% W_diag_p %*% X_mat_p
      XWY_p <- t(X_mat_p) %*% W_diag_p %*% Yp
      if (rcond(XWX_p) > 1e-10) {
        coef_p <- solve(XWX_p, XWY_p)
        B_plus[qi] <- coef_p[s + 1] / factorial(s + 1)
        fitted_p <- X_mat_p %*% coef_p
        V_plus[qi] <- sum(Wp * (Yp - fitted_p)^2) / sum(Wp)
      }
    }
    idx_minus <- (Xc < 0) & (abs(Xc) <= h_pilot)
    if (sum(idx_minus) > (s + 1)) {
      Xm <- Xc[idx_minus]
      Ym <- Qmat[idx_minus, qi]
      Wm <- kernel(Xm / h_pilot)
      X_mat_m <- outer(Xm / h_pilot, 0:s, "^")
      W_diag_m <- diag(Wm)
      XWX_m <- t(X_mat_m) %*% W_diag_m %*% X_mat_m
      XWY_m <- t(X_mat_m) %*% W_diag_m %*% Ym
      if (rcond(XWX_m) > 1e-10) {
        coef_m <- solve(XWX_m, XWY_m)
        B_minus[qi] <- coef_m[s + 1] / factorial(s + 1)
        fitted_m <- X_mat_m %*% coef_m
        V_minus[qi] <- sum(Wm * (Ym - fitted_m)^2) / sum(Wm)
      }
    }
  }
  
  # Denominator in fuzzy case
  if (fuzzy) {
    idx_plus <- (Xc >= 0) & (abs(Xc) <= pilot_h_den)
    if (sum(idx_plus) > (s + 1)) {
      Xp <- Xc[idx_plus]
      Tp <- T[idx_plus]
      Wp <- kernel(Xp / pilot_h_den)
      X_mat_p <- outer(Xp / pilot_h_den, 0:s, "^")
      W_diag_p <- diag(Wp)
      XWX_p <- t(X_mat_p) %*% W_diag_p %*% X_mat_p
      XWY_p <- t(X_mat_p) %*% W_diag_p %*% Tp
      if (rcond(XWX_p) > 1e-10) {
        coef_p <- solve(XWX_p, XWY_p)
        B_plus_den <- coef_p[s + 1] / factorial(s + 1)
        fitted_p <- X_mat_p %*% coef_p
        V_plus_den <- sum(Wp * (Tp - fitted_p)^2) / sum(Wp)
      }
    }
    idx_minus <- (Xc < 0) & (abs(Xc) <= pilot_h_den)
    if (sum(idx_minus) > (s + 1)) {
      Xm <- Xc[idx_minus]
      Tm <- T[idx_minus]
      Wm <- kernel(Xm / pilot_h_den)
      X_mat_m <- outer(Xm / pilot_h_den, 0:s, "^")
      W_diag_m <- diag(Wm)
      XWX_m <- t(X_mat_m) %*% W_diag_m %*% X_mat_m
      XWY_m <- t(X_mat_m) %*% W_diag_m %*% Tm
      if (rcond(XWX_m) > 1e-10) {
        coef_m <- solve(XWX_m, XWY_m)
        B_minus_den <- coef_m[s + 1] / factorial(s + 1)
        fitted_m <- X_mat_m %*% coef_m
        V_minus_den <- sum(Wm * (Tm - fitted_m)^2) / sum(Wm)
      }
    }
  }
  
  # Step 3: Optimal Bandwidth
  if (method == "simple") {
    h_star_num <- numeric(nQ)
    for (qi in seq_len(nQ)) {
      bias_diff <- B_plus[qi] - B_minus[qi]
      var_sum <- (V_plus[qi] + V_minus[qi]) / f_X_hat
      if (abs(bias_diff) > 1e-14) {
        ratio <- (var_sum) / (2 * (s + 1) * bias_diff^2)
        h_star_num[qi] <- ratio^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
      } else {
        h_star_num[qi] <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
      }
    }
  } else {
    dq <- diff(range(q_grid)) / (length(q_grid) - 1)
    A_s <- sum((B_plus - B_minus)^2) * dq
    B_s <- sum(V_plus + V_minus) * dq / f_X_hat
    ratio <- B_s / (2 * (s + 1) * max(A_s, 1e-14))
    h_star_num <- ratio^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
    if (is.na(h_star_num) || h_star_num <= 0) h_star_num <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
  }
  
  if (fuzzy) {
    bias_diff_den <- B_plus_den - B_minus_den
    var_sum_den <- (V_plus_den + V_minus_den) / f_X_hat
    if (abs(bias_diff_den) > 1e-14) {
      ratio_den <- var_sum_den / (2 * (s + 1) * bias_diff_den^2)
      h_star_den <- ratio_den^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
    } else {
      h_star_den <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
    }
  } else {
    h_star_den <- NULL
  }
  
  # Return results
  list(
    method = method,
    q_grid = q_grid,
    h_star_num = h_star_num,
    h_star_den = h_star_den,
    pilot_h_num = pilot_h_num,
    pilot_h_den = if (fuzzy) pilot_h_den else NULL,
    s = s,
    p = p,
    B_plus = B_plus,
    B_minus = B_minus,
    V_plus = V_plus,
    V_minus = V_minus,
    f_X_hat = f_X_hat
  )
}