#' @title r3d_bwselect: Bandwidth Selection for Distributional RDD
#'
#' @description
#' Computes and returns the bandwidth(s) used by \code{\link{r3d}}. 
#' Implements a three-step pilot procedure to find MSE- or IMSE-optimal bandwidths 
#' under the \code{"simple"} or \code{"frechet"} method.
#'
#' @param X Numeric vector of the running variable.
#' @param Y_list A list of numeric vectors, each representing the samples 
#'   from one unit's outcome distribution.
#' @param q_grid A numeric vector of quantiles at which to compute local fits.
#' @param method Either \code{"simple"} for per-quantile MSE bandwidths 
#'   or \code{"frechet"} for a single IMSE bandwidth.
#' @param s An integer. The expansion order in the pilot local polynomial step (usually 1).
#' @param p The final polynomial order to be used in estimation (usually 2).
#' @param kernel Kernel function used for local weighting; default is triangular.
#' @param cutoff Numeric scalar giving the threshold. 
#'   Internally we recenter so that \code{X - cutoff} is used, matching the \code{X=0} theory.
#' @param ... Additional arguments (unused).
#'
#' @details
#'
#' @return 
#' A list with elements:
#' \item{method}{The method used, \code{"simple"} or \code{"frechet"}.}
#' \item{q_grid}{The quantile grid provided.}
#' \item{h_star}{Either a numeric vector of bandwidths (if \code{method="simple"}), 
#'   or a single scalar (if \code{method="frechet"}).}
#' \item{pilot_h}{The pilot bandwidth(s) used.}
#' \item{s}{The pilot polynomial order.}
#' \item{p}{The final polynomial order.}
#' \item{B_plus, B_minus}{Estimates of leading bias terms for each quantile.}
#' \item{V_plus, V_minus}{Estimates of residual variance terms for each quantile.}
#' \item{f_X_hat}{Estimated density of \code{X} at the cutoff.}
#'
#' @export
r3d_bwselect <- function(X, Y_list,
                         q_grid = seq(0.1, 0.9, 0.1),
                         method = c("simple", "frechet"),
                         s = 1, 
                         p = 2,
                         kernel = function(u) pmax(0, 1 - abs(u)),
                         cutoff = 0,
                         ...)
{
  method <- match.arg(method)
  n <- length(X)
  nQ <- length(q_grid)
  
  ## 1) Re-center the running variable so cutoff becomes 0
  Xc <- X - cutoff
  
  # check that there is data on both sides of the cutoff
  if (all(Xc <= 0) || all(Xc >= 0)) {
    stop("Invalid X data: no data on both sides of the cutoff.")
  }
  
  # -----------------------------------------------------
  # Step 1: Preliminary Bandwidth Calculation
  # -----------------------------------------------------
  
  # Estimate density of X at the cutoff c via Silverman's rule
  # We'll do it at X= cutoff => Xc = 0
  sigma_X <- stats::sd(Xc)
  c_n <- 1.06 * sigma_X * n^(-1/5)
  
  # Evaluate kernel at Xc / c_n
  # so that f_X_hat approximates f_X(cutoff)
  f_X_hat <- mean(kernel(Xc / c_n)) / c_n
  
  # Pilot bandwidth for expansion order s (rule-of-thumb)
  pilot_h <- 1.06 * sigma_X * n^(-1/(2*(s+1)+1))
  
  # Compute empirical quantile matrix
  Qmat <- .compute_empirical_qmat(Y_list, q_grid)
  
  # -----------------------------------------------------
  # Step 2: First-Stage Local Polynomial Fits (order s+1)
  # -----------------------------------------------------
  
  # Kernel weights at pilot_h
  w_plus <- kernel(Xc / pilot_h)
  w_plus[Xc < 0] <- 0
  w_minus <- kernel(Xc / pilot_h)
  w_minus[Xc >= 0] <- 0
  
  
  # Arrays to store bias & variance estimates for each quantile
  B_plus <- numeric(nQ)
  B_minus <- numeric(nQ)
  V_plus <- numeric(nQ)
  V_minus <- numeric(nQ)
  
  for(qi in seq_len(nQ)) {
    # + side
    idx_plus <- (Xc >= 0) & (w_plus > 0)
    if(sum(idx_plus) > (s+1)) {
      Xp <- Xc[idx_plus]
      Yp <- Qmat[idx_plus, qi]
      Wp <- w_plus[idx_plus]
      X_mat_p <- outer(Xp / pilot_h, 0:s, "^")
      
      W_diag_p <- diag(Wp)
      XWX_p <- t(X_mat_p) %*% W_diag_p %*% X_mat_p
      XWY_p <- t(X_mat_p) %*% W_diag_p %*% Yp
      
      if(rcond(XWX_p) > 1e-10) {
        coef_p <- solve(XWX_p, XWY_p)
        # Leading bias = coefficient on x^(s+1)
        B_plus[qi] <- coef_p[s+1]
        
        # Weighted residual variance
        fitted_p <- X_mat_p %*% coef_p
        resid_p  <- Yp - fitted_p
        V_plus[qi] <- sum(Wp * resid_p^2) / sum(Wp)
      }
    }
    
    # - side
    idx_minus <- (Xc < 0) & (w_minus > 0)
    if(sum(idx_minus) > (s+1)) {
      Xm <- Xc[idx_minus]
      Ym <- Qmat[idx_minus, qi]
      Wm <- w_minus[idx_minus]
      X_mat_m <- outer(Xm / pilot_h, 0:s, "^")
      
      W_diag_m <- diag(Wm)
      XWX_m <- t(X_mat_m) %*% W_diag_m %*% X_mat_m
      XWY_m <- t(X_mat_m) %*% W_diag_m %*% Ym
      
      if(rcond(XWX_m) > 1e-10) {
        coef_m <- solve(XWX_m, XWY_m)
        B_minus[qi] <- coef_m[s+1]
        
        fitted_m <- X_mat_m %*% coef_m
        resid_m  <- Ym - fitted_m
        V_minus[qi] <- sum(Wm * resid_m^2) / sum(Wm)
      }
    }
  }
  
  # -----------------------------------------------------
  # Step 3: Calculate Optimal Bandwidth
  # -----------------------------------------------------
  
  # Leading factor for bias
  bias_constant <- 1/factorial(s+1)
  
  if(method == "simple") {
    # Per-quantile MSE bandwidth
    h_star <- numeric(nQ)
    
    for(qi in seq_len(nQ)) {
      # Bias^2
      bias_diff <- (B_plus[qi] - B_minus[qi]) * bias_constant
      b2 <- bias_diff^2
      
      # Variance part must be scaled by 1/f_X_hat
      var_sum <- (V_plus[qi] + V_minus[qi]) / f_X_hat
      
      # If near zero => skip
      if(abs(b2) < 1e-14) {
        h_star[qi] <- NA
        next
      }
      
      # MSE formula: 
      #  h* ~ [ var_sum / (2*(s+1)* bias^2 ) ]^(1/(2s+3)) * n^(-1/(2s+3))
      ratio <- (1/(2*(s+1))) * (var_sum / b2)
      
      h_star[qi] <- ratio^(1/(2*s+3)) * n^(-1/(2*s+3))
    }
    
    # If any are NA => fallback
    if(any(is.na(h_star) | h_star <= 0)) {
      warning("Some MSE-optimal bandwidths invalid; replacing with median of valid ones.")
      valid_h <- h_star[!is.na(h_star) & h_star>0]
      if(length(valid_h)==0) {
        h_star <- 1.06 * sigma_X * n^(-1/(2*s+1))
      } else {
        h_star[is.na(h_star) | h_star<=0] <- median(valid_h)
      }
    }
    
  } else {
    # Single IMSE bandwidth
    # integrate bias^2 and variance across q_grid
    dq <- diff(range(q_grid)) / (length(q_grid)-1)
    
    bias_diffs <- (B_plus - B_minus)*bias_constant
    A_s <- sum(bias_diffs^2) * dq
    
    # variance scaled by 1/f_X_hat
    B_s <- sum(V_plus + V_minus) * dq / f_X_hat
    
    # h_star = [ B_s / (2*(p+1)*A_s ) ]^(1/(2p+3)) * n^(-1/(2p+3))
    denom <- max(A_s, 1e-14) # avoid zero
    ratio <- B_s / (2*(s+1)* denom)
    h_ <- ratio^(1/(2*s+3)) * n^(-1/(2*s+3))
    
    if(is.na(h_) || h_ <= 0) {
      warning("IMSE-optimal bandwidth invalid; using rule-of-thumb.")
      h_ <- 1.06*sigma_X* n^(-1/(2*s+1))
    }
    h_star <- h_
  }
  
  # Return
  list(
    method = method,
    q_grid = q_grid,
    h_star = h_star,
    pilot_h = pilot_h,
    s = s,
    p = p,
    B_plus = B_plus,
    B_minus = B_minus,
    V_plus = V_plus,
    V_minus = V_minus,
    f_X_hat = f_X_hat
  )
}