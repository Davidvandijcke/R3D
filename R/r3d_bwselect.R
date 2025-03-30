
#' Bandwidth Selection for R3D
#'
#' Computes bandwidth(s) for local polynomial estimation in a regression discontinuity
#' setting with distribution-valued outcomes. Implements a three-step pilot
#' procedure to find either MSE-optimal (per-quantile) bandwidths or a single IMSE-optimal
#' bandwidth, depending on \code{method}. For details, see the Appendix of \insertCite{vandijcke2025;textual}{R3D}.
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
#' @param coverage Logical indicating whether to apply the coverage correction rule of thumb of
#' \insertCite{calonico2018effect;textual}{R3D}. Default is FALSE.
#' @param ... Additional arguments for future expansions.
#'
#' @details
#' Implements a three-step procedure:
#' \enumerate{
#'   \item Estimates \eqn{f_X(0)} using Silverman’s rule and computes pilot bandwidths via global polynomials.
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
#'   \item{\code{f_X_hat}}{Estimated density of \eqn{X} at cutoff.}
#' }
#'
#' @references
#' \insertAllCited{}
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
                         coverage = FALSE,
                         ...)
{
  method <- match.arg(method)
  n <- length(X)
  nQ <- length(q_grid)
  
  if (fuzzy && is.null(T)) stop("Fuzzy design requires treatment variable T.")
  
  # Re-center running variable
  Xc <- X - cutoff
  if (all(Xc <= 0) || all(Xc >= 0)) stop("Invalid X data: no data on both sides of cutoff.")
  
  # Define kernel type as integer for Fortran call
  kernel_fun <- deparse(substitute(kernel))
  kernel_type <- as.integer(switch(kernel_fun,
                                   "function(u) pmax(0, 1 - abs(u))" = 1, # triangular
                                   "function(u) 0.75 * pmax(0, 1 - u^2)" = 2, # epanechnikov
                                   "function(u) 0.5 * (abs(u) <= 1)" = 3, # uniform
                                   1)) # default to triangular
  
  # ------------------------------------------------
  # Step 1: Preliminary Bandwidth Calculation
  # ------------------------------------------------
  
  # 1(i) Estimate f_X(0) using Silverman's rule
  sigma_X <- stats::sd(Xc)
  c_n <- 1.06 * sigma_X * n^(-1/5)
  f_X_hat <- mean(kernel(Xc / c_n)) / c_n
  
  # Compute empirical quantile matrix for Y
  Qmat <- .compute_empirical_qmat(Y_list, q_grid)
  
  # Calculate kernel matrices needed for bias and variance formulas
  # These are Γ_±,s, Ψ_±,s, Λ_{s,s+1}^± as defined in the paper
  compute_kernel_matrices <- function(s, kernel_fun) {
    # Define integration grid
    grid <- seq(-1, 1, length.out = 1000)
    
    # Create r_s polynomial basis function
    r_s <- function(u, s) {
      outer(u, 0:s, "^")
    }
    
    # Split grid for plus/minus sides
    grid_plus <- grid[grid >= 0]
    grid_minus <- grid[grid < 0]
    
    # Compute kernel values
    k_plus <- sapply(grid_plus, kernel)
    k_minus <- sapply(grid_minus, kernel)
    
    # Compute r_s matrices
    r_plus <- r_s(grid_plus, s)
    r_minus <- r_s(grid_minus, s)
    
    # Compute Gamma matrices (Γ_±,s)
    Gamma_plus <- matrix(0, s+1, s+1)
    Gamma_minus <- matrix(0, s+1, s+1)
    
    for (i in 1:length(grid_plus)) {
      r_i <- r_plus[i,]
      Gamma_plus <- Gamma_plus + k_plus[i] * outer(r_i, r_i) * (2/length(grid))
    }
    
    for (i in 1:length(grid_minus)) {
      r_i <- r_minus[i,]
      Gamma_minus <- Gamma_minus + k_minus[i] * outer(r_i, r_i) * (2/length(grid))
    }
    
    # Compute Lambda vectors (Λ_{s,s+1}^±)
    Lambda_plus <- numeric(s+1)
    Lambda_minus <- numeric(s+1)
    
    for (i in 1:length(grid_plus)) {
      Lambda_plus <- Lambda_plus + grid_plus[i]^(s+1) * k_plus[i] * r_plus[i,] * (2/length(grid))
    }
    
    for (i in 1:length(grid_minus)) {
      Lambda_minus <- Lambda_minus + grid_minus[i]^(s+1) * k_minus[i] * r_minus[i,] * (2/length(grid))
    }
    
    # Compute Psi matrices (Ψ_±,s)
    Psi_plus <- matrix(0, s+1, s+1)
    Psi_minus <- matrix(0, s+1, s+1)
    
    for (i in 1:length(grid_plus)) {
      r_i <- r_plus[i,]
      Psi_plus <- Psi_plus + k_plus[i]^2 * outer(r_i, r_i) * (2/length(grid))
    }
    
    for (i in 1:length(grid_minus)) {
      r_i <- r_minus[i,]
      Psi_minus <- Psi_minus + k_minus[i]^2 * outer(r_i, r_i) * (2/length(grid))
    }
    
    return(list(
      Gamma_plus = Gamma_plus,
      Gamma_minus = Gamma_minus,
      Lambda_plus = Lambda_plus,
      Lambda_minus = Lambda_minus,
      Psi_plus = Psi_plus,
      Psi_minus = Psi_minus
    ))
  }
  
  # Calculate kernel matrices
  kernel_matrices <- compute_kernel_matrices(s, kernel)
  
  # 1(ii) Global polynomial fit approach for preliminary derivative estimates
  fit_global_poly <- function(Y, X, order) {
    fit <- try(lm(Y ~ poly(X, degree = order, raw = TRUE)), silent = TRUE)
    
    if (inherits(fit, "try-error") || is.null(fit)) {
      return(list(derivs = rep(0, order), resid_var = 1))
    }
    
    coefs <- coef(fit)
    derivs <- rep(0, order)
    
    # Calculate derivatives adjusted by factorial terms
    for (i in 1:order) {
      if ((i+1) <= length(coefs) && !is.na(coefs[i+1])) {
        derivs[i] <- coefs[i+1] * factorial(i)
      }
    }
    
    resid_var <- var(resid(fit))
    if (is.na(resid_var)) resid_var <- 1
    
    return(list(derivs = derivs, resid_var = resid_var))
  }
  
  # Initial global polynomial fits for derivatives and variances
  pilot_derivs_plus <- numeric(nQ)
  pilot_derivs_minus <- numeric(nQ)
  pilot_vars_plus <- numeric(nQ)
  pilot_vars_minus <- numeric(nQ)
  
  for (qi in seq_len(nQ)) {
    idx_plus <- Xc >= 0
    idx_minus <- Xc < 0
    
    if (sum(idx_plus) > (s + 1)) {
      fit_plus <- fit_global_poly(Qmat[idx_plus, qi], Xc[idx_plus], s + 1)
      pilot_derivs_plus[qi] <- fit_plus$derivs[s + 1]
      pilot_vars_plus[qi] <- fit_plus$resid_var
    }
    
    if (sum(idx_minus) > (s + 1)) {
      fit_minus <- fit_global_poly(Qmat[idx_minus, qi], Xc[idx_minus], s + 1)
      pilot_derivs_minus[qi] <- fit_minus$derivs[s + 1]
      pilot_vars_minus[qi] <- fit_minus$resid_var
    }
  }
  
  # For fuzzy RDD, get treatment derivatives
  if (fuzzy) {
    idx_plus <- Xc >= 0
    idx_minus <- Xc < 0
    
    fit_T_plus <- fit_global_poly(T[idx_plus], Xc[idx_plus], s + 1)
    fit_T_minus <- fit_global_poly(T[idx_minus], Xc[idx_minus], s + 1)
    
    pilot_deriv_T_plus <- fit_T_plus$derivs[s + 1]
    pilot_deriv_T_minus <- fit_T_minus$derivs[s + 1]
    pilot_var_T_plus <- fit_T_plus$resid_var
    pilot_var_T_minus <- fit_T_minus$resid_var
  }
  
  # Compute e_0 vector (1, 0, ..., 0)
  e_0 <- c(1, rep(0, s))
  
  # Compute Gamma inverse matrices
  Gamma_plus_inv <- try(solve(kernel_matrices$Gamma_plus), silent = TRUE)
  Gamma_minus_inv <- try(solve(kernel_matrices$Gamma_minus), silent = TRUE)
  
  if (inherits(Gamma_plus_inv, "try-error") || inherits(Gamma_minus_inv, "try-error")) {
    warning("Gamma matrices singular. Using regularized inverse.")
    # Use regularized inverse if matrix is singular
    Gamma_plus_inv <- solve(kernel_matrices$Gamma_plus + 1e-8 * diag(s+1))
    Gamma_minus_inv <- solve(kernel_matrices$Gamma_minus + 1e-8 * diag(s+1))
  }
  
  # Compute pilot bandwidths following Appendix A.3
  if (method == "simple") {
    pilot_h_num <- numeric(nQ)
    
    for (qi in seq_len(nQ)) {
      # Compute C_1,0(q) as in the paper
      bias_plus <- as.numeric(
        t(e_0) %*% Gamma_plus_inv %*% kernel_matrices$Lambda_plus * 
          (pilot_derivs_plus[qi] / factorial(s + 1))
      )
      
      bias_minus <- as.numeric(
        t(e_0) %*% Gamma_minus_inv %*% kernel_matrices$Lambda_minus * 
          (pilot_derivs_minus[qi] / factorial(s + 1))
      )
      
      C_1_0 <- bias_plus - bias_minus
      
      # Compute C_1,0'(q) as in the paper
      var_plus <- as.numeric(
        pilot_vars_plus[qi] * t(e_0) %*% Gamma_plus_inv %*% 
          kernel_matrices$Psi_plus %*% Gamma_plus_inv %*% e_0
      )
      
      var_minus <- as.numeric(
        pilot_vars_minus[qi] * t(e_0) %*% Gamma_minus_inv %*% 
          kernel_matrices$Psi_minus %*% Gamma_minus_inv %*% e_0
      )
      
      C_1_0_prime <- (var_plus + var_minus) / f_X_hat
      
      # Compute pilot bandwidth
      if (abs(C_1_0) > 1e-14) {
        ratio <- C_1_0_prime / (2 * (s + 1) * C_1_0^2)
        pilot_h_num[qi] <- ratio^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
      } else {
        # Fallback to Silverman's rule
        pilot_h_num[qi] <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
      }
      
      # Ensure positive bandwidth
      if (is.na(pilot_h_num[qi]) || pilot_h_num[qi] <= 0) {
        pilot_h_num[qi] <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
      }
    }
  } else {
    # For Frechet method, compute IMSE-optimal bandwidth
    dq <- diff(range(q_grid)) / (length(q_grid) - 1)
    
    # Compute bias terms A_s
    A_s_terms <- numeric(nQ)
    for (qi in seq_len(nQ)) {
      bias_plus <- as.numeric(
        t(e_0) %*% Gamma_plus_inv %*% kernel_matrices$Lambda_plus * 
          (pilot_derivs_plus[qi] / factorial(s + 1))
      )
      
      bias_minus <- as.numeric(
        t(e_0) %*% Gamma_minus_inv %*% kernel_matrices$Lambda_minus * 
          (pilot_derivs_minus[qi] / factorial(s + 1))
      )
      
      A_s_terms[qi] <- (bias_plus - bias_minus)^2
    }
    A_s <- sum(A_s_terms) * dq
    
    # Compute variance terms B_s
    B_s_terms <- numeric(nQ)
    for (qi in seq_len(nQ)) {
      var_plus <- as.numeric(
        pilot_vars_plus[qi] * t(e_0) %*% Gamma_plus_inv %*% 
          kernel_matrices$Psi_plus %*% Gamma_plus_inv %*% e_0
      )
      
      var_minus <- as.numeric(
        pilot_vars_minus[qi] * t(e_0) %*% Gamma_minus_inv %*% 
          kernel_matrices$Psi_minus %*% Gamma_minus_inv %*% e_0
      )
      
      B_s_terms[qi] <- (var_plus + var_minus) / f_X_hat
    }
    B_s <- sum(B_s_terms) * dq
    
    # Compute IMSE-optimal pilot bandwidth
    if (A_s > 1e-14) {
      ratio <- B_s / (2 * (s + 1) * A_s)
      pilot_h_num <- ratio^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
    } else {
      pilot_h_num <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
    }
    
    # Ensure positive bandwidth
    if (is.na(pilot_h_num) || pilot_h_num <= 0) {
      pilot_h_num <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
    }
  }
  
  # Compute pilot bandwidth for treatment (fuzzy case)
  if (fuzzy) {
    # Compute bias term for treatment
    bias_T_plus <- as.numeric(
      t(e_0) %*% Gamma_plus_inv %*% kernel_matrices$Lambda_plus * 
        (pilot_deriv_T_plus / factorial(s + 1))
    )
    
    bias_T_minus <- as.numeric(
      t(e_0) %*% Gamma_minus_inv %*% kernel_matrices$Lambda_minus * 
        (pilot_deriv_T_minus / factorial(s + 1))
    )
    
    C_T_0 <- bias_T_plus - bias_T_minus
    
    # Compute variance term for treatment
    var_T_plus <- as.numeric(
      pilot_var_T_plus * t(e_0) %*% Gamma_plus_inv %*% 
        kernel_matrices$Psi_plus %*% Gamma_plus_inv %*% e_0
    )
    
    var_T_minus <- as.numeric(
      pilot_var_T_minus * t(e_0) %*% Gamma_minus_inv %*% 
        kernel_matrices$Psi_minus %*% Gamma_minus_inv %*% e_0
    )
    
    C_T_0_prime <- (var_T_plus + var_T_minus) / f_X_hat
    
    # Compute pilot bandwidth for treatment
    if (abs(C_T_0) > 1e-14) {
      ratio_T <- C_T_0_prime / (2 * (s + 1) * C_T_0^2)
      pilot_h_den <- ratio_T^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
    } else {
      pilot_h_den <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
    }
    
    # Ensure positive bandwidth
    if (is.na(pilot_h_den) || pilot_h_den <= 0) {
      pilot_h_den <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
    }
  }
  
  # ------------------------------------------------
  # Step 2: First-Stage Local Polynomial Fits
  # ------------------------------------------------
  
  # Use Fortran's locweights for efficiency
  # Prepare for Fortran calls
  N <- as.integer(n)
  S_ <- as.integer(s)
  NQ_ <- as.integer(nQ)
  
  # Prepare pilot bandwidths for Fortran call
  if (method == "simple") {
    h_pilot_vec <- pilot_h_num
  } else {
    h_pilot_vec <- rep(pilot_h_num, nQ)
  }
  
  # Call Fortran for first-stage local polynomial fits - plus side
  outPlus <- .Fortran("locweights",
                      X = as.double(Xc),
                      YMAT = as.double(Qmat),
                      N = N,
                      P = S_,
                      H = as.double(h_pilot_vec),
                      SIDE = as.integer(1),  # Plus side
                      KERNEL_TYPE = kernel_type,
                      ALPHA = double((s + 1) * nQ),
                      WINT = double(n * nQ),
                      INFO = integer(1),
                      NQ = NQ_,
                      PACKAGE = "R3D")
  
  # Call Fortran for first-stage local polynomial fits - minus side
  outMinus <- .Fortran("locweights",
                       X = as.double(Xc),
                       YMAT = as.double(Qmat),
                       N = N,
                       P = S_,
                       H = as.double(h_pilot_vec),
                       SIDE = as.integer(0),  # Minus side
                       KERNEL_TYPE = kernel_type,
                       ALPHA = double((s + 1) * nQ),
                       WINT = double(n * nQ),
                       INFO = integer(1),
                       NQ = NQ_,
                       PACKAGE = "R3D")
  
  # Check for issues with Fortran call
  info_plus <- outPlus$INFO
  info_minus <- outMinus$INFO
  if (info_plus != 0) warning("locweights: plus side singular system? info=", info_plus)
  if (info_minus != 0) warning("locweights: minus side singular system? info=", info_minus)
  
  # Extract coefficients and weights
  alpha_plus <- matrix(outPlus$ALPHA, nrow = s + 1, ncol = nQ)
  alpha_minus <- matrix(outMinus$ALPHA, nrow = s + 1, ncol = nQ)
  w_plus <- matrix(outPlus$WINT, nrow = n, ncol = nQ)
  w_minus <- matrix(outMinus$WINT, nrow = n, ncol = nQ)
  
  # For fuzzy design, we need to do the same for T
  if (fuzzy) {
    # Create a matrix for T (just one column)
    T_mat <- matrix(T, nrow = n, ncol = 1)
    
    # Call Fortran for first-stage local polynomial fits - plus side for T
    outTplus <- .Fortran("locweights",
                         X = as.double(Xc),
                         YMAT = as.double(T_mat),
                         N = N,
                         P = S_,
                         H = as.double(pilot_h_den),
                         SIDE = as.integer(1),  # Plus side
                         KERNEL_TYPE = kernel_type,
                         ALPHA = double(s + 1),
                         WINT = double(n),
                         INFO = integer(1),
                         NQ = as.integer(1),
                         PACKAGE = "R3D")
    
    # Call Fortran for first-stage local polynomial fits - minus side for T
    outTminus <- .Fortran("locweights",
                          X = as.double(Xc),
                          YMAT = as.double(T_mat),
                          N = N,
                          P = S_,
                          H = as.double(pilot_h_den),
                          SIDE = as.integer(0),  # Minus side
                          KERNEL_TYPE = kernel_type,
                          ALPHA = double(s + 1),
                          WINT = double(n),
                          INFO = integer(1),
                          NQ = as.integer(1),
                          PACKAGE = "R3D")
    
    # Extract coefficients and weights for T
    alphaT_plus <- matrix(outTplus$ALPHA, nrow = s + 1, ncol = 1)
    alphaT_minus <- matrix(outTminus$ALPHA, nrow = s + 1, ncol = 1)
    wT_plus <- outTplus$WINT
    wT_minus <- outTminus$WINT
  }
  
  # ------------------------------------------------
  # Calculate bias and variance from first-stage fits
  # ------------------------------------------------
  
  # For quantile outcomes
  B_plus <- numeric(nQ)
  B_minus <- numeric(nQ)
  V_plus <- numeric(nQ)
  V_minus <- numeric(nQ)
  
  for (qi in seq_len(nQ)) {
    # Extract the bias term (s+1 derivative coefficient)
    # Properly scaled according to formula in paper
    B_plus[qi] <- alpha_plus[s + 1, qi] / factorial(s+1)
    B_minus[qi] <- alpha_minus[s + 1, qi] / factorial(s+1)
    
    # Calculate variance using weights and residuals
    # Plus side
    idx_plus <- which(w_plus[, qi] > 0)
    if (length(idx_plus) > 0) {
      # Compute fitted values
      h_qi <- if (method == "simple") h_pilot_vec[qi] else h_pilot_vec[1]
      Xp <- Xc[idx_plus] / h_qi
      Xpow_p <- outer(Xp, 0:s, "^")
      fitted_p <- Xpow_p %*% alpha_plus[, qi]
      
      # Compute variance
      residuals_p <- Qmat[idx_plus, qi] - fitted_p
      V_plus[qi] <- sum(w_plus[idx_plus, qi] * residuals_p^2) / sum(w_plus[idx_plus, qi])
    }
    
    # Minus side
    idx_minus <- which(w_minus[, qi] > 0)
    if (length(idx_minus) > 0) {
      # Compute fitted values
      h_qi <- if (method == "simple") h_pilot_vec[qi] else h_pilot_vec[1]
      Xm <- Xc[idx_minus] / h_qi
      Xpow_m <- outer(Xm, 0:s, "^")
      fitted_m <- Xpow_m %*% alpha_minus[, qi]
      
      # Compute variance
      residuals_m <- Qmat[idx_minus, qi] - fitted_m
      V_minus[qi] <- sum(w_minus[idx_minus, qi] * residuals_m^2) / sum(w_minus[idx_minus, qi])
    }
  }
  
  # For treatment variable (fuzzy case)
  if (fuzzy) {
    # Extract bias term
    B_plus_den <- alphaT_plus[s + 1] / factorial(s+1)
    B_minus_den <- alphaT_minus[s + 1] / factorial(s+1)
    
    # Calculate variance using weights and residuals
    # Plus side
    idx_plus <- which(wT_plus > 0)
    if (length(idx_plus) > 0) {
      Xp <- Xc[idx_plus] / pilot_h_den
      Xpow_p <- outer(Xp, 0:s, "^")
      fitted_p <- Xpow_p %*% alphaT_plus
      residuals_p <- T[idx_plus] - fitted_p
      V_plus_den <- sum(wT_plus[idx_plus] * residuals_p^2) / sum(wT_plus[idx_plus])
    } else {
      V_plus_den <- pilot_var_T_plus  # Fallback
    }
    
    # Minus side
    idx_minus <- which(wT_minus > 0)
    if (length(idx_minus) > 0) {
      Xm <- Xc[idx_minus] / pilot_h_den
      Xpow_m <- outer(Xm, 0:s, "^")
      fitted_m <- Xpow_m %*% alphaT_minus
      residuals_m <- T[idx_minus] - fitted_m
      V_minus_den <- sum(wT_minus[idx_minus] * residuals_m^2) / sum(wT_minus[idx_minus])
    } else {
      V_minus_den <- pilot_var_T_minus  # Fallback
    }
  }
  
  # ------------------------------------------------
  # Step 3: Final MSE/IMSE-Optimal Bandwidth
  # ------------------------------------------------
  
  # Simple (per-quantile) MSE-optimal bandwidths
  if (method == "simple") {
    h_star_num <- numeric(nQ)
    for (qi in seq_len(nQ)) {
      bias_diff <- B_plus[qi] - B_minus[qi]
      var_sum <- (V_plus[qi] + V_minus[qi]) / f_X_hat
      
      if (abs(bias_diff) > 1e-14) {
        ratio <- var_sum / (2 * (s + 1) * bias_diff^2)
        h_star_num[qi] <- ratio^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
      } else {
        # Fallback to Silverman's rule if bias difference is effectively zero
        h_star_num[qi] <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
      }
      
      # Ensure positive bandwidth
      if (is.na(h_star_num[qi]) || h_star_num[qi] <= 0) {
        h_star_num[qi] <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
      }
    }
  } else {
    # Frechet method: IMSE-optimal bandwidth
    dq <- diff(range(q_grid)) / (length(q_grid) - 1)
    A_s <- sum((B_plus - B_minus)^2) * dq
    B_s <- sum((V_plus + V_minus) / f_X_hat) * dq
    
    if (A_s > 1e-14) {
      ratio <- B_s / (2 * (s + 1) * A_s)
      h_star_num <- ratio^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
    } else {
      h_star_num <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
    }
    
    # Ensure positive bandwidth
    if (is.na(h_star_num) || h_star_num <= 0) {
      h_star_num <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
    }
  }
  
  # Treatment bandwidth for fuzzy design
  if (fuzzy) {
    bias_diff_den <- B_plus_den - B_minus_den
    var_sum_den <- (V_plus_den + V_minus_den) / f_X_hat
    
    if (abs(bias_diff_den) > 1e-14) {
      ratio_den <- var_sum_den / (2 * (s + 1) * bias_diff_den^2)
      h_star_den <- ratio_den^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
    } else {
      h_star_den <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
    }
    
    # Ensure positive bandwidth
    if (is.na(h_star_den) || h_star_den <= 0) {
      h_star_den <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
    }
  } else {
    h_star_den <- NULL
  }
  
  if (coverage) {
    ## ROT for coverage error (Calonico et al 2020/2018)
    h_star_num <- h_star_num * n^{-s/((2*s+3)*(s+3))}
  
    if (fuzzy) {
      h_star_den <- h_star_den * n^{-s/((2*s+3)*(s+3))}
    }
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