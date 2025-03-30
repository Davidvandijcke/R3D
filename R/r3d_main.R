#' Regression Discontinuity Design with Distribution-Valued Outcomes (R3D)
#'
#' Fits a regression discontinuity design with \emph{distribution-valued} outcomes (random distributions),
#' as developed in \insertCite{vandijcke2025;textual}{R3D},
#' using local polynomial regression on quantiles (\code{method="simple"}) or local Fréchet regression
#' (\code{method="frechet"}). Supports both sharp and fuzzy designs.
#'
#' @param X Numeric vector of running variable values (length \eqn{n}).
#' @param Y_list A list of length \eqn{n}, where each element is a numeric vector
#'   representing the sample from the outcome distribution of one unit.
#' @param T (Optional) For a fuzzy design: numeric or logical vector of length \eqn{n} 
#'   giving partial treatment status (0 or 1). If \code{fuzzy=FALSE}, this can be omitted.
#' @param cutoff Numeric scalar: the threshold for treatment assignment. The design
#'   treats \eqn{X \ge cutoff} as above the cutoff.
#' @param method Either \code{"simple"} or \code{"frechet"}:
#'   \describe{
#'     \item{\code{"simple"}}{Local polynomial regression is done \emph{pointwise} for each quantile,
#'       and then optionally rearranged to maintain monotonicity.}
#'     \item{\code{"frechet"}}{A global local Fréchet regression approach is done, projecting
#'       the entire curve onto quantile functions.}
#'   }
#' @param p Integer specifying the local polynomial order (e.g., 2 for local quadratic).
#' @param q_grid Numeric vector of quantiles at which to estimate the distributional treatment effect.
#' @param fuzzy Logical flag for fuzzy design. Default \code{FALSE} for a sharp design.
#' @param kernel_fun Name of the kernel function. Possible choices are
#'   \code{"epanechnikov"}, \code{"triangular"}, or \code{"uniform"}. 
#'   Defaults to \code{"triangular"}.
#' @param s Integer expansion order in the pilot local polynomial step used for bandwidth
#'   selection. Defaults to 1.
#' @param bandwidths Optional user-specified bandwidths. If \code{NULL} (default), bandwidths are 
#'   automatically selected via \code{\link{r3d_bwselect}}. For a sharp design, can be
#'   either a single numeric value or a vector of length \code{length(q_grid)}. For a fuzzy design,
#'   must be a list with two components: \code{num} (for the outcome equation) and \code{den} (for the 
#'   first-stage equation). Where \code{num} can be a scalar or vector of length \code{length(q_grid)}, 
#'   and \code{den} must be a scalar.
#' @param boot Logical. If \code{TRUE}, automatically runs \code{\link{r3d_bootstrap}}
#'   after estimation to produce uniform confidence bands and tests.
#' @param boot_reps Integer, number of bootstrap draws if \code{boot=TRUE}.
#' @param boot_cores Number of CPU cores for parallelizing the bootstrap. Default 1.
#' @param alpha Significance level for the uniform confidence bands. Default 0.05.
#' @param test Character vector of tests to perform. Options include: \code{"none"}, \code{"nullity"}, 
#'   \code{"homogeneity"}, or \code{"gini"}, or a vector combining multiple test types. When \code{"none"}, 
#'   no tests are performed. When \code{"nullity"}, tests the null hypothesis that the treatment effect 
#'   is zero across all quantiles. When \code{"homogeneity"}, tests the null hypothesis that the 
#'   treatment effect is constant across all quantiles. When \code{"gini"}, tests the null hypothesis 
#'   that there is no difference in the Gini coefficients of the distributions above and below the cutoff.
#' @param test_ranges List of numeric vectors defining the quantile ranges for testing. Each element should be a
#'   vector of length 2 or more defining the ranges. For example, \code{list(c(0.25, 0.75))} to test on 
#'   quantiles between 0.25 and 0.75, or \code{list(c(0.25, 0.5, 0.75))} to test on ranges [0.25, 0.5] and [0.5, 0.75].
#'   If \code{NULL} (default), tests are performed on the entire \code{q_grid}.
#' @param coverage Logical indicating whether to apply the coverage correction rule of thumb of
#' \insertCite{calonico2018effect;textual}{R3D}. Default is FALSE.
#' @param ... Additional arguments passed to the bandwidth selection function \code{\link{r3d_bwselect}}.
#'
#' @details
#' This is the main user-facing function for the R3D approach. Internally, it:
#' \enumerate{
#'   \item Selects bandwidth(s) via \code{\link{r3d_bwselect}} (if \code{bandwidths=NULL}) or 
#'         uses user-provided bandwidths,
#'   \item Computes local polynomial fits on each side of the cutoff for each quantile in \code{q_grid},
#'   \item For \code{method="frechet"}, performs a projection of the fitted quantile curves onto the
#'         space of valid quantile functions (via isotonic regression),
#'   \item (Optionally) runs the multiplier bootstrap for uniform inference.
#' }
#'
#' When \code{test} includes one or more test types and \code{test_ranges} is specified, the function will 
#' perform the specified tests on each range in \code{test_ranges}. For example, if \code{test=c("nullity", "homogeneity")} 
#' and \code{test_ranges=list(c(0.25, 0.75), c(0.1, 0.9))}, four tests will be performed: nullity and homogeneity 
#' tests on the range [0.25, 0.75], and nullity and homogeneity tests on the range [0.1, 0.9].
#' 
#' The Gini test (\code{test="gini"}) examines whether there is a statistically significant difference in 
#' inequality (as measured by the Gini coefficient) between the distributions above and below the cutoff.
#' It uses the estimated conditional quantile functions to calculate the Gini coefficients.
#'
#' @return An S3 object of class \code{"r3d"}, containing (among others):
#' \describe{
#'   \item{\code{tau}}{Estimated distributional treatment effect \eqn{\tau(q)} for each \eqn{q} in \code{q_grid}.}
#'   \item{\code{q_grid}}{The quantile grid used.}
#'   \item{\code{bandwidths}}{The bandwidths used (either user-provided or MSE-/IMSE-optimal).}
#'   \item{\code{w_plus}, \code{w_minus}}{The internal local polynomial weights on the plus/minus side.}
#'   \item{\code{e1_mat}, \code{e2_mat}}{Residual matrices used for the multiplier bootstrap.}
#'   \item{\code{boot_out}}{If \code{boot=TRUE}, a list of bootstrap results from \code{\link{r3d_bootstrap}}.}
#' }
#'
#' @seealso \code{\link{r3d_bootstrap}}, \code{\link{plot.r3d}}, \code{\link{summary.r3d}}
#'
#' @examples
#' \dontrun{
#'   # Simulate data
#'   set.seed(123)
#'   n <- 100
#'   X <- runif(n, -1, 1)
#'   Y_list <- lapply(seq_len(n), function(i) {
#'     # Suppose distribution is Normal with mean depending on X
#'     rnorm(sample(30:50,1), mean=2 + 2*(X[i]>=0))
#'   })
#'   T <- as.numeric(X >= 0) * rbinom(n, 1, 0.8)  # Fuzzy treatment
#'
#'   # Example 1: Test nullity and homogeneity on full range
#'   fit1 <- r3d(X, Y_list, T=T, cutoff=0,
#'               method="frechet", p=2, fuzzy=TRUE,
#'               boot=TRUE, boot_reps=200, alpha=0.05, 
#'               test=c("nullity", "homogeneity"))
#'
#'   # Example 2: Test Gini coefficient difference
#'   fit2 <- r3d(X, Y_list, cutoff=0, method="simple", 
#'               boot=TRUE, test="gini")
#'               
#'   # Example 3: Test on multiple ranges
#'   fit3 <- r3d(X, Y_list, cutoff=0, 
#'               boot=TRUE, test=c("nullity", "homogeneity", "gini"), 
#'               test_ranges=list(c(0.25, 0.5, 0.75)))
#'
#'   # Inspect results
#'   print(fit1)
#'   summary(fit1)
#'   plot(fit1)
#' }
#'
#' @useDynLib R3D, .registration = TRUE
#' 
#' @references
#' \insertAllCited{}
#' 
#' @export
r3d <- function(X, Y_list, T = NULL,
                cutoff = 0,
                method = c("simple", "frechet"),
                p = 2,
                q_grid = seq(0.01, 0.99, 0.01),
                fuzzy = FALSE,
                kernel_fun = c("epanechnikov", "triangular", "uniform"),
                s = 1,
                bandwidths = NULL,
                boot = FALSE,
                boot_reps = 200,
                boot_cores = 1,
                alpha = 0.05,
                test = c("none", "nullity", "homogeneity", "gini"),
                test_ranges = NULL,
                coverage = FALSE,
                ...)
{
  method <- match.arg(method)
  kernel_fun <- match.arg(kernel_fun)
  
  # Handle multiple test options
  if (is.character(test) && length(test) == 1) {
    test <- match.arg(test, choices = c("none", "nullity", "homogeneity", "gini"))
  } else if (is.character(test) && length(test) > 1) {
    valid_tests <- c("none", "nullity", "homogeneity", "gini")
    if (!all(test %in% valid_tests)) {
      stop("Invalid test option(s). Must be one or more of: 'none', 'nullity', 'homogeneity', 'gini'")
    }
    # Remove "none" if other options are present
    if ("none" %in% test && length(test) > 1) {
      test <- setdiff(test, "none")
    }
  } else {
    stop("'test' must be a character vector of test types")
  }
  
  # Process test_ranges
  if (!is.null(test_ranges)) {
    if (!is.list(test_ranges)) {
      stop("'test_ranges' must be a list of numeric vectors")
    }
    # Validate each range
    for (i in seq_along(test_ranges)) {
      range <- test_ranges[[i]]
      if (!is.numeric(range) || length(range) < 2) {
        stop("Each element in 'test_ranges' must be a numeric vector of length >= 2")
      }
      if (any(range < 0 | range > 1)) {
        stop("All values in 'test_ranges' must be between 0 and 1")
      }
      if (min(range) < min(q_grid) | max(range) > max(q_grid)) {
        stop("All values in 'test_ranges' must be within the range of q_grid")
      }
      
      
      # Get min and max of the range
      min_range <- min(range)
      max_range <- max(range)
      
      # Find all q_grid values within this range
      range_values <- q_grid[q_grid >= min_range & q_grid <= max_range]
      
      if (length(range_values) < 2) {
        stop("Range [", min_range, ", ", max_range, "] contains fewer than 2 points from q_grid. Ensure range includes sufficient values.")
      }
      
      # Sort the range but keep original endpoints to maintain user intention
      test_ranges[[i]] <- sort(range)
    }
  }
  
  kernel <- switch(kernel_fun,
                   triangular = function(u) pmax(0, 1 - abs(u)),
                   epanechnikov = function(u) 0.75 * pmax(0, 1 - u^2),
                   uniform = function(u) 0.5 * (abs(u) <= 1),
                   stop("Unknown kernel type."))
  
  # Validate inputs
  validate_r3d_inputs(X, Y_list, T, cutoff, method, as.integer(p), q_grid, fuzzy, 
                      kernel_fun, as.integer(s), boot, as.integer(boot_reps), as.integer(boot_cores), alpha, test, test_ranges)
  
  n <- length(X)
  
  nQ <- length(q_grid)
  
  # 1) Bandwidth selection 
  
  # Check bandwidths if provided
  if (!is.null(bandwidths)) {
    if (is.list(bandwidths)) {
      # For fuzzy RDD: list with num and den components
      if (fuzzy) {
        h_star_num <- bandwidths$num
        h_star_den <- bandwidths$den
        
        # Validate numerator bandwidths
        if (is.null(h_star_num)) {
          stop("For fuzzy RDD, bandwidths$num must be provided.")
        }
        if (!(length(h_star_num) == 1 || length(h_star_num) == nQ)) {
          stop("bandwidths$num must be a scalar or a vector of length equal to q_grid.")
        }
        
        # Validate denominator bandwidth
        if (is.null(h_star_den)) {
          stop("For fuzzy RDD, bandwidths$den must be provided.")
        }
        if (length(h_star_den) != 1) {
          stop("bandwidths$den must be a scalar.")
        }
      } else {
        # Sharp RDD: extract numerator only
        h_star_num <- bandwidths$num
        h_star_den <- NULL
        
        if (is.null(h_star_num)) {
          stop("bandwidths$num must be provided.")
        }
        if (!(length(h_star_num) == 1 || length(h_star_num) == nQ)) {
          stop("bandwidths$num must be a scalar or a vector of length equal to q_grid.")
        }
      }
    } else {
      # Direct bandwidth specification (for sharp RDD)
      if (fuzzy) {
        stop("For fuzzy RDD, bandwidths must be a list with 'num' and 'den' components.")
      }
      
      h_star_num <- bandwidths
      h_star_den <- NULL
      
      if (!(length(h_star_num) == 1 || length(h_star_num) == nQ)) {
        stop("bandwidths must be a scalar or a vector of length equal to q_grid.")
      }
    }
  } else {
    # 1) Use automatic bandwidth selection if no bandwidths provided
    bwres <- r3d_bwselect(X = X,
                          Y_list = Y_list,
                          T = T,  # Pass T for fuzzy design
                          q_grid = q_grid,
                          method = method,
                          s = s,
                          p = p,
                          kernel = kernel,
                          cutoff = cutoff,
                          fuzzy = fuzzy,  # Pass fuzzy flag
                          coverage = coverage,
                          ...)
    
    # Extract bandwidths
    h_star_num <- bwres$h_star_num
    h_star_den <- if (fuzzy) bwres$h_star_den else NULL
  }
  
  # 2) Build matrix of empirical quantiles [n x nQ]
  # check that Y_list is valid (not empty etc)
  
  
  Qmat <- .compute_empirical_qmat(Y_list, q_grid)
  
  # 3) Handle bandwidths based on method -- MODIFIED
  if (method == "simple") {
    h_use_num <- h_star_num  # Vector of bandwidths, length nQ
  } else {
    h_use_num <- h_star_num  # Scalar for frechet
  }
  h_use_den <- if (fuzzy) h_star_den else NULL  # Scalar for T in fuzzy case
  
  # 4) Re-center
  X_centered <- X - cutoff
  
  # Fortran setup
  N <- as.integer(n)
  P_ <- as.integer(p)
  NQ_ <- as.integer(nQ)
  
  # Define kernel type as integer: 1=triangular, 2=epanechnikov, 3=uniform
  kernel_type <- as.integer(switch(kernel_fun,
                                   triangular = 1,
                                   epanechnikov = 2,
                                   uniform = 3))
  
  # 5) Local poly fit on plus side for outcome quantiles 
  outPlus <- .Fortran("locweights",
                      X = as.double(X_centered),
                      YMAT = as.double(Qmat),
                      N = N,
                      P = P_,
                      H = as.double(if (length(h_use_num) == 1) rep(h_use_num, nQ) else h_use_num),
                      SIDE = as.integer(1),
                      KERNEL_TYPE = kernel_type,
                      ALPHA = double((p + 1) * nQ),
                      WINT = double(n * nQ),
                      INFO = integer(1),
                      NQ = NQ_,
                      PACKAGE = "R3D")
  info_plus <- outPlus$INFO
  alpha_plus <- matrix(outPlus$ALPHA, nrow = p + 1, ncol = nQ)
  w_plus <- matrix(outPlus$WINT, nrow = n, ncol = nQ)
  
  # Minus side for outcome quantiles -- MODIFIED
  outMinus <- .Fortran("locweights",
                       X = as.double(X_centered),
                       YMAT = as.double(Qmat),
                       N = N,
                       P = P_,
                       H = as.double(if (length(h_use_num) == 1) rep(h_use_num, nQ) else h_use_num),
                       SIDE = as.integer(0),
                       KERNEL_TYPE = kernel_type,
                       ALPHA = double((p + 1) * nQ),
                       WINT = double(n * nQ),
                       INFO = integer(1),
                       NQ = NQ_,
                       PACKAGE = "R3D")
  info_minus <- outMinus$INFO
  alpha_minus <- matrix(outMinus$ALPHA, nrow = p + 1, ncol = nQ)
  w_minus <- matrix(outMinus$WINT, nrow = n, ncol = nQ)
  
  if (info_plus != 0) warning("locweights: plus side singular system? info=", info_plus)
  if (info_minus != 0) warning("locweights: minus side singular system? info=", info_minus)
  
  # 6) If fuzzy, fit T with h_use_den
  alphaT_plus <- NULL
  alphaT_minus <- NULL
  denom_T <- NA_real_
  
  if (fuzzy) {
    T_mat <- matrix(T, nrow = n, ncol = 1)
    
    # Plus side for T
    outTplus <- .Fortran("locweights",
                         X = as.double(X_centered),
                         YMAT = as.double(T_mat),
                         N = N,
                         P = P_,
                         H = as.double(h_use_den),  # Use h_use_den
                         SIDE = as.integer(1),
                         KERNEL_TYPE = kernel_type,
                         ALPHA = double((p + 1)),
                         WINT = double(n),
                         INFO = integer(1),
                         NQ = as.integer(1),
                         PACKAGE = "R3D")
    alphaT_plus <- matrix(outTplus$ALPHA, nrow = p + 1, ncol = 1)
    
    # Minus side for T
    outTminus <- .Fortran("locweights",
                          X = as.double(X_centered),
                          YMAT = as.double(T_mat),
                          N = N,
                          P = P_,
                          H = as.double(h_use_den),  # Use h_use_den
                          SIDE = as.integer(0),
                          KERNEL_TYPE = kernel_type,
                          ALPHA = double((p + 1)),
                          WINT = double(n),
                          INFO = integer(1),
                          NQ = as.integer(1),
                          PACKAGE = "R3D")
    alphaT_minus <- matrix(outTminus$ALPHA, nrow = p + 1, ncol = 1)
    
    denom_T <- alphaT_plus[1, 1] - alphaT_minus[1, 1]
    if (abs(denom_T) < 1e-14) {
      warning("Denominator in fuzzy RDD is near 0 => effects might be NA.")
    }
  }
  
  # 7) Compute unprojected fits E[Y(q)|X=x_i], row by row 
  Eplus_unproj <- matrix(0, nrow = n, ncol = nQ)  # Storage for plus side fits
  Eminus_unproj <- matrix(0, nrow = n, ncol = nQ) # Storage for minus side fits
  
  for (qi in seq_len(nQ)) {
    # Extract the bandwidth for this quantile
    h_q <- if (length(h_use_num) == 1) h_use_num else h_use_num[qi]
    
    # Scale X_centered for this quantile
    X_scaled_q <- X_centered / h_q
    
    # Compute fits for the plus side
    ap <- alpha_plus[, qi]  # Coefficients for this quantile (plus side)
    idxp <- which(w_plus[, qi] > 0)  # Indices where weights are positive
    if (length(idxp) > 0) {
      Xpow_p <- outer(X_scaled_q[idxp], 0:p, "^")  # Polynomial basis
      Eplus_unproj[idxp, qi] <- Xpow_p %*% ap      # Fitted values
    }
    
    # Compute fits for the minus side
    am <- alpha_minus[, qi]  # Coefficients for this quantile (minus side)
    idxm <- which(w_minus[, qi] > 0)  # Indices where weights are positive
    if (length(idxm) > 0) {
      Xpow_m <- outer(X_scaled_q[idxm], 0:p, "^")  # Polynomial basis
      Eminus_unproj[idxm, qi] <- Xpow_m %*% am     # Fitted values
    }
  }
  
  # 8) If frechet, perform row-wise isotonic regression
  Eplus_final <- Eplus_unproj
  Eminus_final <- Eminus_unproj
  
  if (method == "frechet") {
    for (i in seq_len(n)) {
      if (any(w_plus[i, ] > 0)) {
        rowfit <- Eplus_unproj[i, ]
        isoobj <- stats::isoreg(q_grid, rowfit)
        Eplus_final[i, ] <- as.numeric(isoobj$yf)
      }
      if (any(w_minus[i, ] > 0)) {
        rowfit <- Eminus_unproj[i, ]
        isoobj <- stats::isoreg(q_grid, rowfit)
        Eminus_final[i, ] <- as.numeric(isoobj$yf)
      }
    }
  }
  
  # 9) Build outcome residuals
  e1_mat <- matrix(0, nrow = n, ncol = nQ)
  for (qi in seq_len(nQ)) {
    idxp <- which(w_plus[, qi] > 0)
    if (length(idxp) > 0) {
      e1_mat[idxp, qi] <- Qmat[idxp, qi] - Eplus_final[idxp, qi]
    }
    idxm <- which(w_minus[, qi] > 0)
    if (length(idxm) > 0) {
      e1_mat[idxm, qi] <- Qmat[idxm, qi] - Eminus_final[idxm, qi]
    }
  }
  
  # 10) If fuzzy, build T residuals
  e2_mat <- NULL
  if (fuzzy) {
    e2_mat <- matrix(0, nrow = n, ncol = nQ)
    
    # Plus side for T
    idxp <- which(w_plus[, 1] > 0)  # Using first column as proxy
    if (length(idxp) > 0) {
      X_scaled_den <- X_centered[idxp] / h_use_den
      Xpow_p <- outer(X_scaled_den, 0:p, "^")
      fitpT <- Xpow_p %*% alphaT_plus
      residp <- T[idxp] - fitpT
      e2_mat[idxp, ] <- matrix(rep(residp, nQ), ncol = nQ)
    }
    
    # Minus side for T
    idxm <- which(w_minus[, 1] > 0)
    if (length(idxm) > 0) {
      X_scaled_den <- X_centered[idxm] / h_use_den
      Xpow_m <- outer(X_scaled_den, 0:p, "^")
      fitmT <- Xpow_m %*% alphaT_minus
      residm <- T[idxm] - fitmT
      e2_mat[idxm, ] <- matrix(rep(residm, nQ), ncol = nQ)
    }
  }
  
  # 11) Compute final LAQTE from intercepts
  int_plus <- alpha_plus[1, ]
  int_minus <- alpha_minus[1, ]
  
  if (method == "frechet") {
    iso_p <- stats::isoreg(q_grid, int_plus)
    int_plus <- as.numeric(iso_p$yf)
    iso_m <- stats::isoreg(q_grid, int_minus)
    int_minus <- as.numeric(iso_m$yf)
  } else {
    int_plus <- Rearrangement::rearrangement(list(q_grid), int_plus)
    int_minus <- Rearrangement::rearrangement(list(q_grid), int_minus)
  }
  
  tau_vec <- int_plus - int_minus
  if (fuzzy) {
    if (abs(denom_T) < 1e-14) {
      tau_vec[] <- NA_real_
    } else {
      tau_vec <- tau_vec / denom_T
    }
  }
  
  # 12) Build the output 
  out <- list(
    results = list(
      tau = tau_vec,
      q_grid = q_grid,
      method = method,
      fuzzy = fuzzy,
      p = p,
      bandwidths = list(h_star_num = h_star_num, h_star_den = h_star_den),  
      h_used = list(h_use_num = h_use_num, h_use_den = h_use_den)  
    ),
    coefficients = list(
      alpha_plus = alpha_plus,
      alpha_minus = alpha_minus,
      alphaT_plus = alphaT_plus,
      alphaT_minus = alphaT_minus
    ),
    bootstrap = list(
      w_plus = w_plus,
      w_minus = w_minus,
      e1_mat = e1_mat,
      e2_mat = e2_mat
    ),
    conditional_means = list(
      plus = int_plus,
      minus = int_minus
    ),
    inputs = list(
      X = X,
      Y_list = Y_list,
      T = T,
      cutoff = cutoff,
      call = match.call()
    ),
    diagnostics = list(
      info_plus = info_plus,
      info_minus = info_minus,
      denominator = if (fuzzy) denom_T else NULL
    )
  )
  
  # Replicate top-level items for convenience 
  out$tau <- tau_vec
  out$q_grid <- q_grid
  out$method <- method
  out$fuzzy <- fuzzy
  out$p <- p
  out$bandwidths <- list(h_star_num = h_star_num, h_star_den = h_star_den)
  out$kernel <- kernel
  out$w_plus <- w_plus
  out$w_minus <- w_minus
  out$alpha_plus <- alpha_plus
  out$alpha_minus <- alpha_minus
  out$alphaT_plus <- alphaT_plus
  out$alphaT_minus <- alphaT_minus
  out$int_plus <- int_plus
  out$int_minus <- int_minus
  out$w_plus <- w_plus
  out$w_minus <- w_minus
  out$e1_mat <- e1_mat
  out$e2_mat <- e2_mat
  out$X <- X
  out$Y_list <- Y_list
  out$T <- T
  out$cutoff <- cutoff
  out$call <- match.call()
  
  class(out) <- "r3d"
  
  # Optional bootstrap
  if (boot) {
    boot_out <- r3d_bootstrap(
      object = out,
      X = X,
      Y_list = Y_list,
      T = T,
      B = boot_reps,
      alpha = alpha,
      test = test,
      test_ranges = test_ranges,
      cores = boot_cores
    )
    out$bootstrap$boot_out <- boot_out
    out$boot_out <- boot_out
  }
  
  out
}