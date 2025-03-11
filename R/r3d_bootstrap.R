#' Multiplier Bootstrap for Distributional RDD
#'
#' Performs a multiplier bootstrap to obtain uniform confidence bands 
#' and (optionally) conduct hypothesis tests for an \emph{R3D} or \emph{F3D} design.
#' It reuses the partial-sum intercept weights and residuals from a fitted \code{r3d} object
#' to avoid re-estimating local polynomial regressions within each bootstrap loop.
#'
#' @param object An S3 object of class \code{"r3d"}, typically the output of \code{\link{r3d}}.
#' @param X Numeric vector of the running variable (the same one used in \code{\link{r3d}}).
#' @param Y_list A list of numeric vectors; each element in this list is the sample from the
#'   outcome distribution of one unit (same data passed to \code{\link{r3d}}).
#' @param T (Optional) numeric or logical vector of treatment statuses for the fuzzy design;
#'   used only if the original \code{r3d} call was fuzzy. Otherwise, can be \code{NULL}.
#' @param B Integer, number of bootstrap draws. Defaults to 200.
#' @param alpha Significance level for uniform confidence bands. Defaults to 0.05.
#' @param test Character indicating which hypothesis test to conduct:
#'   \describe{
#'     \item{\code{"none"}}{No test, only compute confidence bands.}
#'     \item{\code{"nullity"}}{Tests the null \eqn{H_0: \tau(q) = 0 \text{ for all } q}.}
#'     \item{\code{"homogeneity"}}{Tests the null \eqn{H_0: \tau(q) \text{ is constant in } q}.}
#'   }
#' @param cores Number of CPU cores used for parallel computation of bootstrap draws (default 1).
#' @param seed Optional integer to set a random seed for the multiplier draws (for reproducibility).
#' @param ... Unused additional arguments (for compatibility).
#'
#' @details
#' This function implements the multiplier bootstrap approach for distributional RD:
#' it draws i.i.d. normal multipliers \eqn{\xi_i}, re-scales the stored residual partial sums,
#' and reconstructs approximate realizations of the limiting process. 
#' The maximum deviation across quantiles in each realization
#' is then used to form uniform confidence bands and test statistics.
#'
#' For the fuzzy design, the ratio-based estimator is handled similarly, using
#' the same partial-sum logic for the treatment variable.
#'
#' @return A list with the elements:
#' \describe{
#'   \item{\code{cb_lower}, \code{cb_upper}}{Numeric vectors giving the lower and upper
#'     uniform confidence bands at each quantile in \code{object$q_grid}.}
#'   \item{\code{boot_taus}}{A matrix of bootstrap-draw realizations of the entire \eqn{\tau(q)}.}
#'   \item{\code{supvals}}{For each bootstrap draw, the supremum (max) absolute deviation.}
#'   \item{\code{crit_val}}{The critical value (e.g., the \code{(1 - alpha)} quantile of \code{supvals}).}
#'   \item{\code{test_stat}, \code{test_crit_val}, \code{p_value}}{If \code{test != "none"},
#'     these store the test statistic, critical value, and a simple bootstrap-based p-value.}
#' }
#'
#' @seealso \code{\link{r3d}}, \code{\link{plot.r3d}}, \code{\link{summary.r3d}}
#'
#' @export
r3d_bootstrap <- function(object, X, Y_list, T = NULL,
                          B = 200, alpha = 0.05,
                          test = c("none", "nullity", "homogeneity"),
                          cores = 1, seed = NULL, ...)
{
  test <- match.arg(test)
  if (!inherits(object, "r3d")) stop("Need r3d object from r3d()")
  
  # Extract components from object
  n <- length(X)
  q_grid <- object$q_grid
  nQ <- length(q_grid)
  e1_mat <- object$e1_mat
  e2_mat <- object$e2_mat
  tauhat <- object$tau
  fuzzy <- object$fuzzy
  method <- object$method
  w_plus <- object$w_plus
  w_minus <- object$w_minus
  h_star_num <- object$bandwidths$h_star_num
  h_star_den <- if (fuzzy) object$bandwidths$h_star_den else NULL
  
  # Check if h_star_num is a vector (for "simple" method)
  is_vector_h_num <- length(h_star_num) == nQ && method == "simple"
  
  # For fuzzy RDD, extract additional components
  if (fuzzy) {
    alphaT_plus <- object$alphaT_plus
    alphaT_minus <- object$alphaT_minus
    int_plus <- object$int_plus
    int_minus <- object$int_minus
    denomT <- alphaT_plus[1, 1] - alphaT_minus[1, 1]
    if (abs(denomT) < 1e-14) {
      message("Bootstrap: fuzzy denominator near 0 => might blow up")
    }
    num_diff <- int_plus - int_minus
  }
  
  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  # Center X values
  X_centered <- X - object$cutoff
  
  # Estimate f_X(0) using Silverman's rule
  sigma_X <- stats::sd(X_centered)
  h_bw <- 1.06 * sigma_X * n^(-1/5)
  kernel <- object$kernel
  f_X_hat <- mean(kernel(X_centered / h_bw)) / h_bw
  
  # Pre-compute matrices for efficiency
  e1_w_plus <- e1_mat * w_plus
  e1_w_minus <- e1_mat * w_minus
  if (fuzzy) {
    e2_w_plus <- e2_mat * w_plus
    e2_w_minus <- e2_mat * w_minus
  }
  
  # Define function to process one bootstrap draw
  doOneDraw <- function(bi) {
    xi <- rnorm(n)
    
    if (is_vector_h_num) {
      # "simple" method: vector of bandwidths (one per quantile)
      plus_sums <- numeric(nQ)
      minus_sums <- numeric(nQ)
      for (q in 1:nQ) {
        h_q <- h_star_num[q]
        scaling_q <- 1 / (sqrt(n * h_q) * f_X_hat)
        plus_sums[q] <- sum(xi * e1_w_plus[, q]) * scaling_q
        minus_sums[q] <- sum(xi * e1_w_minus[, q]) * scaling_q
      }
      out_sharp <- plus_sums - minus_sums
    } else {
      # "frechet" method: single bandwidth
      h_num <- if (length(h_star_num) == 1) h_star_num else mean(h_star_num, na.rm = TRUE)
      scaling_num <- 1 / (sqrt(n * h_num) * f_X_hat)
      plus_sums <- colSums(xi * e1_w_plus) * scaling_num
      minus_sums <- colSums(xi * e1_w_minus) * scaling_num
      out_sharp <- plus_sums - minus_sums
    }
    
    if (!fuzzy) {
      return(out_sharp)
    } else {
      # Fuzzy RDD: denominator uses h_star_den
      scaling_den <- 1 / (sqrt(n * h_star_den) * f_X_hat)
      plus_sums2 <- sum(xi * e2_w_plus[, 1]) * scaling_den
      minus_sums2 <- sum(xi * e2_w_minus[, 1]) * scaling_den
      top <- denomT * out_sharp - num_diff * (plus_sums2 - minus_sums2)
      return(top / (denomT^2))
    }
  }
  
  # Run bootstrap in parallel or serially
  if (cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cl))
      boot_list <- parallel::parLapply(cl, 1:B, function(i) doOneDraw(i))
    } else {
      boot_list <- parallel::mclapply(1:B, doOneDraw, mc.cores = cores)
    }
  } else {
    boot_list <- lapply(1:B, doOneDraw)
  }
  
  # Convert list to matrix (columns are bootstrap draws)
  boot_mat <- do.call(cbind, boot_list)
  
  # Calculate uniform confidence bands
  supvals <- apply(boot_mat, 2, function(colb) max(abs(colb), na.rm = TRUE))
  cval <- stats::quantile(supvals, probs = 1 - alpha, na.rm = TRUE)
  cb_lower <- tauhat - cval
  cb_upper <- tauhat + cval
  
  # Test statistics
  test_stat <- NA_real_
  test_crit <- NA_real_
  p_val <- NA_real_
  
  if (test == "nullity") {
    test_stat <- max(abs(tauhat), na.rm = TRUE)
    supvals_null <- apply(boot_mat, 2, function(colb) max(abs(colb), na.rm = TRUE))
    test_crit <- stats::quantile(supvals_null, 1 - alpha, na.rm = TRUE)
    p_val <- mean(supvals_null >= test_stat, na.rm = TRUE)
  } else if (test == "homogeneity") {
    mbar <- mean(tauhat, na.rm = TRUE)
    test_stat <- max(abs(tauhat - mbar), na.rm = TRUE)
    boot_means <- colMeans(boot_mat, na.rm = TRUE)
    boot_centered <- sweep(boot_mat, 2, boot_means)
    supvals_homo <- apply(boot_centered, 2, function(col) max(abs(col), na.rm = TRUE))
    test_crit <- stats::quantile(supvals_homo, 1 - alpha, na.rm = TRUE)
    p_val <- mean(supvals_homo >= test_stat, na.rm = TRUE)
  }
  
  # Return results
  list(cb_lower = cb_lower, 
       cb_upper = cb_upper,
       boot_taus = boot_mat,
       supvals = supvals,
       crit_val = cval,
       test_stat = test_stat,
       test_crit_val = test_crit,
       p_value = p_val)
}