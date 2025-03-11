# File: R3D/R/r3d_bootstrap.R

#' @title r3d_bootstrap: Multiplier Bootstrap for Uniform Inference
#'
#' @description
#' Performs a multiplier bootstrap to obtain uniform confidence bands 
#' and conduct hypothesis tests in \emph{distributional} RD settings. 
#' Reuses the partial-sum intercept weights and residuals from a fitted \code{r3d} object,
#' thus avoiding recomputing local polynomial regressions inside the bootstrap loop.
#'
#' @param object An S3 object of class \code{"r3d"}, typically the output of \code{\link{r3d}}.
#' @param X Numeric vector (same as in the original call).
#' @param Y_list List of numeric vectors (same as in the original call).
#' @param T (Optional) numeric or logical vector for the fuzzy design; same data used in the original \code{r3d} call.
#' @param B Integer, number of bootstrap draws (default 200).
#' @param alpha Significance level for uniform confidence bands (default 0.05).
#' @param test Character, either \code{"none"}, \code{"nullity"}, or \code{"homogeneity"}. 
#'   \itemize{
#'     \item \code{"none"}: no test, just compute confidence bands.
#'     \item \code{"nullity"}: test the null that \eqn{\tau(q) = 0} for all \eqn{q}.
#'     \item \code{"homogeneity"}: test the null that \eqn{\tau(q)} is constant across \eqn{q}.
#'   }
#' @param cores Number of CPU cores for parallel computing of bootstrap draws (default 1).
#' @param seed Optional integer to set random seed for the multiplier draws.
#' @param ... Unused additional arguments.
#'
#' @details
#' \bold{Approach}:
#' A multiplier bootstrap draws i.i.d. normal multipliers \eqn{\xi_i}, 
#' and re-scales the residual partial sums to generate approximate realizations 
#' of the limiting process. See the references for theoretical details.
#'
#' The resulting object can be used to construct uniform confidence bands and test statistics.
#'
#' \bold{Tests}:
#' \describe{
#'   \item{\code{test = "nullity"}}{Tests \eqn{H_0: \tau(q)=0~\forall q} vs. \eqn{H_a: \tau(q)\neq 0} for some \eqn{q}.}
#'   \item{\code{test = "homogeneity"}}{Tests \eqn{H_0: \tau(q) ~\text{is constant in}~ q} vs. \eqn{H_a: \tau(q) \text{ not constant}.}
#' }
#'
#' @return A list with components:
#' \item{cb_lower, cb_upper}{Numeric vectors of lower/upper confidence band values, same length as \code{object$q_grid}.}
#' \item{boot_taus}{A matrix of the B bootstrap draws of the entire \eqn{\tau(q)} function. 
#'    (Rows or columns, depending on your internal structure.)}
#' \item{supvals}{For each bootstrap sample, the max absolute deviation from the point estimate.}
#' \item{crit_val}{The critical value for the uniform band, e.g. the (1-\code{alpha}) quantile of \code{supvals}.}
#' \item{test_stat, test_crit_val, p_value}{If \code{test != "none"}, these store the chosen test statistic, 
#'    the critical value, and the resulting p-value.}
#'
#' @references
#' Van Dijcke, D. (2025). \emph{Regression Discontinuity Design with Distributional Outcomes (R3D).} 
#' Working paper. 
#' \cr
#' \cite{chiang2019robust}, \cite{calonico2014robust}, among others, discuss multiplier bootstraps in RD contexts.
#'
#' @seealso 
#' \code{\link{r3d}}, \code{\link{plot.r3d}}, \code{\link{summary.r3d}}
#'
#' @examples
#' \dontrun{
#'   # Suppose you already fit r3d:
#'   fit <- r3d(X, Y_list, boot=FALSE)  # no bootstrap initially
#'   
#'   # Then post hoc, you can run:
#'   bootout <- r3d_bootstrap(fit, X, Y_list, B=300, alpha=0.05, test="homogeneity")
#'   names(bootout)
#'   # you get cb_lower, cb_upper, etc.
#' }
#'
#' @export
r3d_bootstrap <- function(object, X, Y_list, method, T=NULL,
                          B=200, alpha=0.05,
                          test=c("none","nullity","homogeneity"),
                          cores=1, seed=NULL, ...)
{
  test <- match.arg(test)
  if(!inherits(object,"r3d")) stop("Need r3d object from r3d()")
  
  # Extract components from object
  n        <- length(X)
  q_grid   <- object$q_grid
  nQ       <- length(q_grid)
  e1_mat   <- object$e1_mat
  e2_mat   <- object$e2_mat
  tauhat   <- object$tau
  fuzzy    <- object$fuzzy
  method   <- object$method
  w_plus   <- object$w_plus
  w_minus  <- object$w_minus
  
  # For fuzzy RDD, extract additional components
  if(fuzzy) {
    alpha_plus   <- object$alpha_plus
    alpha_minus  <- object$alpha_minus
    alphaT_plus  <- object$alphaT_plus
    alphaT_minus <- object$alphaT_minus
    int_plus <- object$int_plus
    int_minus <- object$int_minus
    denomT       <- alphaT_plus[1,1] - alphaT_minus[1,1]
    
    if(abs(denomT) < 1e-14) {
      message("Bootstrap: fuzzy denominator near 0 => might blow up")
    }
  }
  
  # Set seed if provided
  if(!is.null(seed)) set.seed(seed)
  
  # Pre-compute matrices for efficiency
  # These calculations are independent of the bootstrap draws
  e1_w_plus <- e1_mat * w_plus
  e1_w_minus <- e1_mat * w_minus
  
  if(fuzzy) {
    e2_w_plus <- e2_mat * w_plus
    e2_w_minus <- e2_mat * w_minus
    num_diff <- int_plus - int_minus # alpha_plus[1,] - alpha_minus[1,]
  }
  
  # Define function to process one bootstrap draw
  doOneDraw <- function(bi) {
    # Generate normal random variables
    xi <- rnorm(n)
    
    # For sharp RDD, compute plus and minus partial sums directly
    plus_sums <- colSums(xi * e1_w_plus)
    minus_sums <- colSums(xi * e1_w_minus)
    out_sharp <- plus_sums - minus_sums
    
    if(!fuzzy) {
      return(out_sharp)
    } else {
      # For fuzzy RDD, compute treatment partial sums
      plus_sums2 <- colSums(xi * e2_w_plus)
      minus_sums2 <- colSums(xi * e2_w_minus)
      
      # Calculate numerator for ratio
      top <- denomT * out_sharp - num_diff * (plus_sums2 - minus_sums2)
      
      # Return ratio (handle possible division by zero)
      return(top / (denomT^2))
    }
  }
  
  # Run bootstrap in parallel
  if(cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    # Use parallel processing if cores > 1 and parallel package available
    if(.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cl))
      boot_list <- parallel::parLapply(cl, 1:B, function(i) doOneDraw(i))
    } else {
      boot_list <- parallel::mclapply(1:B, doOneDraw, mc.cores = cores)
    }
  } else {
    # Serial processing
    boot_list <- lapply(1:B, doOneDraw)
  }
  
  # Convert list to matrix (columns are bootstrap draws)
  boot_mat <- do.call(cbind, boot_list)
  
  # Calculate uniform confidence bands
  # Max absolute deviation across quantiles for each bootstrap sample
  supvals <- apply(boot_mat, 2, function(colb) max(abs(colb - tauhat), na.rm = TRUE))
  cval <- stats::quantile(supvals, probs = 1 - alpha, na.rm = TRUE)
  cb_lower <- tauhat - cval
  cb_upper <- tauhat + cval
  
  # Test statistics
  test_stat <- NA_real_
  test_crit <- NA_real_
  p_val <- NA_real_
  
  if(test == "nullity") {
    # Test H0: tau(q) = 0 for all q
    test_stat <- max(abs(tauhat), na.rm = TRUE)
    supvals_null <- apply(boot_mat, 2, function(colb) max(abs(colb), na.rm = TRUE))
    test_crit <- stats::quantile(supvals_null, 1 - alpha, na.rm = TRUE)
    p_val <- mean(supvals_null >= test_stat)
  } else if(test == "homogeneity") {
    # Test H0: tau(q) = constant for all q
    mbar <- mean(tauhat, na.rm = TRUE)
    test_stat <- max(abs(tauhat - mbar), na.rm = TRUE)
    
    # Vectorized calculation for homogeneity test
    # For each bootstrap column, calculate mean, subtract from column, get max abs deviation
    boot_means <- colMeans(boot_mat, na.rm = TRUE)
    boot_centered <- sweep(boot_mat, 2, boot_means)
    supvals_homo <- apply(boot_centered, 2, function(col) max(abs(col), na.rm = TRUE))
    
    test_crit <- stats::quantile(supvals_homo, 1 - alpha, na.rm = TRUE)
    p_val <- mean(supvals_homo >= test_stat)
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