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
#' @param test Character vector indicating which hypothesis tests to conduct:
#'   \describe{
#'     \item{\code{"none"}}{No test, only compute confidence bands.}
#'     \item{\code{"nullity"}}{Tests the null \eqn{H_0: \tau(q) = 0 \text{ for all } q}.}
#'     \item{\code{"homogeneity"}}{Tests the null \eqn{H_0: \tau(q) \text{ is constant in } q}.}
#'     \item{\code{"gini"}}{Tests the null \eqn{H_0: \text{Gini coefficient above} = \text{Gini coefficient below}}.}
#'   }
#'   Multiple test types can be specified as a vector, e.g., \code{c("nullity", "homogeneity")}.
#' @param test_ranges List of numeric vectors defining the quantile ranges for testing. Each element should be a
#'   vector of length 2 or more defining the ranges. For example, \code{list(c(0.25, 0.75))} to test on 
#'   quantiles between 0.25 and 0.75, or \code{list(c(0.25, 0.5, 0.75))} to test on ranges \\[0.25, 0.5\\] and \\[0.5, 0.75\\].
#'   If \code{NULL} (default), tests are performed on the entire \code{q_grid}.
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
#' When multiple tests and/or multiple range specifications are provided, the function
#' will perform each test on each specified range, returning results for all combinations.
#'
#' The Gini test examines whether there is a statistically significant difference between
#' the Gini coefficients of the estimated quantile functions above and below the cutoff.
#'
#' @return A list with the elements:
#' \describe{
#'   \item{\code{cb_lower}, \code{cb_upper}}{Numeric vectors giving the lower and upper
#'     uniform confidence bands at each quantile in \code{object$q_grid}.}
#'   \item{\code{boot_taus}}{A matrix of bootstrap-draw realizations of the entire \eqn{\tau(q)}.}
#'   \item{\code{supvals}}{For each bootstrap draw, the supremum (max) absolute deviation.}
#'   \item{\code{crit_val}}{The critical value (e.g., the \code{(1 - alpha)} quantile of \code{supvals}).}
#'   \item{\code{test_results}}{If \code{test} is not \code{"none"}, a list of test results for each test
#'     and range combination, containing the test statistic, critical value, and p-value for each test.}
#' }
#'
#' @seealso \code{\link{r3d}}, \code{\link{plot.r3d}}, \code{\link{summary.r3d}}
#'
#' @export
r3d_bootstrap <- function(object, X, Y_list, T = NULL,
                          B = 200, alpha = 0.05,
                          test = c("none", "nullity", "homogeneity", "gini"),
                          test_ranges = NULL,
                          cores = 1, seed = NULL, ...)
{
  # Validate test parameter
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
  
  # Conditional means (estimated quantile functions)
  int_plus <- object$int_plus
  int_minus <- object$int_minus
  
  # Check if h_star_num is a vector (for "simple" method)
  is_vector_h_num <- length(h_star_num) == nQ && method == "simple"
  
  # For fuzzy RDD, extract additional components
  if (fuzzy) {
    alphaT_plus <- object$alphaT_plus
    alphaT_minus <- object$alphaT_minus
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
  
  # Check if we need to store the separate bootstrap components for Gini test
  need_separate_components <- "gini" %in% test
  
  # Define function to process one bootstrap draw
  # Modified to store nu_plus and nu_minus if needed
  doOneDraw <- function(bi) {
    xi <- rnorm(n)
    
    # Initialize variables to store nu_plus and nu_minus if needed
    if (need_separate_components) {
      nu_plus <- numeric(nQ)
      nu_minus <- numeric(nQ)
    }
    
    if (is_vector_h_num) {
      # "simple" method: vector of bandwidths (one per quantile)
      plus_sums <- numeric(nQ)
      minus_sums <- numeric(nQ)
      for (q in 1:nQ) {
        h_q <- h_star_num[q]
        scaling_q <- 1
        plus_sums[q] <- sum(xi * e1_w_plus[, q]) * scaling_q
        minus_sums[q] <- sum(xi * e1_w_minus[, q]) * scaling_q
      }
      
      if (need_separate_components) {
        nu_plus <- plus_sums
        nu_minus <- minus_sums
      }
      
      out_sharp <- plus_sums - minus_sums
    } else {
      # "frechet" method: single bandwidth
      h_num <- if (length(h_star_num) == 1) h_star_num else mean(h_star_num, na.rm = TRUE)
      scaling_num <- 1
      plus_sums <- colSums(xi * e1_w_plus) * scaling_num
      minus_sums <- colSums(xi * e1_w_minus) * scaling_num
      
      if (need_separate_components) {
        nu_plus <- plus_sums
        nu_minus <- minus_sums
      }
      
      out_sharp <- plus_sums - minus_sums
    }
    
    if (!fuzzy) {
      # For sharp design, return the result and optionally the components
      if (need_separate_components) {
        return(list(
          tau = out_sharp,
          nu_plus = nu_plus,
          nu_minus = nu_minus
        ))
      } else {
        return(out_sharp)
      }
    } else {
      # Fuzzy RDD: denominator uses h_star_den
      scaling_den <- 1
      plus_sums2 <- sum(xi * e2_w_plus[, 1]) * scaling_den
      minus_sums2 <- sum(xi * e2_w_minus[, 1]) * scaling_den
      denom_term <- plus_sums2 - minus_sums2
      
      # For fuzzy design, adjust the calculation
      top <- denomT * out_sharp - num_diff * denom_term
      tau_fuzzy <- top / (denomT^2)
      
      # Return the result and optionally the components
      if (need_separate_components) {
        return(list(
          tau = tau_fuzzy,
          nu_plus = nu_plus,
          nu_minus = nu_minus,
          denom_term = denom_term
        ))
      } else {
        return(tau_fuzzy)
      }
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
  
  # Process the bootstrap results
  if (need_separate_components) {
    # Extract tau, nu_plus, and nu_minus from each bootstrap draw
    boot_tau_list <- lapply(boot_list, function(x) x$tau)
    boot_mat <- do.call(cbind, boot_tau_list)
    
    # Store nu_plus and nu_minus for Gini test
    nu_plus_mat <- do.call(cbind, lapply(boot_list, function(x) x$nu_plus))
    nu_minus_mat <- do.call(cbind, lapply(boot_list, function(x) x$nu_minus))
  } else {
    # For regular bootstrap without Gini test
    boot_mat <- do.call(cbind, boot_list)
  }
  
  # Calculate uniform confidence bands
  supvals <- apply(boot_mat, 2, function(colb) max(abs(colb), na.rm = TRUE))
  supvals_sorted <- sort(supvals)
  k <- ceiling((1 - alpha) * (B + 1))
  cval <- supvals_sorted[k]  
  cb_lower <- tauhat -  cval
  cb_upper <- tauhat + cval
  
  # Initialize test_results list
  test_results <- list()
  
  # Process test ranges
  if (is.null(test_ranges)) {
    # If no ranges provided, use the full q_grid
    test_ranges <- list(range(q_grid))
  }
  
  # Expand ranges with multiple points into pairs
  processed_ranges <- list()
  for (range in test_ranges) {
    if (length(range) == 2) {
      processed_ranges[[length(processed_ranges) + 1]] <- range
    } else if (length(range) > 2) {
      # Convert a vector like c(0.25, 0.5, 0.75) into pairs [(0.25, 0.5), (0.5, 0.75)]
      for (i in 1:(length(range) - 1)) {
        processed_ranges[[length(processed_ranges) + 1]] <- range[i:(i+1)]
      }
    }
  }
  
  # For each test type and each range, compute test statistics
  if (length(test) > 0 && !all(test == "none")) {
    for (test_type in test) {
      if (test_type == "none") next
      
      range_results <- list()
      
      # Special case for gini test
      if (test_type == "gini") {
        # Calculate Gini coefficients from the estimated conditional quantile functions
        gini_above <- calculate_gini_from_quantile(q_grid, int_plus)
        gini_below <- calculate_gini_from_quantile(q_grid, int_minus)
        gini_diff <- gini_above - gini_below
        
        # Bootstrap the Gini difference using the stored nu_plus and nu_minus
        gini_diffs_boot <- numeric(B)
        
        for (b in 1:B) {
          # Bootstrap the quantile functions using the nu components
          boot_plus <- int_plus + nu_plus_mat[, b]
          boot_minus <- int_minus + nu_minus_mat[, b]
          
          # Calculate Gini coefficients for this bootstrap sample
          boot_gini_above <- calculate_gini_from_quantile(q_grid, boot_plus)
          boot_gini_below <- calculate_gini_from_quantile(q_grid, boot_minus)
          gini_diffs_boot[b] <- boot_gini_above - boot_gini_below
        }
        
        # Calculate critical value and p-value correctly
        test_stat <- abs(gini_diff)
        test_crit <- stats::quantile(gini_diffs_boot, 1 - alpha, na.rm = TRUE)
        
        
        # We don't need to modify the bootstrap distribution
        # since we're testing if the difference is 0
        p_val <- mean(abs(gini_diffs_boot) >= test_stat, na.rm = TRUE)
        
        range_results[["full_sample"]] <- list(
          gini_above = gini_above,
          gini_below = gini_below,
          gini_diff = gini_diff,
          test_stat = test_stat,
          test_crit_val = test_crit,
          p_value = p_val,
          bootstrap_diffs = gini_diffs_boot,
          description = "Gini Test: H0: No difference in Gini coefficients between treatment and control"
        )
      } else {
        # Standard tests (nullity and homogeneity) with ranges
        for (i in seq_along(processed_ranges)) {
          range <- processed_ranges[[i]]
          # Find indices in q_grid that correspond to the current range
          range_idx <- which(q_grid >= min(range) & q_grid <= max(range))
          
          if (length(range_idx) == 0) {
            warning(sprintf("No quantiles in q_grid fall within range [%f, %f]", min(range), max(range)))
            next
          }
          
          range_str <- sprintf("[%.2f, %.2f]", min(range), max(range))
          
          if (test_type == "nullity") {
            # Test for nullity: H0: tau(q) = 0 for all q in the range
            test_stat <- max(abs(tauhat[range_idx]), na.rm = TRUE)
            supvals_null <- apply(boot_mat[range_idx, , drop = FALSE], 2, 
                                  function(colb) max(abs(colb), na.rm = TRUE))
            test_crit <- stats::quantile(supvals_null, 1 - alpha, na.rm = TRUE)
            p_val <- mean(supvals_null >= test_stat, na.rm = TRUE)
            
            range_results[[range_str]] <- list(
              range = range,
              test_stat = test_stat,
              test_crit_val = test_crit,
              p_value = p_val,
              description = "Nullity Test: H0: tau(q) = 0 for all q in the range"
            )
            
          } else if (test_type == "homogeneity") {
            # Test for homogeneity: H0: tau(q) is constant for all q in the range
            range_tau <- tauhat[range_idx]
            mbar <- mean(range_tau, na.rm = TRUE)
            test_stat <- max(abs(range_tau - mbar), na.rm = TRUE)
            
            # Calculate centered bootstrap samples for this range
            range_boot <- boot_mat[range_idx, , drop = FALSE]
            boot_means <- colMeans(range_boot, na.rm = TRUE)
            boot_centered <- sweep(range_boot, 2, boot_means)
            
            supvals_homo <- apply(boot_centered, 2, 
                                  function(col) max(abs(col), na.rm = TRUE))
            test_crit <- stats::quantile(supvals_homo, 1 - alpha, na.rm = TRUE)
            p_val <- mean(supvals_homo >= test_stat, na.rm = TRUE)
            
            range_results[[range_str]] <- list(
              range = range,
              test_stat = test_stat,
              test_crit_val = test_crit,
              p_value = p_val,
              description = "Homogeneity Test: H0: tau(q) is constant for all q in the range"
            )
          }
        }
      }
      
      test_results[[test_type]] <- range_results
    }
  }
  
  # Return results
  result <- list(
    cb_lower = cb_lower, 
    cb_upper = cb_upper,
    boot_taus = boot_mat,
    supvals = supvals,
    crit_val = cval
  )
  
  # If we performed tests, add them to the results
  if (length(test_results) > 0) {
    result$test_results <- test_results
  }
  
  # Add attributes for backward compatibility
  if (length(test) == 1) {
    attr(result, "test_type") <- test
    if (test != "none" && length(test_results) > 0) {
      # For backward compatibility, add the first test results to the main list
      first_test <- test_results[[test]]
      if (length(first_test) > 0) {
        first_range <- first_test[[1]]
        result$test_stat <- first_range$test_stat
        result$test_crit_val <- first_range$test_crit_val
        result$p_value <- first_range$p_value
      }
    }
  }
  
  return(result)
}