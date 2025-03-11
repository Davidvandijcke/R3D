# File: R3D/R/r3d_utils.R

#' Compute Empirical Quantiles for Each Unit
#'
#' \strong{Internal function.} Constructs a matrix of empirical quantiles for a list of distributions.
#' Each row corresponds to one unit's distribution, and each column a quantile in \code{q_grid}.
#'
#' @param Y_list A list of numeric vectors, each giving observed samples from one unit's distribution.
#' @param q_grid A numeric vector of quantiles, e.g. \code{seq(0.01, 0.99, 0.01)}.
#'
#' @return A numeric matrix of dimension \code{length(Y_list)} \eqn{\times} \code{length(q_grid)}.
#'   Row \eqn{i} is \eqn{(\hat{Q}_{Y_i}(q_1), \ldots, \hat{Q}_{Y_i}(q_m))}.
#'
#' @keywords internal
.compute_empirical_qmat <- function(Y_list, q_grid) {
  n <- length(Y_list)
  nQ <- length(q_grid)
  Qmat <- matrix(NA, nrow=n, ncol=nQ)
  
  for(i in seq_len(n)) {
    y_i <- Y_list[[i]]
    if(length(y_i) == 0) {
      warning("Empty Y_list entry at index ", i)
      Qmat[i,] <- NA
    } else {
      Qmat[i,] <- stats::quantile(y_i, probs=q_grid, na.rm=TRUE, names=FALSE)
    }
  }
  
  return(Qmat)
}


# File: R3D/R/utils.R


#' Cross-Platform Parallel Lapply Helper
#'
#' \strong{Internal function.} A wrapper around \code{\link[parallel]{mclapply}} that
#' uses \code{\link[parallel]{parLapply}} on Windows (where \code{mclapply} is not supported),
#' and falls back to \code{\link{lapply}} if \code{cores=1} or the \pkg{parallel} package
#' is unavailable. 
#'
#' @param X A list (or vector) over which to iterate.
#' @param FUN The function to apply to each element of \code{X}.
#' @param mc.cores Number of CPU cores requested. Defaults to 1.
#' @param ... Additional arguments passed to \code{FUN}.
#'
#' @return A list of the same length as \code{X}, containing the results of \code{FUN}.
#'
#' @keywords internal
mclapply.hack <- function(X, FUN, mc.cores = 1, ...) {
  if(mc.cores <= 1 || !requireNamespace("parallel", quietly = TRUE)) {
    # Sequential processing
    return(lapply(X, FUN, ...))
  }
  
  # Check OS and use appropriate parallel method
  if(.Platform$OS.type == "windows") {
    # On Windows: use parLapply with explicit cluster creation
    cl <- parallel::makeCluster(mc.cores)
    on.exit(parallel::stopCluster(cl))
    return(parallel::parLapply(cl, X, FUN, ...))
  } else {
    # On Unix-like systems: use mclapply
    return(parallel::mclapply(X, FUN, mc.cores = mc.cores, ...))
  }
}

#' Dot Product of Two Numeric Vectors
#'
#' \strong{Internal function.} Computes the dot product (sum of elementwise products)
#' of two numeric vectors, ignoring NAs.
#'
#' @param x A numeric vector.
#' @param y A numeric vector of the same length as \code{x}.
#'
#' @return A single numeric value, \code{sum(x * y, na.rm=TRUE)}.
#'
#' @keywords internal
.dot_product <- function(x, y) {
  sum(x * y, na.rm = TRUE)
}



# File: r3d_utils.R

#' Validate Inputs for r3d Function
#'
#' This function checks the validity of inputs provided to the \code{r3d} function.
#' It ensures that all inputs are of the correct type, length, and within expected ranges.
#'
#' @param X Numeric vector of the running variable.
#' @param Y_list List of numeric vectors representing outcome distributions.
#' @param T Optional numeric or logical vector for treatment status in fuzzy designs.
#' @param cutoff Numeric scalar for the treatment threshold.
#' @param method Character string specifying the method ("simple" or "frechet").
#' @param p Integer specifying the polynomial order.
#' @param q_grid Numeric vector of quantiles.
#' @param fuzzy Logical indicating if the design is fuzzy.
#' @param kernel_fun Character string specifying the kernel function.
#' @param s Integer specifying the expansion order for pilot bandwidths.
#' @param boot Logical indicating if bootstrap should be performed.
#' @param boot_reps Integer specifying the number of bootstrap repetitions.
#' @param boot_cores Integer specifying the number of CPU cores for bootstrap.
#' @param alpha Numeric scalar for the significance level.
#' @param test Character string specifying the type of test to perform.
#'
#' @return NULL if all checks pass. Otherwise, stops with an informative error message.
#'
#' @keywords internal
validate_r3d_inputs <- function(X, Y_list, T = NULL, cutoff, method, p, q_grid, fuzzy, kernel_fun, s, boot, boot_reps, boot_cores, alpha, test) {
  
  # Check X
  if (!is.numeric(X)) {
    stop("X must be a numeric vector.")
  }
  if (any(is.na(X))) {
    stop("X contains missing values.")
  }
  if (length(X) < 5) {
    warning("X has fewer than 5 observations. Results may be unreliable.")
  }
  
  # Check Y_list
  if (!is.list(Y_list)) {
    stop("Y_list must be a list.")
  }
  if (length(Y_list) != length(X)) {
    stop(paste("Length of Y_list (", length(Y_list), ") does not match length of X (", length(X), ").", sep = ""))
  }
  for (i in seq_along(Y_list)) {
    if (!is.numeric(Y_list[[i]])) {
      stop(paste("Element", i, "in Y_list is not a numeric vector."))
    }
    if (length(Y_list[[i]]) == 0) {
      stop(paste("Element", i, "in Y_list is an empty vector. This may indicate missing data."))
    }
  }
  
  # Check T (if fuzzy = TRUE)
  if (fuzzy) {
    if (is.null(T)) {
      stop("For fuzzy designs, T must be provided.")
    }
    if (!is.numeric(T) && !is.logical(T)) {
      stop("T must be a numeric or logical vector.")
    }
    if (length(T) != length(X)) {
      stop(paste("Length of T (", length(T), ") does not match length of X (", length(X), ").", sep = ""))
    }
  }
  
  # Check cutoff
  if (!is.numeric(cutoff) || length(cutoff) != 1) {
    stop("cutoff must be a numeric scalar.")
  }
  
  # Check method
  if (!method %in% c("simple", "frechet")) {
    stop("method must be either 'simple' or 'frechet'.")
  }
  
  # Check p
  if (!is.integer(p) || p < 1) {
    stop("p must be a positive integer.")
  }
  
  # Check q_grid
  if (!is.numeric(q_grid)) {
    stop("q_grid must be a numeric vector.")
  }
  if (any(q_grid <= 0 | q_grid >= 1)) {
    stop("All values in q_grid must be between 0 and 1.")
  }
  
  # Check fuzzy
  if (!is.logical(fuzzy) || length(fuzzy) != 1) {
    stop("fuzzy must be a logical scalar.")
  }
  
  # Check kernel_fun
  allowed_kernels <- c("epanechnikov", "triangular", "uniform")
  if (!kernel_fun %in% allowed_kernels) {
    stop(paste("kernel_fun must be one of:", paste(allowed_kernels, collapse = ", ")))
  }
  
  # Check s
  if (!is.integer(s) || s < 1) {
    stop("s must be a positive integer.")
  }
  
  # Check boot
  if (!is.logical(boot) || length(boot) != 1) {
    stop("boot must be a logical scalar.")
  }
  
  # Check boot_reps
  if (boot) {
    if (!is.integer(boot_reps) || boot_reps < 1) {
      stop("boot_reps must be a positive integer when boot = TRUE.")
    }
  }
  
  # Check boot_cores
  if (!is.integer(boot_cores) || boot_cores < 1) {
    stop("boot_cores must be a positive integer.")
  }
  
  # Check alpha
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a numeric scalar between 0 and 1.")
  }
  
  # Check test
  allowed_tests <- c("none", "nullity", "homogeneity")
  if (!test %in% allowed_tests) {
    stop(paste("test must be one of:", paste(allowed_tests, collapse = ", ")))
  }
  
  # If all checks pass, return NULL (no error)
  return(NULL)
}

