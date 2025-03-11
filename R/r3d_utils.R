# File: R3D/R/r3d_utils.R

#' Internal function to compute empirical quantile matrix
#' 
#' @param Y_list List of numeric vectors representing distributions
#' @param q_grid Vector of quantiles to compute
#' @return Matrix with each row containing quantiles for one unit
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
      Qmat[i,] <- stats::quantile(y_i, probs=q_grid, na.rm=TRUE)
    }
  }
  
  return(Qmat)
}


# File: R3D/R/utils.R


#' @title Internal: Cross-Platform Parallel mclapply Hack
#'
#' @description
#' A wrapper around \code{\link[parallel]{mclapply}} that uses \code{parLapply} on Windows. 
#' If \code{mc.cores=1} or the parallel package is unavailable, 
#' it falls back to a simple \code{\link{lapply}} call.
#'
#' @param X A vector (list) to iterate over.
#' @param FUN The function to apply to each element of \code{X}.
#' @param mc.cores Integer, number of cores to use.
#' @param ... Additional arguments to FUN.
#'
#' @details
#' Internal convenience function to unify parallel calls across different OS. 
#' 
#' @return A list of the same length as \code{X}, each element the result of \code{FUN}.
#' @keywords internal
#' @noRd
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

#' @title Internal: Dot Product Helper
#'
#' @description 
#' Computes the dot product of two vectors efficiently, ignoring \code{NA} values.
#'
#' @param x A numeric vector.
#' @param y A numeric vector of the same length as \code{x}.
#'
#' @return A scalar, the sum of elementwise products of \code{x} and \code{y}.
#'
#' @keywords internal
#' @noRd
.dot_product <- function(x, y) {
  sum(x * y, na.rm = TRUE)
}