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
summary.r3d <- function(object, samples=c(0.25,0.5,0.75), ...) {
  cat("Call:\n")
  print(object$call)
  cat("\nMethod:", object$method, "\n")
  cat("Polynomial order p:", object$p, "\n")
  cat("Fuzzy:", object$fuzzy,"\n")
  cat("Sample size:", length(object$X), "\n")
  
  cat("Bandwidth(s):\n")
  if(object$method=="simple") {
    cat(" Per-quantile MSE. Range:", 
        round(min(object$bandwidths, na.rm=TRUE), 4), "to", 
        round(max(object$bandwidths, na.rm=TRUE), 4), "\n")
  } else {
    cat(" Single IMSE bandwidth:", round(object$results$h_used, 4),"\n")
  }
  
  cat("\nQuantile treatment effects:\n")
  outdf <- data.frame(
    Quantile = object$q_grid, 
    Effect = round(object$tau, 4),
    stringsAsFactors = FALSE
  )
  print(outdf)
  
  # If we have bootstrap results, we can show uniform CB or test results
  if(!is.null(object$boot_out)) {
    cat("\nUniform Confidence Bands (alpha=", 
        format(attr(object$boot_out, "alpha", exact=TRUE) %||% 0.05, digits=2), "):\n", sep="")
    cb_l <- object$boot_out$cb_lower
    cb_u <- object$boot_out$cb_upper
    
    cb_df <- data.frame(
      Quantile = object$q_grid,
      Lower = round(cb_l, 4), 
      Upper = round(cb_u, 4),
      stringsAsFactors = FALSE
    )
    print(cb_df)
    
    # test
    ts <- object$boot_out$test_stat
    if(!is.na(ts)) {
      cat("\nHypothesis test results:\n")
      cat("Test statistic:", round(ts, 4),"\n")
      cat("Critical value:", round(object$boot_out$test_crit_val, 4),"\n")
      cat("P-value:", round(object$boot_out$p_value, 4),"\n")
      
      if(!is.null(attr(object$boot_out, "test_type"))) {
        test_type <- attr(object$boot_out, "test_type")
      } else {
        # Try to infer from object$call
        if(!is.null(object$call$test)) {
          test_type <- as.character(object$call$test)
        } else {
          test_type <- "unknown"
        }
      }
      
      if(test_type == "nullity") {
        cat("Test: H0: tau(q) = 0 for all q\n")
      } else if(test_type == "homogeneity") {
        cat("Test: H0: tau(q) is constant across q\n")
      } else {
        cat("Test type:", test_type, "\n")
      }
    }
    
    # Aggregated effect: partition q into intervals
    s0 <- c(0, samples, 1)
    cat("\nAggregated distributional effects:\n")
    agg_mat <- NULL
    for(i in seq_len(length(s0)-1)){
      from <- s0[i]
      to <- s0[i+1]
      sel <- which(object$q_grid >= from & object$q_grid <= to)
      if(length(sel) > 0) {
        mean_eff <- mean(object$tau[sel])
        ci_l <- mean(object$boot_out$cb_lower[sel])
        ci_u <- mean(object$boot_out$cb_upper[sel])
        
        agg_mat <- rbind(agg_mat, 
                         data.frame(
                           Quantile_Range = paste(from, "-", to),
                           Average_Effect = round(mean_eff, 4),
                           CI_Lower = round(ci_l, 4),
                           CI_Upper = round(ci_u, 4),
                           stringsAsFactors = FALSE
                         ))
      }
    }
    print(agg_mat, row.names=FALSE)
  } else {
    cat("\nNo bootstrap results. Set boot=TRUE in r3d() to get inference.\n")
  }
  invisible(object)
}

# Helper for summary method
`%||%` <- function(x, y) if(is.null(x)) y else x

#' @title Plot Method for r3d Objects
#'
#' @description
#' Plots the estimated distributional RD effects, \eqn{\tau(q)} versus \eqn{q}, 
#' optionally with uniform confidence bands (if available).
#'
#' @param x An \code{r3d} object from \code{\link{r3d}}.
#' @param main Main title for the plot.
#' @param xlab X-axis label (defaults to "Quantile").
#' @param ylab Y-axis label (defaults to "Treatment Effect").
#' @param col Color for the main line (defaults to "blue").
#' @param lwd Line width for the main line (defaults to 1.5).
#' @param ci_col Color for confidence-band lines (defaults to "gray").
#' @param ci_lty Line type for confidence-band lines (defaults to 2).
#' @param ref_line Logical, if \code{TRUE}, draws a horizontal reference line at 0.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}.
#'
#' @details
#' This plots the function \eqn{\tau(q)} stored in \code{x$tau} across the grid \code{x$q_grid}. 
#' If \code{x$boot_out} is non-null and includes \code{cb_lower, cb_upper}, 
#' then those uniform confidence bands are overlaid.
#'
#' @return Invisibly returns \code{x}. Called for its side effects (plot).
#'
#' @seealso \code{\link{r3d}}, \code{\link{summary.r3d}}
#'
#' @examples
#' \dontrun{
#'   fit <- r3d(X, Y_list, boot=TRUE)
#'   plot(fit, main="Distributional RD Effects", ref_line=TRUE)
#' }
#'
#' @export
plot.r3d <- function(x, main=NULL, xlab="Quantile", ylab="Treatment Effect", 
                     col="blue", lwd=1.5, ci_col="gray", ci_lty=2, 
                     ref_line=TRUE, ...) {
  obj <- x
  tau <- obj$tau
  qq <- obj$q_grid
  if(is.null(main)) {
    main <- paste0("Distributional RDD Effects (", 
                   ifelse(obj$fuzzy, "Fuzzy", "Sharp"), ", ", 
                   obj$method, ")")
  }
  
  # Create plot with sensible y-axis limits
  if(!is.null(obj$boot_out)) {
    cb_l <- obj$boot_out$cb_lower
    cb_u <- obj$boot_out$cb_upper
    ylim <- range(c(tau, cb_l, cb_u), na.rm=TRUE)
  } else {
    ylim <- range(tau, na.rm=TRUE)
  }
  
  # Add small margin to y-limits
  y_margin <- 0.1 * diff(ylim)
  ylim <- ylim + c(-y_margin, y_margin)
  
  # Create the plot
  plot(qq, tau, type="l", main=main, xlab=xlab, ylab=ylab, 
       ylim=ylim, col=col, lwd=lwd, ...)
  
  # Add reference line at zero if requested
  if(ref_line) {
    abline(h=0, col="darkgray", lty=3)
  }
  
  # Add confidence bands if available
  if(!is.null(obj$boot_out)) {
    graphics::lines(qq, cb_l, col=ci_col, lty=ci_lty)
    graphics::lines(qq, cb_u, col=ci_col, lty=ci_lty)
    graphics::legend("topleft",
                     c("Point estimate", "Uniform CI"),
                     col=c(col, ci_col), lty=c(1, ci_lty),
                     lwd=c(lwd, 1))
  }
  
  invisible(x)
}

#' @title Print Method for r3d Objects
#'
#' @description
#' Provides a concise printed display for \code{r3d} objects. 
#' Shows the method (sharp/fuzzy, polynomial order, sample size, whether bootstrap was used, etc.), 
#' plus a short numeric summary of the estimated treatment effect \code{tau}.
#'
#' @param x An \code{r3d} object from \code{\link{r3d}}.
#' @param ... Additional parameters (not used).
#'
#' @details
#' If you need a more detailed summary, use \code{\link{summary.r3d}}.
#'
#' @seealso 
#' \code{\link{r3d}}, \code{\link{summary.r3d}}, \code{\link{plot.r3d}}
#'
#' @examples
#' \dontrun{
#'   fit <- r3d(X, Y_list)
#'   print(fit)
#' }
#'
#' @return Invisibly returns \code{x}.
#' @export
print.r3d <- function(x, ...) {
  cat("R3D: Regression Discontinuity with Distributional Outcomes\n")
  cat("-------------------------------------------------------\n")
  cat("Method:", x$method, "\n")
  cat("Type:", ifelse(x$fuzzy, "Fuzzy RDD", "Sharp RDD"), "\n")
  cat("Polynomial order:", x$p, "\n")
  
  # Bandwidth info
  if(x$method == "simple") {
    cat("Bandwidths: MSE-optimal (varies by quantile)\n")
  } else {
    cat("Bandwidth: IMSE-optimal =", format(x$results$h_used, digits=4), "\n")
  }
  
  # Sample sizes
  cat("Sample size:", length(x$X), "\n")
  cat("Quantiles evaluated:", length(x$q_grid), "\n")
  
  # Bootstrap info
  if(!is.null(x$boot_out)) {
    cat("Bootstrap: YES (", ncol(x$boot_out$boot_taus), "replications)\n", sep="")
  } else {
    cat("Bootstrap: NO\n")
  }
  
  # Treatment effect summary
  effect_summary <- summary(x$tau)
  cat("\nTreatment Effect Summary:\n")
  print(effect_summary)
  
  # Hint for more detailed information
  cat("\nUse summary() for more details and plot() for visualization\n")
  
  invisible(x)
}