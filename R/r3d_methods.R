#' Summarize an r3d Object
#'
#' Displays a summary of the local polynomial (or Fréchet) regression discontinuity model
#' fitted by \code{\link{r3d}}. Includes basic information about the method, polynomial order,
#' fuzzy vs. sharp design, sample size, bandwidth choice(s), and a numeric summary of
#' the estimated distributional treatment effect \eqn{\tau(q)}.
#' If bootstrap inference was performed, also reports uniform confidence bands
#' and test results.
#'
#' @param object An \code{r3d} object produced by \code{\link{r3d}}.
#' @param samples A numeric vector of quantile cut-points at which to display aggregated
#'   effects. Defaults to \code{c(0.25,0.5,0.75)}. Used only for convenient partial summaries.
#' @param ... Additional arguments (not used).
#'
#' @details
#' By default, prints:
#' \itemize{
#'   \item The call used to create the \code{r3d} object,
#'   \item Basic setup and bandwidth details,
#'   \item A table of \eqn{\tau(q)} for all quantiles in \code{object$q_grid},
#'   \item If available, uniform confidence bands and test statistics from the bootstrap,
#'   \item Optionally, average effects on sub-ranges of \eqn{q}.
#' }
#'
#' @return The \code{object} is returned invisibly. This function is primarily for printing.
#'
#' @seealso \code{\link{r3d}}, \code{\link{plot.r3d}}, \code{\link{print.r3d}}
#'
#' @examples
#' \dontrun{
#'   fit <- r3d(X, Y_list, boot=TRUE)
#'   summary(fit, samples=c(0.25, 0.75))
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
        round(min(unlist(object$bandwidths), na.rm=TRUE), 4), "to", 
        round(max(unlist(object$bandwidths), na.rm=TRUE), 4), "\n")
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

#' Plot an r3d Object
#'
#' Produces a simple plot of the estimated distributional RD effect \eqn{\tau(q)} as a function
#' of the quantile \eqn{q}, optionally with uniform confidence bands if they are available in
#' the fitted \code{r3d} object.
#'
#' @param x An \code{r3d} object from \code{\link{r3d}}.
#' @param main An overall title for the plot. If \code{NULL}, a default title is constructed.
#' @param xlab Label of the x-axis, defaults to \code{"Quantile"}.
#' @param ylab Label of the y-axis, defaults to \code{"Treatment Effect"}.
#' @param col Color for the main line. Defaults to \code{"blue"}.
#' @param lwd Line width for the main curve. Defaults to 1.5.
#' @param ci_col Color for the confidence-band boundaries. Defaults to \code{"gray"}.
#' @param ci_lty Line type for the confidence-band boundaries. Defaults to 2 (dashed).
#' @param ref_line Logical indicating whether to draw a horizontal reference line at 0. Default \code{TRUE}.
#' @param ... Additional plotting parameters passed to \code{\link[graphics]{plot}}.
#'
#' @details
#' If bootstrap inference was performed (\code{boot=TRUE} in \code{\link{r3d}}), the object
#' contains \code{cb_lower} and \code{cb_upper} for uniform confidence bands.
#' These bands are added to the plot if available.
#'
#' @return Returns the \code{x} object invisibly.
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
plot.r3d <- function(x, main=NULL, ylim=NULL, xlab="Quantile", ylab="Treatment Effect", 
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
  if(!is.null(obj$boot_out)) {
    cb_l <- obj$boot_out$cb_lower
    cb_u <- obj$boot_out$cb_upper
  }
  
  # Create plot with sensible y-axis limits
  if (is.null(ylim)) { 
    if (!is.null(obj$boot_out)) {
      ylim <- range(c(tau, cb_l, cb_u), na.rm=TRUE)
    } else {
      ylim <- range(tau, na.rm=TRUE)
    }
    
    # Add small margin to y-limits
    y_margin <- 0.1 * diff(ylim)
    ylim <- ylim + c(-y_margin, y_margin)
  }

  
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

#' Print Method for r3d Objects
#'
#' Gives a concise overview of an \code{r3d} object's main properties, including the design type
#' (sharp or fuzzy), local polynomial order, sample size, and bandwidth choice. It also shows
#' a numeric summary (min, median, max) of the estimated distributional RD effect \eqn{\tau(q)}.
#'
#' @param x An \code{r3d} object returned by \code{\link{r3d}}.
#' @param ... Additional arguments (not used).
#'
#' @details
#' This function is invoked automatically when an \code{r3d} object is printed on the console,
#' e.g., simply by typing its name. For a more detailed summary, use \code{\link{summary.r3d}}.
#'
#' @return Returns the \code{x} object invisibly.
#'
#' @examples
#' \dontrun{
#'   fit <- r3d(X, Y_list, boot=TRUE)
#'   print(fit)
#' }
#'
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