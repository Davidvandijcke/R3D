#' Main Estimator for Distributional RD (R3D)
#'
#' Fits a regression discontinuity design with \emph{distributional} outcomes (random distributions)
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
#' @param boot Logical. If \code{TRUE}, automatically runs \code{\link{r3d_bootstrap}}
#'   after estimation to produce uniform confidence bands and tests.
#' @param boot_reps Integer, number of bootstrap draws if \code{boot=TRUE}.
#' @param boot_cores Number of CPU cores for parallelizing the bootstrap. Default 1.
#' @param alpha Significance level for the uniform confidence bands. Default 0.05.
#' @param test Either \code{"none"}, \code{"nullity"}, or \code{"homogeneity"} (see \code{\link{r3d_bootstrap}}).
#' @param ... Additional arguments passed to the bandwidth selection function \code{\link{r3d_bwselect}}.
#'
#' @details
#' This is the main user-facing function for the R3D approach. Internally, it:
#' \enumerate{
#'   \item Selects bandwidth(s) via \code{\link{r3d_bwselect}},
#'   \item Computes local polynomial fits on each side of the cutoff for each quantile in \code{q_grid},
#'   \item For \code{method="frechet"}, performs a projection of the fitted quantile curves onto the
#'         space of valid quantile functions (via isotonic regression),
#'   \item (Optionally) runs the multiplier bootstrap for uniform inference.
#' }
#'
#' @return An S3 object of class \code{"r3d"}, containing (among others):
#' \describe{
#'   \item{\code{tau}}{Estimated distributional treatment effect \eqn{\tau(q)} for each \eqn{q} in \code{q_grid}.}
#'   \item{\code{q_grid}}{The quantile grid used.}
#'   \item{\code{bandwidths}}{The MSE- or IMSE-optimal bandwidth(s).}
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
#'   fit <- r3d(X, Y_list, T=T, cutoff=0,
#'              method="frechet", p=2, fuzzy=TRUE,
#'              boot=TRUE, boot_reps=200, alpha=0.05, test="nullity")
#'
#'   # Inspect results
#'   print(fit)
#'   summary(fit)
#'   plot(fit)
#' }
#'
#' @useDynLib R3D, .registration = TRUE
#' @export
r3d <- function(X, Y_list, T = NULL,
                cutoff = 0,
                method = c("simple", "frechet"),
                p = 2,
                q_grid = seq(0.01, 0.99, 0.01),
                fuzzy = FALSE,
                kernel_fun = c("epanechnikov", "triangular", "uniform"),
                s = 1,
                boot = FALSE,
                boot_reps = 200,
                boot_cores = 1,
                alpha = 0.05,
                test = c("none", "nullity", "homogeneity"),
                ...)
{
  method <- match.arg(method)
  kernel_fun <- match.arg(kernel_fun)
  test <- match.arg(test)
  
  kernel <- switch(kernel_fun,
                   triangular = function(u) pmax(0, 1 - abs(u)),
                   epanechnikov = function(u) 0.75 * pmax(0, 1 - u^2),
                   uniform = function(u) 0.5 * (abs(u) <= 1),
                   stop("Unknown kernel type."))
  
  # Validate inputs
  validate_r3d_inputs(X, Y_list, T, cutoff, method, as.integer(p), q_grid, fuzzy, 
                      kernel_fun, as.integer(s), boot, as.integer(boot_reps), as.integer(boot_cores), alpha, test)
  
  n <- length(X)
  
  nQ <- length(q_grid)
  
  # 1) Bandwidth selection -- MODIFIED
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
                        ...)
  
  # Extract bandwidths -- MODIFIED
  h_star_num <- bwres$h_star_num
  h_star_den <- if (fuzzy) bwres$h_star_den else NULL
  
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
      cores = boot_cores
    )
    out$bootstrap$boot_out <- boot_out
    out$boot_out <- boot_out
  }
  
  out
}#' Summarize an r3d Object
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
}#' Bandwidth Selection for Distributional RDD
#'
#' Computes bandwidth(s) for local polynomial estimation in a regression discontinuity
#' setting with distributional (Fréchet-valued) outcomes. Implements a three-step pilot
#' procedure to find either MSE-optimal (per-quantile) bandwidths or a single IMSE-optimal
#' bandwidth, depending on \code{method}.
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
#' @param ... Additional arguments for future expansions.
#'
#' @details
#' Implements a three-step procedure:
#' \enumerate{
#'   \item Estimates \( f_X(0) \) using Silverman’s rule and computes pilot bandwidths via global polynomials.
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
#'   \item{\code{f_X_hat}}{Estimated density of \( X \) at cutoff.}
#' }
#'
#' @references
#' Van Dijcke, D. (2025). \emph{Regression Discontinuity Design with Distributional Outcomes (R3D)}.
#' Working paper.
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
                         ...)
{
  method <- match.arg(method)
  n <- length(X)
  nQ <- length(q_grid)
  
  if (fuzzy && is.null(T)) stop("Fuzzy design requires treatment variable T.")
  
  # Re-center running variable
  Xc <- X - cutoff
  if (all(Xc <= 0) || all(Xc >= 0)) stop("Invalid X data: no data on both sides of cutoff.")
  
  # Step 1: Preliminary Bandwidth Calculation
  # Estimate f_X(0) using Silverman’s rule
  sigma_X <- stats::sd(Xc)
  c_n <- 1.06 * sigma_X * n^(-1/5)
  f_X_hat <- mean(kernel(Xc / c_n)) / c_n
  
  # Compute empirical quantile matrix for Y
  Qmat <- .compute_empirical_qmat(Y_list, q_grid)
  
  # Global polynomial fits for pilot bandwidths
  fit_global_poly <- function(Y, X, order) {
    fit <- lm(Y ~ poly(X, degree = order, raw = TRUE))
    coefs <- coef(fit)
    derivs <- coefs[2:(order + 1)] * factorial(0:(order - 1))
    resid_var <- var(resid(fit))
    list(derivs = derivs, resid_var = resid_var)
  }
  
  # Numerator (outcome)
  B_plus_num <- numeric(nQ)
  B_minus_num <- numeric(nQ)
  V_plus_num <- numeric(nQ)
  V_minus_num <- numeric(nQ)
  
  for (qi in seq_len(nQ)) {
    idx_plus <- Xc >= 0
    idx_minus <- Xc < 0
    if (sum(idx_plus) > (s + 1)) {
      fit_plus <- fit_global_poly(Qmat[idx_plus, qi], Xc[idx_plus], s + 1)
      B_plus_num[qi] <- fit_plus$derivs[s + 1] / factorial(s + 1)
      V_plus_num[qi] <- fit_plus$resid_var
    }
    if (sum(idx_minus) > (s + 1)) {
      fit_minus <- fit_global_poly(Qmat[idx_minus, qi], Xc[idx_minus], s + 1)
      B_minus_num[qi] <- fit_minus$derivs[s + 1] / factorial(s + 1)
      V_minus_num[qi] <- fit_minus$resid_var
    }
  }
  
  # Denominator (treatment) in fuzzy case
  if (fuzzy) {
    fit_T_plus <- fit_global_poly(T[Xc >= 0], Xc[Xc >= 0], s + 1)
    fit_T_minus <- fit_global_poly(T[Xc < 0], Xc[Xc < 0], s + 1)
    B_plus_den <- fit_T_plus$derivs[s + 1] / factorial(s + 1)
    B_minus_den <- fit_T_minus$derivs[s + 1] / factorial(s + 1)
    V_plus_den <- fit_T_plus$resid_var
    V_minus_den <- fit_T_minus$resid_var
  }
  
  # Compute pilot bandwidths
  if (method == "simple") {
    pilot_h_num <- numeric(nQ)
    for (qi in seq_len(nQ)) {
      bias_diff <- B_plus_num[qi] - B_minus_num[qi]
      var_sum <- V_plus_num[qi] + V_minus_num[qi]
      if (abs(bias_diff) > 1e-14) {
        ratio <- (var_sum) / (2 * (s + 1) * bias_diff^2)
        pilot_h_num[qi] <- ratio^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
      } else {
        pilot_h_num[qi] <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
      }
    }
  } else {
    dq <- diff(range(q_grid)) / (length(q_grid) - 1)
    A_s_num <- sum((B_plus_num - B_minus_num)^2) * dq
    B_s_num <- sum(V_plus_num + V_minus_num) * dq
    ratio <- B_s_num / (2 * (s + 1) * max(A_s_num, 1e-14))
    pilot_h_num <- ratio^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
    if (is.na(pilot_h_num) || pilot_h_num <= 0) pilot_h_num <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
  }
  
  if (fuzzy) {
    bias_diff_den <- B_plus_den - B_minus_den
    var_sum_den <- V_plus_den + V_minus_den
    if (abs(bias_diff_den) > 1e-14) {
      ratio_den <- var_sum_den / (2 * (s + 1) * bias_diff_den^2)
      pilot_h_den <- ratio_den^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
    } else {
      pilot_h_den <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
    }
  }
  
  # Step 2: First-Stage Local Polynomial Fits
  B_plus <- numeric(nQ)
  B_minus <- numeric(nQ)
  V_plus <- numeric(nQ)
  V_minus <- numeric(nQ)
  
  for (qi in seq_len(nQ)) {
    h_pilot <- if (method == "simple") pilot_h_num[qi] else pilot_h_num
    idx_plus <- (Xc >= 0) & (abs(Xc) <= h_pilot)
    if (sum(idx_plus) > (s + 1)) {
      Xp <- Xc[idx_plus]
      Yp <- Qmat[idx_plus, qi]
      Wp <- kernel(Xp / h_pilot)
      X_mat_p <- outer(Xp / h_pilot, 0:s, "^")
      W_diag_p <- diag(Wp)
      XWX_p <- t(X_mat_p) %*% W_diag_p %*% X_mat_p
      XWY_p <- t(X_mat_p) %*% W_diag_p %*% Yp
      if (rcond(XWX_p) > 1e-10) {
        coef_p <- solve(XWX_p, XWY_p)
        B_plus[qi] <- coef_p[s + 1] / factorial(s + 1)
        fitted_p <- X_mat_p %*% coef_p
        V_plus[qi] <- sum(Wp * (Yp - fitted_p)^2) / sum(Wp)
      }
    }
    idx_minus <- (Xc < 0) & (abs(Xc) <= h_pilot)
    if (sum(idx_minus) > (s + 1)) {
      Xm <- Xc[idx_minus]
      Ym <- Qmat[idx_minus, qi]
      Wm <- kernel(Xm / h_pilot)
      X_mat_m <- outer(Xm / h_pilot, 0:s, "^")
      W_diag_m <- diag(Wm)
      XWX_m <- t(X_mat_m) %*% W_diag_m %*% X_mat_m
      XWY_m <- t(X_mat_m) %*% W_diag_m %*% Ym
      if (rcond(XWX_m) > 1e-10) {
        coef_m <- solve(XWX_m, XWY_m)
        B_minus[qi] <- coef_m[s + 1] / factorial(s + 1)
        fitted_m <- X_mat_m %*% coef_m
        V_minus[qi] <- sum(Wm * (Ym - fitted_m)^2) / sum(Wm)
      }
    }
  }
  
  # Denominator in fuzzy case
  if (fuzzy) {
    idx_plus <- (Xc >= 0) & (abs(Xc) <= pilot_h_den)
    if (sum(idx_plus) > (s + 1)) {
      Xp <- Xc[idx_plus]
      Tp <- T[idx_plus]
      Wp <- kernel(Xp / pilot_h_den)
      X_mat_p <- outer(Xp / pilot_h_den, 0:s, "^")
      W_diag_p <- diag(Wp)
      XWX_p <- t(X_mat_p) %*% W_diag_p %*% X_mat_p
      XWY_p <- t(X_mat_p) %*% W_diag_p %*% Tp
      if (rcond(XWX_p) > 1e-10) {
        coef_p <- solve(XWX_p, XWY_p)
        B_plus_den <- coef_p[s + 1] / factorial(s + 1)
        fitted_p <- X_mat_p %*% coef_p
        V_plus_den <- sum(Wp * (Tp - fitted_p)^2) / sum(Wp)
      }
    }
    idx_minus <- (Xc < 0) & (abs(Xc) <= pilot_h_den)
    if (sum(idx_minus) > (s + 1)) {
      Xm <- Xc[idx_minus]
      Tm <- T[idx_minus]
      Wm <- kernel(Xm / pilot_h_den)
      X_mat_m <- outer(Xm / pilot_h_den, 0:s, "^")
      W_diag_m <- diag(Wm)
      XWX_m <- t(X_mat_m) %*% W_diag_m %*% X_mat_m
      XWY_m <- t(X_mat_m) %*% W_diag_m %*% Tm
      if (rcond(XWX_m) > 1e-10) {
        coef_m <- solve(XWX_m, XWY_m)
        B_minus_den <- coef_m[s + 1] / factorial(s + 1)
        fitted_m <- X_mat_m %*% coef_m
        V_minus_den <- sum(Wm * (Tm - fitted_m)^2) / sum(Wm)
      }
    }
  }
  
  # Step 3: Optimal Bandwidth
  if (method == "simple") {
    h_star_num <- numeric(nQ)
    for (qi in seq_len(nQ)) {
      bias_diff <- B_plus[qi] - B_minus[qi]
      var_sum <- (V_plus[qi] + V_minus[qi]) / f_X_hat
      if (abs(bias_diff) > 1e-14) {
        ratio <- (var_sum) / (2 * (s + 1) * bias_diff^2)
        h_star_num[qi] <- ratio^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
      } else {
        h_star_num[qi] <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
      }
    }
  } else {
    dq <- diff(range(q_grid)) / (length(q_grid) - 1)
    A_s <- sum((B_plus - B_minus)^2) * dq
    B_s <- sum(V_plus + V_minus) * dq / f_X_hat
    ratio <- B_s / (2 * (s + 1) * max(A_s, 1e-14))
    h_star_num <- ratio^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
    if (is.na(h_star_num) || h_star_num <= 0) h_star_num <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
  }
  
  if (fuzzy) {
    bias_diff_den <- B_plus_den - B_minus_den
    var_sum_den <- (V_plus_den + V_minus_den) / f_X_hat
    if (abs(bias_diff_den) > 1e-14) {
      ratio_den <- var_sum_den / (2 * (s + 1) * bias_diff_den^2)
      h_star_den <- ratio_den^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
    } else {
      h_star_den <- 1.06 * sigma_X * n^(-1 / (2 * s + 1))
    }
  } else {
    h_star_den <- NULL
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
}#' @title R3D: Regression Discontinuity with Distributional Outcomes
#' @name R3D-package
#' @docType package
#'
#' @description
#' The \pkg{R3D} package provides methods to estimate, infer, and visualize 
#' a new class of regression discontinuity designs where the outcome is 
#' a distribution rather than a single scalar. It includes:
#' \itemize{
#'   \item \code{\link{r3d}}: main estimation function
#'   \item \code{\link{r3d_bwselect}}: bandwidth selection
#'   \item \code{\link{r3d_bootstrap}}: multiplier bootstrap for uniform inference
#'   \item S3 methods: \code{\link{summary.r3d}}, \code{\link{plot.r3d}}, \code{\link{print.r3d}}
#' }
#'
#' @author 
#' Your Name (David Van Dijcke) \email{dvdijcke@umich.edu}
#'
#' @references
#' Van Dijcke, D. (2025). \emph{Regression Discontinuity Design with Distributional Outcomes.}
#' Working paper. 
#'
#' @keywords package
NULL# File: R3D/R/r3d_utils.R

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