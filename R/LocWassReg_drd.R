#' @title Local Wasserstein Regression with Optional Density Preprocessing
#' @noRd
#' @description Modified Local Fréchet regression to estimate from the left and right of a threshold in a regression discontinuity design, with optional density preprocessing.
#'
#' @param xin A n by p matrix holding the n observations of the predictor.
#' @param qin An n by m matrix with values of quantile functions, or a list of quantile or density functions if `compute_density` is TRUE.
#' @param xout A k by p matrix holding the k output predictor values.
#' @param optns A list of control parameters, including bandwidths, kernel type, support bounds, and threshold `c` for D-RDD.
#' @param compute_density Boolean; if TRUE, converts density inputs to quantiles. Default is FALSE, so qin is expected to be quantiles.
#' @details Available control options:
#' \describe{
#' \item{bw}{Bandwidths for each predictor dimension.}
#' \item{ker}{Kernel function type.}
#' \item{lower}{Lower bound of the distribution's support.}
#' \item{upper}{Upper bound of the distribution's support.}
#' \item{c}{Threshold for D-RDD to apply left and right estimation.}
#' \item{qSup}{Support for quantiles, required if converting densities to quantiles.}
#' \item{denLowerThreshold}{Threshold for density preprocessing (if compute_density is TRUE).}
#' }
#' @importFrom osqp solve_osqp osqpSettings

LocWassRegDRDD = function(xin, qin, xout, optns = list(), compute_density = FALSE) {
  
  # Initial checks for matrix/vector input and dimension consistency
  if(!is.matrix(xin) & !is.vector(xin)){
    stop('xin must be a matrix or vector')
  }
  if(is.vector(xin)){
    xin <- matrix(xin, length(xin))
  }
  if(is.null(xout)){
    xout <- xin
  }
  if(!is.matrix(xout) & !is.vector(xout)){
    stop('xout must be a matrix or vector')
  }
  if(is.vector(xout)){
    xout <- matrix(xout, length(xout))
  }
  if(ncol(xin) != ncol(xout)){
    stop('xin and xout must have the same number of columns')
  }
  
  # Validate bandwidth
  if(is.null(optns$bw)){
    stop("optns$bw has no default values and must be input by user.")
  }
  if(!is.numeric(optns$bw) | (length(optns$bw) != ncol(xin))){
    stop("optns$bw should be a numerical vector of length p.")
  }
  
  # Set kernel type, defaulting to Gaussian
  if(is.null(optns$ker)){
    optns$ker <- 'gauss'
  }
  ker <- kerFctn(optns$ker)
  
  # Define threshold
  threshold <- optns$c
  if (is.null(threshold)) {
    stop("Threshold 'c' must be specified in options for D-RDD.")
  }
  
  # Optional: convert density functions to quantiles if compute_density is TRUE
  if (compute_density) {
    # Validate support for quantiles
    if (is.null(optns$qSup)) {
      stop("optns$qSup must be specified for density to quantile conversion.")
    }
    if (is.list(qin)) {
      # Convert densities to quantiles
      den <- lapply(qin, function(deni) fdadensity::dens2quantile(dens = deni$y, dSup = deni$x, qSup = optns$qSup))
      qin <- t(sapply(den, identity))
    } else if (is.matrix(qin)) {
      # Threshold densities if required
      if (!is.null(optns$denLowerThreshold)) {
        den <- apply(qin, 1, function(d) {
          lower <- optns$denLowerThreshold / diff(range(d$x))
          if (sum(d$y < lower) > 0) {
            d$y[d$y < lower] <- lower
            d$y <- d$y / pracma::trapz(d$x, d$y)
          }
          fdadensity::dens2quantile(dens = d$y, dSup = d$x, qSup = optns$qSup)
        })
        qin <- t(den)
      }
    }
  }
  
  # LFR weights with indicators for left and right estimation
  getLFRweights = function(x0, side){
    indicator = if (side == "left") as.numeric(xin[, 1] < threshold) else as.numeric(xin[, 1] >= threshold)
    aux = indicator * K(xin - matrix(t(x0), nrow = n, ncol = length(x0), byrow = TRUE), optns$bw)
    mu0 = mean(aux)
    mu1 = colMeans(aux * (xin - matrix(t(x0), nrow = n, ncol = length(x0), byrow = TRUE)))
    mu2 = 0
    for(i in 1:n){
      mu2 = mu2 + aux[i] * (xin[i, ] - x0) %*% t(xin[i, ] - x0) / n
    }
    sL = array(0, n)
    for(i in 1:n){
      sL[i] = aux[i] * (1 - t(mu1) %*% solve(mu2) %*% (xin[i, ] - x0))
    }
    s = sum(sL)
    return(sL / s)  # Normalized weights
  }
  
  # Set up constraints for quantile bounds
  A <- cbind(diag(m), rep(0, m)) + cbind(rep(0, m), -diag(m))
  if (!is.null(optns$upper) & !is.null(optns$lower)) {
    b0 <- c(optns$lower, rep(0, m - 1), -optns$upper)
  } else if (!is.null(optns$upper)) {
    A <- A[, -1]
    b0 <- c(rep(0, m - 1), -optns$upper)
  } else if (!is.null(optns$lower)) {
    A <- A[, -ncol(A)]
    b0 <- c(optns$lower, rep(0, m - 1))
  } else {
    A <- A[, -c(1, ncol(A))]
    b0 <- rep(0, m - 1)
  }
  
  # Sparse matrices for osqp solver
  Pmat <- as(diag(m), "sparseMatrix")
  Amat <- as(t(A), "sparseMatrix")
  
  # Compute quantiles and return weights for left and right estimations
  results <- lapply(c("left", "right"), function(side) {
    qout <- sapply(1:k, function(j){
      s = getLFRweights(xout[j, ], side)
      s = as.vector(s)
      gx <- colMeans(qin * s) * n
      res <- do.call(osqp::solve_osqp, list(P = Pmat, q = -gx, A = Amat, l = b0, pars = osqp::osqpSettings(verbose = FALSE)))
      return(sort(res$x))
    })
    list(qout = t(qout), weights = sapply(1:k, function(j) getLFRweights(xout[j, ], side)))
  })
  
  return(results)
}
