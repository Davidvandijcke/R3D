# File: R3D/R/r3d_utils.R

#' @title Internal Utilities for R3D
#' @description
#'  - .compute_empirical_qmat: convert Y_list into a matrix of quantiles
#'  - ...
#' @keywords internal

.compute_empirical_qmat <- function(Y_list, q_grid) {
  lapply(Y_list, function(y) stats::quantile(y, probs=q_grid, type=7, na.rm=TRUE))
}


#' @title mclapply.hack
#' @description This function mimics forking (mclapply) also for Windows.
#' @param ... arguments as for lapply
#' @param verbose bool
#' @param mc.cores integer number of cores
#' @keywords internal
mclapply.hack <- function(..., verbose=FALSE, mc.cores=1) {
  if (mc.cores == 1) {
    return(lapply(...))
  }
  
  if (Sys.info()[['sysname']] == 'Windows') {
    # Windows => parLapply
    if (is.null(mc.cores)) {
      size.of.list <- length(list(...)[[1]])
      mc.cores <- min(size.of.list, parallel::detectCores())
    }
    cl <- parallel::makeCluster(mc.cores, outfile="")
    on.exit(parallel::stopCluster(cl), add=TRUE)
    
    # Copy environment
    this.env <- environment()
    while (!identical(this.env, globalenv())) {
      parallel::clusterExport(cl, ls(all.names=TRUE, envir=this.env), envir=this.env)
      this.env <- parent.env(this.env)
    }
    parallel::clusterExport(cl, ls(all.names=TRUE, envir=globalenv()), envir=globalenv())
    
    # Load packages on clusters if needed ...
    # run parLapply
    out <- parallel::parLapply(cl, ...)
    
    if (verbose) {
      message("mclapply.hack: running on Windows with parLapply")
    }
    out
  } else {
    # Non-windows => real mclapply
    parallel::mclapply(..., mc.cores=mc.cores)
  }
}
