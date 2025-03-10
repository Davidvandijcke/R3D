# File: R3D/R/r3d_methods.R

#' @title Summary method for r3d objects
#' @description prints basic info, bandwidths, and possibly aggregated distribution effects
#' @param object an r3d object
#' @param samples numeric vector of quantile sub-ranges for aggregated summary. Defaults c(0.25, 0.5, 0.75).
#' @param ... not used
#' @export
summary.r3d <- function(object, samples=c(0.25,0.5,0.75), ...) {
  cat("Call:\n")
  print(object$call)
  cat("\nMethod:", object$method, "\n")
  cat("Polynomial order p:", object$call$p, "\n")
  cat("Fuzzy:", object$fuzzy,"\n")
  cat("Bandwidth(s):\n")
  if(object$method=="simple") {
    cat(" Per-quantile MSE. E.g. first few:\n")
    print(head(object$bw))
  } else {
    cat(" Single IMSE bandwidth:", object$bw,"\n")
  }
  cat("\nQuantile treatment effects:\n")
  outdf<- data.frame(q=object$q_grid, tau=object$tau)
  print(head(outdf,10))
  
  # If we have bootstrap results, we can show uniform CB or test results
  if(!is.null(object$boot_out)) {
    cat("\nUniform Confidence Bands at alpha=", object$call$alpha," => \n")
    cb_l<- object$boot_out$cb_lower
    cb_u<- object$boot_out$cb_upper
    cat(" first few:\n")
    tmp<- data.frame(q= object$q_grid, 
                     lower=cb_l, 
                     upper=cb_u)
    print(head(tmp,10))
    
    # test
    ts<- object$boot_out$test_stat
    if(!is.na(ts)) {
      cat("\nTest statistic:", ts,"\n")
      cat("Crit Value:", object$boot_out$test_crit_val,"\n")
      cat("p-value:", object$boot_out$p_value,"\n")
    }
    
    # aggregated effect: e.g. partition q in (0, samples[1]), (samples[1], samples[2]), ...
    # a small example
    s0<- c(0, samples, 1)
    cat("\nAggregated distributional effects:\n")
    agg_mat<- NULL
    for(i in seq_len(length(s0)-1)){
      from<- s0[i]
      to<- s0[i+1]
      sel<- which(object$q_grid>=from & object$q_grid<=to)
      if(length(sel)>0) {
        mean_eff<- mean(object$tau[sel])
        agg_mat<- rbind(agg_mat, 
                        data.frame(q_from=from, q_to=to, avg_effect=mean_eff))
      }
    }
    print(agg_mat, row.names=FALSE)
  } else {
    cat("\nNo bootstrap results. Set boot=TRUE in r3d() to get inference.\n")
  }
  invisible(object)
}

#' @title Plot method for r3d objects
#' @description Plots \tau(q) vs q with optional confidence bands
#' @param x r3d object
#' @param main Title
#' @param xlab x label
#' @param ylab y label
#' @param ... not used
#' @export
plot.r3d <- function(x, main=NULL, xlab="Quantile", ylab="Treatment Effect", ...) {
  obj<- x
  tau<- obj$tau
  qq<- obj$q_grid
  if(is.null(main)) main<- paste("RDD Distributional Effects -", obj$method)
  plot(qq, tau, type="l", main=main, xlab=xlab, ylab=ylab, ...)
  if(!is.null(obj$boot_out)) {
    cb_l<- obj$boot_out$cb_lower
    cb_u<- obj$boot_out$cb_upper
    graphics::lines(qq, cb_l, col="gray", lty=2)
    graphics::lines(qq, cb_u, col="gray", lty=2)
    graphics::legend("topleft",
                     c("point estimate","uniform CB"),
                     col=c("black","gray"), lty=c(1,2))
  }
}
