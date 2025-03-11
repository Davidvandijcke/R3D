r3d <- function(X, Y_list, T=NULL,
                cutoff=0,
                method=c("simple","frechet"),
                p=2,
                q_grid=seq(0.01,0.99,0.01),
                fuzzy=FALSE,
                kernel_fun = c("epanechnikov", "triangular", "uniform"),
                s=1,
                boot=FALSE,
                boot_reps=200,
                boot_cores=1,
                alpha=0.05,
                test=c("none","nullity","homogeneity"),
                ...)
{
  method <- match.arg(method)
  kernel_fun <- match.arg(kernel_fun)
  test   <- match.arg(test)
  
  kernel <- switch(kernel_fun,
                       triangular = function(u) pmax(0, 1 - abs(u)),
                       epanechnikov = function(u) 0.75 * pmax(0, 1 - u^2),
                       uniform = function(u) 0.5 * (abs(u) <= 1),
                       gaussian = function(u) dnorm(u),
                       stop("Unknown kernel type."))
  
  n <- length(X)
  if(fuzzy && is.null(T)) stop("fuzzy=TRUE requires T.")
  if(n < 5) warning("r3d: Very small sample size? Check if meaningful.")
  
  nQ <- length(q_grid)
  
  # 1) Bandwidth selection
  bwres <- r3d_bwselect(X       = X,
                        Y_list  = Y_list,
                        q_grid  = q_grid,
                        method  = method,
                        s       = s,
                        p       = p,
                        kernel  = kernel,
                        cutoff  = cutoff,
                        ...)
  h_star <- bwres$h_star
  
  # 2) Build matrix of empirical quantiles [n x nQ]
  Qmat <- .compute_empirical_qmat(Y_list, q_grid)
  
  # 3) Single "actual" bandwidth if multiple
  if(length(h_star) > 1 && method=="simple"){
    h_use <- mean(h_star, na.rm=TRUE)
  } else {
    h_use <- ifelse(length(h_star)==1, h_star, mean(h_star, na.rm=TRUE))
  }
  # h_use <- 0.5 
  
  # 4) Re-center
  X_centered <- X - cutoff
  
  # 5) Kernel weights for plus/minus: length n each
  side_plus_w <- kernel(X_centered / h_use)
  side_plus_w[X_centered < 0] <- 0
  
  side_minus_w <- kernel(X_centered / h_use)
  side_minus_w[X_centered >= 0] <- 0
  
  # Fortran setup
  N   <- as.integer(n)
  P_  <- as.integer(p)
  NQ_ <- as.integer(nQ)
  
  # 6) Local poly fit on plus side
  outPlus <- .Fortran("locweights",
                      X       = as.double(X_centered),
                      YMAT    = as.double(Qmat),  # shape n x nQ
                      N       = N,
                      P       = P_,
                      H       = as.double(h_use),
                      SIDE    = as.integer(1),
                      KERNELW = as.double(side_plus_w), # length n
                      ALPHA   = double((p+1)*nQ),
                      WINT    = double(n*nQ),
                      INFO    = integer(1),
                      NQ      = NQ_,
                      PACKAGE = "R3D")
  info_plus   <- outPlus$INFO
  alpha_plus  <- matrix(outPlus$ALPHA, nrow=p+1, ncol=nQ)
  w_plus      <- matrix(outPlus$WINT,  nrow=n, ncol=nQ)
  
  # minus side
  outMinus <- .Fortran("locweights",
                       X       = as.double(X_centered),
                       YMAT    = as.double(t(Qmat)),
                       N       = N,
                       P       = P_,
                       H       = as.double(h_use),
                       SIDE    = as.integer(0),
                       KERNELW = as.double(side_minus_w), # length n
                       ALPHA   = double((p+1)*nQ),
                       WINT    = double(n*nQ),
                       INFO    = integer(1),
                       NQ      = NQ_,
                       PACKAGE = "R3D")
  info_minus   <- outMinus$INFO
  alpha_minus  <- matrix(outMinus$ALPHA, nrow=p+1, ncol=nQ)
  w_minus      <- matrix(outMinus$WINT,  nrow=n, ncol=nQ)
  
  if(info_plus != 0)  warning("locweights: plus side singular system? info=", info_plus)
  if(info_minus!=0)  warning("locweights: minus side singular system? info=",info_minus)
  
  # 8) If fuzzy => do T
  alphaT_plus  <- NULL
  alphaT_minus <- NULL
  denom_T <- NA_real_
  
  if(fuzzy){
    T_mat <- matrix(T, nrow=n, ncol=1)
    
    outTplus <- .Fortran("locweights",
                         X       = as.double(X_centered),
                         YMAT    = as.double(T_mat),
                         N       = N,
                         P       = P_,
                         H       = as.double(h_use),
                         SIDE    = as.integer(1),
                         KERNELW = as.double(side_plus_w),
                         ALPHA   = double((p+1)),
                         WINT    = double(n),
                         INFO    = integer(1),
                         NQ      = as.integer(1),
                         PACKAGE="R3D")
    alphaT_plus <- matrix(outTplus$ALPHA, nrow=p+1, ncol=1)
    
    outTminus <- .Fortran("locweights",
                          X       = as.double(X_centered),
                          YMAT    = as.double(T_mat),
                          N       = N,
                          P       = P_,
                          H       = as.double(h_use),
                          SIDE    = as.integer(0),
                          KERNELW = as.double(side_minus_w),
                          ALPHA   = double((p+1)),
                          WINT    = double(n),
                          INFO    = integer(1),
                          NQ      = as.integer(1),
                          PACKAGE="R3D")
    alphaT_minus <- matrix(outTminus$ALPHA, nrow=p+1, ncol=1)
    
    denom_T <- alphaT_plus[1,1] - alphaT_minus[1,1]
    if(abs(denom_T) < 1e-14){
      warning("Denominator in fuzzy RDD is near 0 => effects might be NA.")
    }
  }
  
  ## ----------------------------------------------------
  ## (A) Compute unprojected fits E[Y(q)|X=x_i], row by row
  ## ----------------------------------------------------
  X_scaled <- X_centered / h_use
  Eplus_unproj  <- matrix(0, nrow=n, ncol=nQ)
  Eminus_unproj <- matrix(0, nrow=n, ncol=nQ)
  
  for(qi in seq_len(nQ)){
    ap <- alpha_plus[,qi]   # (p+1) coefs
    am <- alpha_minus[,qi]
    
    # plus side: pick all i for which side_plus_w>0
    idxp <- which(side_plus_w > 0)
    if(length(idxp)>0){
      Xpow_p <- outer(X_scaled[idxp], 0:p, "^")  # n_i x (p+1)
      Eplus_unproj[idxp, qi] <- Xpow_p %*% ap    # n_i x 1
    }
    
    # minus side
    idxm <- which(side_minus_w > 0)
    if(length(idxm)>0){
      Xpow_m <- outer(X_scaled[idxm], 0:p, "^")
      Eminus_unproj[idxm, qi] <- Xpow_m %*% am
    }
  }
  
  ## ----------------------------------------------------
  ## (B) If frechet => row-wise isoreg across quantiles
  ## ----------------------------------------------------
  Eplus_final  <- Eplus_unproj
  Eminus_final <- Eminus_unproj
  
  if(method=="frechet"){
    for(i in seq_len(n)){
      # plus
      if(side_plus_w[i] > 0){ 
        rowfit <- Eplus_unproj[i,]  # length nQ
        isoobj <- stats::isoreg(q_grid, rowfit)
        Eplus_final[i,] <- as.numeric(isoobj$yf)
      }
      # minus
      if(side_minus_w[i] > 0){
        rowfit <- Eminus_unproj[i,]
        isoobj <- stats::isoreg(q_grid, rowfit)
        Eminus_final[i,] <- as.numeric(isoobj$yf)
      }
    }
  }
  
  ## ----------------------------------------------------
  ## (C) Build outcome residuals
  ## ----------------------------------------------------
  e1_mat <- matrix(0, nrow=n, ncol=nQ)
  for(qi in seq_len(nQ)){
    idxp <- which(side_plus_w>0)
    if(length(idxp)>0){
      e1_mat[idxp, qi] <- Qmat[idxp, qi] - Eplus_final[idxp, qi]
    }
    idxm <- which(side_minus_w>0)
    if(length(idxm)>0){
      e1_mat[idxm, qi] <- Qmat[idxm, qi] - Eminus_final[idxm, qi]
    }
  }
  
  ## (D) If fuzzy => T residuals
  e2_mat <- NULL
  if(fuzzy){
    e2_mat <- matrix(0, nrow=n, ncol=nQ)
    
    # plus side
    idxp <- which(side_plus_w>0)
    if(length(idxp)>0){
      Xpow_p <- outer(X_scaled[idxp], 0:p, "^")
      fitpT  <- Xpow_p %*% alphaT_plus  # length idxp
      residp <- T[idxp] - fitpT
      # replicate across nQ
      e2_mat[idxp,] <- matrix(rep(residp, nQ), ncol=nQ)
    }
    # minus side
    idxm <- which(side_minus_w>0)
    if(length(idxm)>0){
      Xpow_m <- outer(X_scaled[idxm], 0:p, "^")
      fitmT  <- Xpow_m %*% alphaT_minus
      residm <- T[idxm] - fitmT
      e2_mat[idxm,] <- matrix(rep(residm, nQ), ncol=nQ)
    }
  }
  
  ## ----------------------------------------------------
  ## (E) Compute final LAQTE from intercepts alpha_+(1,qi) minus alpha_-(1,qi)
  ##     If method="frechet", we do isoreg on the intercept vector
  ## ----------------------------------------------------
  int_plus  <- alpha_plus[1,]  # length nQ
  int_minus <- alpha_minus[1,]
  
  if(method=="frechet"){
    # monotonic rearrangement of these intercepts
    iso_p <- stats::isoreg(q_grid, int_plus)
    int_plus <- as.numeric(iso_p$yf)
    iso_m <- stats::isoreg(q_grid, int_minus)
    int_minus <- as.numeric(iso_m$yf)
  }
  
  tau_vec <- int_plus - int_minus
  if(fuzzy){
    if(abs(denom_T)<1e-14){
      tau_vec[] <- NA_real_
    } else {
      tau_vec <- tau_vec / denom_T
    }
  }
  
  ## ----------------------------------------------------
  ## Build the output
  out <- list(
    results = list(
      tau        = tau_vec,
      q_grid     = q_grid,
      method     = method,
      fuzzy      = fuzzy,
      p          = p,
      bandwidths = bwres$h_star,
      h_used     = h_use
    ),
    
    coefficients = list(
      alpha_plus    = alpha_plus,
      alpha_minus   = alpha_minus,
      alphaT_plus   = alphaT_plus,
      alphaT_minus  = alphaT_minus
    ),
    
    bootstrap = list(
      w_plus  = w_plus,
      w_minus = w_minus,
      e1_mat  = e1_mat,
      e2_mat  = e2_mat
    ),
    
    # Return the row-by-row fitted means (projected if frechet)
    conditional_means = list(
      plus  = int_plus,
      minus = int_minus
    ),
    
    inputs = list(
      X      = X,
      Y_list = Y_list,
      T      = T,
      cutoff = cutoff,
      call   = match.call()
    ),
    
    diagnostics = list(
      info_plus   = info_plus,
      info_minus  = info_minus,
      denominator = if(fuzzy) denom_T else NULL
    )
  )
  
  # Also replicate top-level items for convenience
  out$tau         <- tau_vec
  out$q_grid      <- q_grid
  out$method      <- method
  out$fuzzy       <- fuzzy
  out$p           <- p
  out$bandwidths  <- bwres$h_star
  
  out$alpha_plus  <- alpha_plus
  out$alpha_minus <- alpha_minus
  out$alphaT_plus <- alphaT_plus
  out$alphaT_minus<- alphaT_minus
  out$int_plus <- int_plus
  out$int_minus <- int_minus
  
  out$w_plus      <- w_plus
  out$w_minus     <- w_minus
  out$e1_mat      <- e1_mat
  out$e2_mat      <- e2_mat
  
  out$X           <- X
  out$Y_list      <- Y_list
  out$T           <- T
  out$cutoff      <- cutoff
  out$call        <- match.call()
  
  class(out) <- "r3d"
  
  # optional bootstrap
  if(boot){
    boot_out <- r3d_bootstrap(
      object   = out,
      X        = X,
      Y_list   = Y_list,
      T        = T,
      B        = boot_reps,
      alpha    = alpha,
      test     = test,
      cores    = boot_cores
    )
    out$bootstrap$boot_out <- boot_out
    out$boot_out <- boot_out
  }
  
  out
}