# File: R3D/R/r3d_bwselect.R

#' @title Bandwidth Selection for R3D
#'
#' @description
#' Provides \code{bwselect()} to compute MSE-optimal (per-quantile) or IMSE-optimal 
#' (single) bandwidth(s) for the distributional RDD, without doing a full R3D estimation.
#'
#' @param X Numeric vector of running variables (\eqn{n}).
#' @param Y_list List of length \eqn{n}, each containing draws from unit i's distribution.
#' @param T (Optional) if fuzzy design (not used here, but included for consistency).
#' @param cutoff Numeric cutoff (default 0).
#' @param method "simple" => MSE bandwidths, "frechet" => IMSE single bandwidth.
#' @param p Polynomial order (1 or 2).
#' @param q_grid Numeric vector of quantiles in (0,1).
#' @param pilot_h (optional) pilot bandwidth for derivative/variance. If NULL, defaults to 0.1*n^{-1/(2(p+1)+3)}.
#' @param kernel Kernel function, triangular by default.
#'
#' @return A list with elements:
#' \item{h_star}{ The recommended bandwidth(s). For "simple", a vector \code{length(q_grid)}. 
#'    For "frechet", a single scalar.}
#' \item{pilot_h}{ The pilot bandwidth used.}
#' \item{method}{ "simple" or "frechet".}
#' \item{q_grid}{ The same quantile grid.}
#'
#' @export
bwselect <- function(X, Y_list, T=NULL,
                     cutoff=0,
                     method=c("simple","frechet"),
                     p=2,
                     q_grid=seq(0.05,0.95,by=0.05),
                     pilot_h=NULL,
                     kernel=function(u) pmax(0,1-abs(u))) 
{
  method <- match.arg(method)
  n <- length(X)
  if (length(Y_list)!=n) stop("bwselect: mismatch Y_list vs X length.")
  
  # compute empirical quantiles
  Qmat <- .compute_empirical_qmat(Y_list, q_grid)
  nq <- length(q_grid)
  
  if(is.null(pilot_h)) {
    pilot_h <- 0.1* n^(-1/(2*(p+1)+3))
  }
  
  # shift X for the pilot
  Xs <- X - cutoff
  Bplus <- numeric(nq)
  Bminus<- numeric(nq)
  Vplus <- numeric(nq)
  Vminus<- numeric(nq)
  
  # fill pilot estimates for each q
  for(j in seq_len(nq)) {
    qvals_j <- vapply(Qmat, `[`, 1.0, j)
    leftp  <- .pilot_side(Xs, qvals_j, kernel, pilot_h, side="left",  p=p)
    rightp <- .pilot_side(Xs, qvals_j, kernel, pilot_h, side="right", p=p)
    Bminus[j] <- leftp$bias
    Bplus[j]  <- rightp$bias
    Vminus[j] <- leftp$var
    Vplus[j]  <- rightp$var
  }
  
  if(method=="simple") {
    # MSE => vector
    h_star_q <- .lp_bw_mse_eachq(Bplus,Bminus,Vplus,Vminus,p,n)
    return(list(h_star=h_star_q, pilot_h=pilot_h, method="simple", q_grid=q_grid))
  } else {
    # frechet => single IMSE
    h_star <- .lp_bw_imse(q_grid, Bplus,Bminus,Vplus,Vminus,p,n)
    return(list(h_star=h_star, pilot_h=pilot_h, method="frechet", q_grid=q_grid))
  }
}


#=================== Internal pilot code ===========================

.pilot_side <- function(Xs, Qvals, kernel, pilot_h, side=c("left","right"), p) {
  side <- match.arg(side)
  n <- length(Xs)
  side_ok <- if(side=="left") (Xs<0) else (Xs>=0)
  U <- Xs/pilot_h
  kv <- numeric(n)
  for(i in seq_len(n)){
    if(side_ok[i]){
      kv[i] = kernel(U[i])
    } else {
      kv[i] = 0
    }
  }
  # local poly of order p+1 => dimension p+1
  # The last coefficient is ~ the derivative. We'll do a small matrix approach
  coefs <- .lp_boundary_coefs_pilot(Xs, Qvals, kv, pilot_h, ord=p+1)
  bias_approx <- coefs[length(coefs)]
  
  # variance => var of Qvals for |Xs|< pilot_h & side_ok
  sel_small <- side_ok & (abs(Xs)<pilot_h)
  var_est <- stats::var(Qvals[sel_small], na.rm=TRUE)
  if(is.na(var_est)) var_est<-0
  list(bias=bias_approx, var=var_est)
}

.lp_boundary_coefs_pilot <- function(Xs, Y, kv, h, ord) {
  n <- length(Xs)
  d <- ord
  XWX <- matrix(0,d,d)
  XWy <- numeric(d)
  for(i in seq_len(n)){
    if(kv[i]==0) next
    double w = kv[i];
    # R[i,:]
    double u = Xs[i]/h;
    for(a in seq_len(d)){
      for(b in seq_len(d)){
        XWX[a,b] = XWX[a,b] + (u^(a-1))*(u^(b-1))* w
      }
      XWy[a] = XWy[a] + (u^(a-1))* Y[i]* w
    }
  }
  sol <- .inv_smallR(XWX)
  if(is.null(sol)) return(rep(0,d))
  coefs <- as.vector(sol %*% XWy)
  coefs
}

.inv_smallR <- function(M) {
  d <- nrow(M)
  if(d!=ncol(M)) return(NULL)
  if(d<1 || d>3) return(NULL)
  dettol <- 1e-14
  if(d==1){
    if(abs(M[1,1])<dettol) return(NULL)
    return(matrix(1/M[1,1],1,1))
  } else if(d==2){
    detv <- M[1,1]*M[2,2] - M[1,2]*M[2,1]
    if(abs(detv)<dettol) return(NULL)
    invd <- 1/detv
    out <- matrix(0,2,2)
    out[1,1] =  M[2,2]*invd
    out[1,2] = -M[1,2]*invd
    out[2,1] = -M[2,1]*invd
    out[2,2] =  M[1,1]*invd
    return(out)
  } else {
    # 3x3
    a00 <- M[1,1]; a01 <- M[1,2]; a02 <- M[1,3]
    a10 <- M[2,1]; a11 <- M[2,2]; a12 <- M[2,3]
    a20 <- M[3,1]; a21 <- M[3,2]; a22 <- M[3,3]
    detv <- a00*(a11*a22 - a12*a21) 
    - a01*(a10*a22 - a12*a20)
    + a02*(a10*a21 - a11*a20)
    if(abs(detv)<dettol) return(NULL)
    invd <- 1/detv
    out <- matrix(0,3,3)
    out[1,1] = (a11*a22 - a12*a21)*invd
    out[1,2] = (a02*a21 - a01*a22)*invd
    out[1,3] = (a01*a12 - a02*a11)*invd
    out[2,1] = (a12*a20 - a10*a22)*invd
    out[2,2] = (a00*a22 - a02*a20)*invd
    out[2,3] = (a02*a10 - a00*a12)*invd
    out[3,1] = (a10*a21 - a11*a20)*invd
    out[3,2] = (a01*a20 - a00*a21)*invd
    out[3,3] = (a00*a11 - a01*a10)*invd
    return(out)
  }
}


# MSE bandwidth per q
.lp_bw_mse_eachq <- function(Bplus,Bminus,Vplus,Vminus,p,n) {
  nq <- length(Bplus)
  out <- numeric(nq)
  for(j in seq_len(nq)) {
    num <- Vplus[j] + Vminus[j]
    diffB <- (Bplus[j]-Bminus[j])^2
    den <- 2*(p+1)* diffB
    if(abs(den)<1e-14){
      out[j]<-0
    } else {
      ratio <- num/den
      out[j] <- ( ratio )^(1/(2*p+3)) * n^(-1/(2*p+3))
      if(!is.finite(out[j])|| out[j]<0) out[j]<-0
    }
  }
  out
}

# IMSE bandwidth single
.lp_bw_imse <- function(q_grid, Bplus,Bminus,Vplus,Vminus,p,n) {
  if(length(q_grid)<2) return(0)
  dq <- diff(range(q_grid))/(length(q_grid)-1)
  diffsB <- (Bplus - Bminus)
  Ap <- sum(diffsB^2)*dq
  Bp <- sum(Vplus + Vminus)*dq
  if(abs(Ap)<1e-14) return(0)
  ratio <- Bp /(2*(p+1)* Ap)
  hstar <- (ratio)^(1/(2*p+3)) * n^(-1/(2*p+3))
  if(!is.finite(hstar)|| hstar<0) hstar=0
  hstar
}
