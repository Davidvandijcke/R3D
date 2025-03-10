# File: R3D/R/r3d_frechet.R

#' @title Internal: Frechet R3D (local polynomial + single IMSE bandwidth + isotonic projection)
#' @description
#' Usually called by r3d() if method="frechet".
#' @keywords internal
frechet_r3d_impl <- function(X, Y_list, T, cutoff, p, q_grid, h_star, fuzzy, kernel) {
  n <- length(X)
  Qmat <- .compute_empirical_qmat(Y_list, q_grid)
  nq <- length(q_grid)
  
  U <- (X - cutoff)/h_star
  kv <- kernel(U)
  w_minus <- lpweights_rcpp(X, cutoff, TRUE, p, kv, h_star)
  w_plus  <- lpweights_rcpp(X, cutoff, FALSE,p, kv, h_star)
  
  raw_minus <- numeric(nq)
  raw_plus  <- numeric(nq)
  for(j in seq_len(nq)) {
    qvals_j <- vapply(Qmat, `[`, 1.0, j)
    raw_minus[j] <- sum(w_minus*qvals_j)
    raw_plus[j]  <- sum(w_plus *qvals_j)
  }
  if(fuzzy){
    T_minus <- sum(w_minus*T)
    T_plus  <- sum(w_plus *T)
    denom   <- T_plus - T_minus
  }
  
  # monotone rearrangement
  proj_m_minus <- .monotone_projection(raw_minus, q_grid)
  proj_m_plus  <- .monotone_projection(raw_plus,  q_grid)
  
  if(!fuzzy){
    tau <- proj_m_plus - proj_m_minus
  } else {
    tau <- (proj_m_plus - proj_m_minus)/denom
  }
  list(tau=tau,
       raw_m_minus=raw_minus, raw_m_plus=raw_plus,
       proj_m_minus=proj_m_minus, proj_m_plus=proj_m_plus)
}

.monotone_projection <- function(vals, x) {
  if (!requireNamespace("isoreg", quietly=TRUE)) {
    iso <- stats::isoreg(x, vals)
    return(as.numeric(iso$yf))
  } else {
    out <- isoreg::isoreg(x, vals)
    return(as.numeric(out$yf))
  }
}
