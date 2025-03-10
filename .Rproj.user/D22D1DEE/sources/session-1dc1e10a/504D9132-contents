# File: R3D/R/r3d_simple.R
simple_r3d_impl <- function(X, Y_list, T, cutoff, p, q_grid, h_star_q, fuzzy, kernel) {
  n <- length(X)
  Qmat <- .compute_empirical_qmat(Y_list, q_grid)
  nq <- length(q_grid)
  m_minus_hat<- numeric(nq)
  m_plus_hat <- numeric(nq)
  if(fuzzy){
    T_minus_hat<- numeric(nq)
    T_plus_hat <- numeric(nq)
  }
  
  for(j in seq_len(nq)){
    qvals_j<- vapply(Qmat, `[`, 1.0, j)
    bw_j<- h_star_q[j]
    U<- (X-cutoff)/bw_j
    kv<- kernel(U)
    w_minus<- lpweights_rcpp(X, cutoff, TRUE, p, kv, bw_j)
    w_plus <- lpweights_rcpp(X, cutoff, FALSE,p, kv, bw_j)
    m_minus_hat[j]<- sum(w_minus* qvals_j)
    m_plus_hat[j] <- sum(w_plus * qvals_j)
    if(fuzzy){
      T_minus_hat[j]<- sum(w_minus*T)
      T_plus_hat[j] <- sum(w_plus *T)
    }
  }
  if(fuzzy){
    denom<- T_plus_hat - T_minus_hat
    tau <- (m_plus_hat - m_minus_hat)/denom
  } else {
    tau <- m_plus_hat - m_minus_hat
  }
  list(tau=tau, m_minus=m_minus_hat, m_plus=m_plus_hat)
}
