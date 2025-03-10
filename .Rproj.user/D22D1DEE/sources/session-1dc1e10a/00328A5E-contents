// File: lpweights.cpp

#include <Rcpp.h>
using namespace Rcpp;

/*
  We compute local polynomial boundary weights for extracting the intercept 
  at x=cutoff, side=left or right, p=1 or p=2. 
  Steps:
    1) Identify side (X<cutoff vs X>=cutoff).
    2) Build R[i,:] for each obs: [1, U, U^2] if p=2, with U=(X[i]-cutoff)/h.
    3) Accumulate X'WX and X'W, invert, multiply to get the intercept weights.
    4) Return an n-vector of weights for each observation i.
*/

// inline helpers to invert 2x2 or 3x3 quickly
inline bool inv2(const double *M, double *Minv) {
  double det = M[0]*M[3] - M[1]*M[2];
  if(std::fabs(det)<1e-14) return false;
  double invdet = 1.0/det;
  Minv[0] =  M[3]*invdet;
  Minv[1] = -M[1]*invdet;
  Minv[2] = -M[2]*invdet;
  Minv[3] =  M[0]*invdet;
  return true;
}
inline bool inv3(const double *M, double *Minv) {
  double det = 
    M[0]*(M[4]*M[8] - M[5]*M[7]) 
    - M[1]*(M[3]*M[8] - M[5]*M[6])
    + M[2]*(M[3]*M[7] - M[4]*M[6]);
  if(std::fabs(det)<1e-14) return false;
  double invd = 1.0/det;
  Minv[0] = ( M[4]*M[8] - M[5]*M[7]) * invd;
  Minv[1] = ( M[2]*M[7] - M[1]*M[8]) * invd;
  Minv[2] = ( M[1]*M[5] - M[2]*M[4]) * invd;
  Minv[3] = ( M[5]*M[6] - M[3]*M[8]) * invd;
  Minv[4] = ( M[0]*M[8] - M[2]*M[6]) * invd;
  Minv[5] = ( M[2]*M[3] - M[0]*M[5]) * invd;
  Minv[6] = ( M[3]*M[7] - M[4]*M[6]) * invd;
  Minv[7] = ( M[1]*M[6] - M[0]*M[7]) * invd;
  Minv[8] = ( M[0]*M[4] - M[1]*M[3]) * invd;
  return true;
}

// [[Rcpp::export]]
NumericVector lpweights_rcpp(NumericVector X, double cutoff,
                             bool left_side, int p,
                             NumericVector kernvals,
                             double h) 
{
  // X: running var
  // cutoff
  // left_side => side=TRUE => X<cutoff
  // p=1 or p=2
  // kernvals: kernel( (X-cutoff)/h ) user-provided
  // h: bandwidth
  // returns n-vector of boundary weights for the intercept.

  int n = X.size();
  int dim = p+1; // 2 if p=1, 3 if p=2
  NumericVector wts(n, 0.0);

  // We'll accumulate XWX (dimxdim) row-major + XW (dimxn).
  std::vector<double> XWX(dim*dim, 0.0);
  std::vector<double> XW(dim*n, 0.0);

  for(int i=0;i<n;i++){
    bool side_ok = left_side ? (X[i]<cutoff) : (X[i]>=cutoff);
    double w_i = side_ok ? kernvals[i] : 0.0;
    if(w_i==0.0) continue;
    double u = (X[i] - cutoff)/h;
    // R[i,:]
    double r0 = 1.0;
    double r1 = (p>=1 ? u : 0.0);
    double r2 = (p>=2 ? u*u : 0.0);
    double r[3]; 
    r[0]=r0; r[1]=r1; r[2]=r2;

    // accumulate
    for(int a=0;a<dim;a++){
      for(int b=0;b<dim;b++){
        XWX[a*dim+b] += r[a]*r[b]* w_i;
      }
      XW[a*n + i] += (r[a]* w_i);
    }
  }

  // invert XWX
  std::vector<double> invXWX(dim*dim, 0.0);
  bool ok=false;
  if(dim==2) {
    ok=inv2(&XWX[0], &invXWX[0]);
  } else {
    ok=inv3(&XWX[0], &invXWX[0]);
  }
  if(!ok) {
    // singular => all wts=0
    return wts;
  }

  // alphaInv => first row of invXWX => intercept row
  // alphaInv[j] = invXWX[0*dim + j]
  std::vector<double> alphaInv(dim,0.0);
  for(int j=0;j<dim;j++){
    alphaInv[j] = invXWX[j];
  }

  // compute wts[i] => sum_j alphaInv[j]* XW[j,i]
  for(int i=0;i<n;i++){
    double sumv=0.0;
    for(int j=0;j<dim;j++){
      sumv += alphaInv[j]* XW[j*n + i];
    }
    wts[i]=sumv;
  }

  return wts;
}
