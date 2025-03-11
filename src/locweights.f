C ====================================================================
      subroutine locweights(X, YMAT, N, P, H, SIDE, KERNELW,
     +                      ALPHA, WINT, INFO, NQ)
C
C  Purpose:
C    Compute local-polynomial regression coefficients and intercept
C    weights for a set of NQ quantiles, for E[Y(q) | X=0] using
C    polynomial order P (one-sided "SIDE").
C
C  Inputs:
C    X(N), KERNELW(N): the data for X and kernel weights K(...) at X_i/h
C    YMAT(N*NQ): each row i is the vector [Q_{Y_i}(q=1),..., Q_{Y_i}(q=NQ)]
C                stored in column-major order
C    P: polynomial order (0.. e.g. 2)
C    H: bandwidth
C    SIDE: 1 => plus side (X>=0), 0 => minus side (X<0)
C    NQ: number of quantiles
C
C  Outputs:
C    ALPHA((p+1)*NQ): the local-poly regression coefficients for each q, flat array
C    WINT(N*NQ): intercept weights for each observation i and quantile q, flat array
C    INFO: =0 if success, !=0 if singular
C
      implicit none
C     Input/output variables
      integer N, P, SIDE, NQ, INFO
      double precision X(N)
      double precision YMAT(N*NQ)  ! Flat array in column-major order
      double precision H
      double precision KERNELW(N)
      double precision ALPHA((P+1)*NQ)  ! Flat output
      double precision WINT(N*NQ)       ! Flat output
C
C     Fixed size working arrays
      integer MAXDIM
      parameter (MAXDIM = 10)  ! Max polynomial degree
C
C     Local variables
      integer i, j, k, q, idx, jdx, dim
      double precision MAT(MAXDIM, MAXDIM)
      double precision RHS(MAXDIM)
      double precision basis(MAXDIM)
      double precision w_i, xx
      integer IPIV(MAXDIM)
      double precision ZERO, ONE
      parameter (ZERO = 0.0d0, ONE = 1.0d0)
      character TRANS
      parameter (TRANS = 'N')
C
C     Initialize
      INFO = 0
      dim = P + 1
      
C     Basic dimension check
      if (dim .gt. MAXDIM) then
         INFO = -99
         return
      endif
C
C     Zero out result arrays
      do i = 1, (P+1)*NQ
         ALPHA(i) = ZERO
      enddo
      
      do i = 1, N*NQ
         WINT(i) = ZERO
      enddo
C
C     Process one quantile at a time
      do q = 1, NQ
C        1. Reset the matrix and RHS for each quantile
         do i = 1, dim
            do j = 1, dim
               MAT(i,j) = ZERO
            enddo
            RHS(i) = ZERO
         enddo
         
C        2. Accumulate X'WX and X'WY for this quantile q
         do i = 1, N
C           Skip if on wrong side or zero weight
            if ((SIDE .eq. 1 .and. X(i) .lt. ZERO) .or.
     +          (SIDE .eq. 0 .and. X(i) .ge. ZERO) .or.
     +          (KERNELW(i) .eq. ZERO)) then
               continue
            else
C              Compute basis functions
               xx = X(i) / H
               basis(1) = ONE
               do k = 2, dim
                  basis(k) = basis(k-1) * xx
               enddo
               
C              Weight
               w_i = KERNELW(i)
               
C              Accumulate X'WX
               do j = 1, dim
                  do k = 1, dim
                     MAT(j,k) = MAT(j,k) + w_i * basis(j) * basis(k)
                  enddo
               enddo
               
C              Accumulate X'WY using column-major indexing
               idx = (q-1)*N + i
               do j = 1, dim
                  RHS(j) = RHS(j) + w_i * basis(j) * YMAT(idx)
               enddo
            endif
         enddo
         
C        3. Solve the system for this quantile
C           Use LAPACK DGETRF instead of LINPACK DGEFA
         call DGETRF(dim, dim, MAT, MAXDIM, IPIV, INFO)
         if (INFO .ne. 0) then
            return  ! Singular matrix
         endif
         
C           Use LAPACK DGETRS instead of LINPACK DGESL
         call DGETRS(TRANS, dim, 1, MAT, MAXDIM, IPIV, RHS, dim, INFO)
         
C        4. Store the coefficients
         do j = 1, dim
            jdx = (q-1)*dim + j
            ALPHA(jdx) = RHS(j)
         enddo
         
C        5. Compute intercept weights using inverse row
C           Prepare unit vector e1 = (1,0,0,...)
         do j = 1, dim
            RHS(j) = ZERO
         enddo
         RHS(1) = ONE
         
C           Solve to get first row of inverse using LAPACK
         call DGETRS(TRANS, dim, 1, MAT, MAXDIM, IPIV, RHS, dim, INFO)
         
C           Compute intercept weights for each observation
         do i = 1, N
C              Skip if on wrong side or zero weight
            if ((SIDE .eq. 1 .and. X(i) .lt. ZERO) .or.
     +          (SIDE .eq. 0 .and. X(i) .ge. ZERO) .or.
     +          (KERNELW(i) .eq. ZERO)) then
               continue
            else
C              Recompute basis
               xx = X(i) / H
               basis(1) = ONE
               do k = 2, dim
                  basis(k) = basis(k-1) * xx
               enddo
               
C              Dot product with first row of inverse
               w_i = ZERO
               do k = 1, dim
                  w_i = w_i + RHS(k) * basis(k)
               enddo
               
C              Scale by kernel weight
               w_i = w_i * KERNELW(i)
               
C              Store in output array
               idx = (q-1)*N + i
               WINT(idx) = w_i
            endif
         enddo
      enddo
C
      return
      end