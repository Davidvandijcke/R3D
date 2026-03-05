      subroutine locweights(X, YMAT, N, P, H, SIDE, KERNEL_TYPE,
     +                      ALPHA, WINT, INFO, NQ)
C
C  Purpose:
C    Compute local-polynomial regression coefficients and intercept
C    weights for a set of NQ quantiles, for E[Y(q) | X=0] using
C    polynomial order P (one-sided "SIDE"), with kernel weights
C    computed internally based on KERNEL_TYPE and H(q).
C
C  Inputs:
C    X(N): the data for X (centered at cutoff)
C    YMAT(N*NQ): each row i is the vector [Q_{Y_i}(q=1),..., Q_{Y_i}(q=NQ)],
C                stored in column-major order
C    P: polynomial order (0.. e.g. 2)
C    H(NQ): bandwidth vector, one value per quantile
C    SIDE: 1 => plus side (X>=0), 0 => minus side (X<0)
C    KERNEL_TYPE: integer (1=triangular, 2=epanechnikov, 3=uniform)
C    NQ: number of quantiles
C
C  Outputs:
C    ALPHA((P+1)*NQ): the local-poly regression coefficients for each q, flat array
C    WINT(N*NQ): intercept weights for each observation i and quantile q, flat array
C    INFO: =0 if success, <0 if error:
C           -94 N < dim (too few observations for polynomial order)
C           -95 NQ < 1 (no quantiles requested)
C           -96 invalid SIDE (must be 0 or 1)
C           -97 zero or negative bandwidth
C           -98 unknown kernel type
C           -99 polynomial degree exceeds MAXDIM
C           positive: singular system in DGETRF/DGETRS (per-quantile cycle)
C
      implicit none
C     Input/output variables
      integer N, P, SIDE, NQ, INFO, KERNEL_TYPE
      double precision X(N)
      double precision YMAT(N*NQ)  ! Flat array in column-major order
      double precision H(NQ)       ! Array of bandwidths
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
      double precision w_i, xx, u
      integer IPIV(MAXDIM)
      double precision w_int
      double precision ZERO, ONE
      parameter (ZERO = 0.0d0, ONE = 1.0d0)
      character TRANS
      parameter (TRANS = 'N')
      double precision W(N)   ! F4: kernel weight cache to avoid recomputation
      double precision FZERO  ! For runtime NaN (avoids compile-time constant error)
C
C     Initialize
      INFO = 0
      dim = P + 1
      FZERO = 0.0d0

C     Basic dimension check
      if (dim .gt. MAXDIM) then
         INFO = -99
         return
      endif
C
C     F5: Input validation
      if (SIDE .ne. 0 .and. SIDE .ne. 1) then
         INFO = -96
         return
      endif
      if (NQ .lt. 1) then
         INFO = -95
         return
      endif
      if (N .lt. dim) then
         INFO = -94
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
C
C        F2: Guard against zero or negative bandwidth (divide-by-zero)
         if (H(q) .le. ZERO) then
            INFO = -97
            return
         endif
C
C        2. Accumulate X'WX and X'WY for this quantile q
         do i = 1, N
C           Compute kernel weight based on X(i)/H(q) and kernel type
            u = X(i) / H(q)
            if (KERNEL_TYPE .eq. 1) then
C              Triangular: K(u) = max(0, 1 - |u|)
               if (abs(u) .le. ONE) then
                  w_i = ONE - abs(u)
               else
                  w_i = ZERO
               endif
            else if (KERNEL_TYPE .eq. 2) then
C              Epanechnikov: K(u) = 0.75 * max(0, 1 - u^2)
               if (abs(u) .le. ONE) then
                  w_i = 0.75d0 * (ONE - u*u)
               else
                  w_i = ZERO
               endif
            else if (KERNEL_TYPE .eq. 3) then
C              Uniform: K(u) = 0.5 if |u| <= 1, 0 otherwise
               if (abs(u) .le. ONE) then
                  w_i = 0.5d0
               else
                  w_i = ZERO
               endif
            else
               INFO = -98  ! Unknown kernel type
               return
            endif
C
C           Apply side restriction
            if (SIDE .eq. 1 .and. X(i) .lt. ZERO) then
               w_i = ZERO
            else if (SIDE .eq. 0 .and. X(i) .ge. ZERO) then
               w_i = ZERO
            endif
C
C           F4: Cache final kernel weight for reuse in intercept pass
            W(i) = w_i
C
C           Skip if zero weight
            if (w_i .eq. ZERO) then
               continue
            else
C              Compute basis functions using bandwidth for this quantile
               xx = X(i) / H(q)
               basis(1) = ONE
               do k = 2, dim
                  basis(k) = basis(k-1) * xx
               enddo

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
         call DGETRF(dim, dim, MAT, MAXDIM, IPIV, INFO)
         if (INFO .ne. 0) then
C           F3: Singular system -- mark this quantile NaN and cycle to next
            INFO = 0
            do j = 1, dim
               jdx = (q-1)*dim + j
               ALPHA(jdx) = FZERO/FZERO
            enddo
            do i = 1, N
               idx = (q-1)*N + i
               WINT(idx) = FZERO/FZERO
            enddo
            cycle
         endif

         call DGETRS(TRANS, dim, 1, MAT, MAXDIM, IPIV, RHS, dim, INFO)

C        4. Store the coefficients
         do j = 1, dim
            jdx = (q-1)*dim + j
            ALPHA(jdx) = RHS(j)
         enddo

C        5. Compute intercept weights using inverse row
         do j = 1, dim
            RHS(j) = ZERO
         enddo
         RHS(1) = ONE

         call DGETRS(TRANS, dim, 1, MAT, MAXDIM, IPIV, RHS, dim, INFO)

C           Compute intercept weights for each observation
         do i = 1, N
C           F4: Reuse cached kernel weight from first pass
            w_i = W(i)
C
C           Skip if zero weight
            if (w_i .eq. ZERO) then
               continue
            else
C              Recompute basis using bandwidth for this quantile
               xx = X(i) / H(q)
               basis(1) = ONE
               do k = 2, dim
                  basis(k) = basis(k-1) * xx
               enddo

C              Dot product with first row of inverse
               w_int = ZERO
               do k = 1, dim
                  w_int = w_int + RHS(k) * basis(k)
               enddo

C              Scale by kernel weight
               w_int = w_int * w_i

C              Store in output array
               idx = (q-1)*N + i
               WINT(idx) = w_int
            endif
         enddo
      enddo
C
      return
      end
