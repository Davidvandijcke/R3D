module r3d_locweights
  implicit none
  
contains

  subroutine locweights(X, YMAT, N, P, H, SIDE, KERNEL_TYPE, &
                        ALPHA, WINT, INFO, NQ)
  !
  ! Purpose:
  !   Compute local-polynomial regression coefficients and intercept
  !   weights for a set of NQ quantiles, for E[Y(q) | X=0] using
  !   polynomial order P (one-sided "SIDE"), with kernel weights
  !   computed internally based on KERNEL_TYPE and H(q).
  !
  ! Inputs:
  !   X(N): the data for X (centered at cutoff)
  !   YMAT(N,NQ): each row i is the vector [Q_{Y_i}(q=1),..., Q_{Y_i}(q=NQ)]
  !   P: polynomial order (0.. e.g. 2)
  !   H(NQ): bandwidth vector, one value per quantile
  !   SIDE: 1 => plus side (X>=0), 0 => minus side (X<0)
  !   KERNEL_TYPE: integer (1=triangular, 2=epanechnikov, 3=uniform)
  !   NQ: number of quantiles
  !
  ! Outputs:
  !   ALPHA(P+1,NQ): the local-poly regression coefficients for each q
  !   WINT(N,NQ): intercept weights for each observation i and quantile q
  !   INFO: =0 if success, !=0 if singular
  !
    implicit none
    ! Input/output variables
    integer, intent(in) :: N, P, SIDE, NQ, KERNEL_TYPE
    real(8), intent(in) :: X(N)
    real(8), intent(in) :: YMAT(N,NQ)
    real(8), intent(in) :: H(NQ)
    real(8), intent(out) :: ALPHA(P+1,NQ)
    real(8), intent(out) :: WINT(N,NQ)
    integer, intent(out) :: INFO
    
    ! Local variables
    integer :: i, j, k, q, dim
    real(8) :: MAT(P+1,P+1)
    real(8) :: RHS(P+1)
    real(8) :: basis(P+1)
    real(8) :: w_i, xx, u, w_int
    integer :: IPIV(P+1)
    real(8), parameter :: ZERO = 0.0d0, ONE = 1.0d0
    
    ! Initialize
    INFO = 0
    dim = P + 1
    
    ! Zero out result arrays
    ALPHA = ZERO
    WINT = ZERO
    
    ! Process one quantile at a time
    do q = 1, NQ
      ! 1. Reset the matrix and RHS for each quantile
      MAT = ZERO
      RHS = ZERO
      
      ! 2. Accumulate X'WX and X'WY for this quantile q
      do i = 1, N
        ! Compute kernel weight based on X(i)/H(q) and kernel type
        u = X(i) / H(q)
        
        select case (KERNEL_TYPE)
        case (1)  ! Triangular: K(u) = max(0, 1 - |u|)
          if (abs(u) <= ONE) then
            w_i = ONE - abs(u)
          else
            w_i = ZERO
          endif
        case (2)  ! Epanechnikov: K(u) = 0.75 * max(0, 1 - u^2)
          if (abs(u) <= ONE) then
            w_i = 0.75d0 * (ONE - u*u)
          else
            w_i = ZERO
          endif
        case (3)  ! Uniform: K(u) = 0.5 if |u| <= 1, 0 otherwise
          if (abs(u) <= ONE) then
            w_i = 0.5d0
          else
            w_i = ZERO
          endif
        case default
          INFO = -98  ! Unknown kernel type
          return
        end select
        
        ! Apply side restriction
        if (SIDE == 1 .and. X(i) < ZERO) then
          w_i = ZERO
        else if (SIDE == 0 .and. X(i) >= ZERO) then
          w_i = ZERO
        endif
        
        ! Skip if zero weight
        if (w_i /= ZERO) then
          ! Compute basis functions using bandwidth for this quantile
          xx = X(i) / H(q)
          basis(1) = ONE
          do k = 2, dim
            basis(k) = basis(k-1) * xx
          enddo
          
          ! Accumulate X'WX
          do j = 1, dim
            do k = 1, dim
              MAT(j,k) = MAT(j,k) + w_i * basis(j) * basis(k)
            enddo
          enddo
          
          ! Accumulate X'WY
          do j = 1, dim
            RHS(j) = RHS(j) + w_i * basis(j) * YMAT(i,q)
          enddo
        endif
      enddo
      
      ! 3. Solve the system for this quantile
      call DGETRF(dim, dim, MAT, dim, IPIV, INFO)
      if (INFO /= 0) return  ! Singular matrix
      
      call DGETRS('N', dim, 1, MAT, dim, IPIV, RHS, dim, INFO)
      
      ! 4. Store the coefficients
      ALPHA(:,q) = RHS
      
      ! 5. Compute intercept weights using inverse row
      RHS = ZERO
      RHS(1) = ONE
      
      call DGETRS('N', dim, 1, MAT, dim, IPIV, RHS, dim, INFO)
      
      ! Compute intercept weights for each observation
      do i = 1, N
        ! Initialize weight
        w_i = ZERO
        
        ! Recompute kernel weight (same as above)
        u = X(i) / H(q)
        
        select case (KERNEL_TYPE)
        case (1)  ! Triangular
          if (abs(u) <= ONE) then
            w_i = ONE - abs(u)
          else
            w_i = ZERO
          endif
        case (2)  ! Epanechnikov
          if (abs(u) <= ONE) then
            w_i = 0.75d0 * (ONE - u*u)
          else
            w_i = ZERO
          endif
        case (3)  ! Uniform
          if (abs(u) <= ONE) then
            w_i = 0.5d0
          else
            w_i = ZERO
          endif
        case default
          w_i = ZERO
        end select
        
        ! Apply side restriction
        if (SIDE == 1 .and. X(i) < ZERO) then
          w_i = ZERO
        else if (SIDE == 0 .and. X(i) >= ZERO) then
          w_i = ZERO
        endif
        
        ! Skip if zero weight
        if (w_i /= ZERO) then
          ! Recompute basis using bandwidth for this quantile
          xx = X(i) / H(q)
          basis(1) = ONE
          do k = 2, dim
            basis(k) = basis(k-1) * xx
          enddo
          
          ! Dot product with first row of inverse
          w_int = ZERO
          do k = 1, dim
            w_int = w_int + RHS(k) * basis(k)
          enddo
          
          ! Scale by kernel weight and store
          WINT(i,q) = w_int * w_i
        endif
      enddo
    enddo
    
  end subroutine locweights

end module r3d_locweights

! C-callable wrapper for STATA plugin interface
subroutine stata_locweights(x, ymat, n, p, h, side, kernel_type, &
                           alpha, wint, info, nq) bind(C, name="stata_locweights_")
  use iso_c_binding
  use r3d_locweights
  implicit none
  
  ! C-style arguments (arrays passed as pointers)
  integer(c_int), intent(in), value :: n, p, side, kernel_type, nq
  real(c_double), intent(in) :: x(n), ymat(n*nq), h(nq)
  real(c_double), intent(out) :: alpha((p+1)*nq), wint(n*nq)
  integer(c_int), intent(out) :: info
  
  ! Local Fortran arrays with proper shape
  real(8) :: ymat_2d(n, nq), alpha_2d(p+1, nq), wint_2d(n, nq)
  integer :: i, j, idx
  
  ! Convert flat arrays to 2D arrays (column-major order)
  do j = 1, nq
    do i = 1, n
      idx = (j-1)*n + i
      ymat_2d(i,j) = ymat(idx)
    end do
  end do
  
  ! Call the main routine
  call locweights(x, ymat_2d, n, p, h, side, kernel_type, alpha_2d, wint_2d, info, nq)
  
  ! Convert 2D arrays back to flat arrays
  do j = 1, nq
    do i = 1, p+1
      idx = (j-1)*(p+1) + i
      alpha(idx) = alpha_2d(i,j)
    end do
  end do
  
  do j = 1, nq
    do i = 1, n
      idx = (j-1)*n + i
      wint(idx) = wint_2d(i,j)
    end do
  end do
  
end subroutine stata_locweights