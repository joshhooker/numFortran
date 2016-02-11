module matrixSolver
  implicit none

  private
  integer, parameter :: dp = kind(0.d0)

  public :: tridiagonalSolver

contains

  subroutine tridiagonalSolver(a,b,c,f,u)
    !! Tridiagonal matrix solver that solves
    !! a tridiagonal system for n unknowns may be written as
    !! \[a_i u_{i-1} + b_i u_i + c_i u_{i+1} = f_i\]
    !! where \(a_1 = 0\)  and \(c_n = 0\)
    !! \[ \]
    !! The tridiagonal matrix looks like:
    !! \[ \left[ \begin{array}{ccccccc}
    !! b_1 & c_1 & 0 & \dots & & & 0 \\\
    !! a_2 & b_2 & c_2 & \dots & & \unicode{x22f0} & \vdots \\\
    !! & & & \dots & & & \\\
    !! \vdots & \unicode{x22f0} & & \dots & a_{N-1} & b_{N-1} & c_{N-1} \\\
    !! 0 & & & \dots & 0 & a_N & b_N \\
    !! \end{array} \right] \cdot
    !! \left[ \begin{array}{c} u_1 \\\ u_2 \\\ \vdots \\\ u_{N-1} \\\ u_N \end{array} \right]
    !! =
    !! \left[ \begin{array}{c} f_1 \\\ f_2 \\\ \vdots \\\ f_{N-1} \\\ f_N \end{array} \right] \]

    integer :: i, n
    real(dp) :: a(:)
      !! input: array a of length n for tridiagonal matrix
    real(dp) :: b(:)
      !! input: array b of length n for tridiagonal matrix
    real(dp) :: c(:)
      !! input: array c of length n for tridiagonal matrix
    real(dp) :: f(:)
      !! input: array f of length n for tridiagonal matrix
    real(dp) :: u(:)
      !! output: array u of length n of the solution to the tridiagonal matrix
    real(dp) :: bp
    real(dp), allocatable :: cp(:)

    n = size(a,1)

    allocate(cp(n))

    bp = b(1)
    u(1) = f(1)/bp
    do i=2,n
      cp(i) = c(i-1)/bp
      bp = b(i) - a(i)*cp(i)
      u(i) = (f(i) - a(i)*u(i-1))/bp
    end do

    do i=n-1,1,-1
      u(i) = u(i) - cp(i+1)*u(i+1)
    end do

    deallocate(cp)
  end subroutine tridiagonalSolver

end module