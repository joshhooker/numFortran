!*********************************************************
!>
! Functions to compute natural cubic spline of the form
! \[S_i(x) = A_i + B_i*(x-x_i) + C_i*(x-x_i)^2 + D_i*(x-x_i)^3\]
! where \(A\), \(B\), \(C\) and \(D\) are arrays thats are computed from
! the data. For each step, i.e. between \(x_{i+1}\) and \(x_i\), a set of
! \(A_i\), \(B_i\), \(C_i\) and \(D_i\) is computed.

module cubicSpline
  implicit none

  private :: dp
  integer, parameter :: dp = kind(0.d0)

  contains

  subroutine tridiagonalSolver(a,b,c,f,u)
    !! date: January 8, 2016
    !! version: v0.1
    !!
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

  subroutine naturalCubicSpline(xi, yi, splineA, splineB, splineC, splineD)
    !! date: January 8, 2016
    !! version: v0.1
    !!
    !! Computes parameters used in the natural cubic spline
    !! of the form:
    !! \[S_i(x) = A_i + B_i*(x-x_i) + C_i*(x-x_i)^2 + D_i*(x-x_i)^3\]

    integer :: n
    real(dp) :: xi(:)
      !! input: array x, of length n, to compute spline variables
    real(dp) :: yi(:)
      !! input: array y, of length n, to compute spline variables
    real(dp) :: splineA(:)
      !! output: array splineA, of length n, to compute y at any x
    real(dp) :: splineB(:)
      !! output: array splineB, of length n, to compute y at any x
    real(dp) :: splineC(:)
      !! output: array splineC, of length n, to compute y at any x
    real(dp) :: splineD(:)
      !! output: array splineD, of length n, to compute y at any x
    real(dp), dimension(:), allocatable :: a,b,c,f,h

    n = size(xi,1)
    allocate(a(n), b(n), c(n), f(n), h(n))

    a(1) = 0.d0
    b(1) = 1.d0
    c(1) = 0.d0
    f(1) = 0.d0
    splineA(1:n) = yi(1:n)
    splineB(1:n) = 0.d0
    splineC(1:n) = 0.d0
    splineD(1:n) = 0.d0
    h(1:n) = abs(xi(2:n+1) - xi(1:n))
    a(2:n) = h(1:n-1)
    b(2:n) = 2.d0*(h(1:n-1) + h(2:n))
    c(2:n) = h(2:n)
    f(2:n) = abs(3.d0*(yi(3:n+1)-yi(2:n))/h(2:n) - 3.d0*(yi(2:n)-yi(1:n-1))/h(1:n-1))
    a(n) = 0.d0
    b(n) = 1.d0
    c(n) = 0.d0
    f(n) = 0.d0

    call tridiagonalSolver(a,b,c,f,splineC)

    splineD(1:n-1) = (splineC(2:n) - splineC(1:n-1))/(3.d0*h(1:n-1))
    splineB(1:n-1) = (yi(2:n) - yi(1:n-1))/h(1:n-1) - splineC(1:n-1)*h(1:n-1) -&
                      splineD(1:n-1)*h(1:n-1)*h(1:n-1)

    deallocate(a,b,c,f,h)
  end subroutine naturalCubicSpline

  real(dp) function computeNaturalCubicSpline(xi, x, splineA, splineB, splineC, splineD)
    !! date: January 8, 2016
    !! version: v0.1
    !!
    !! Computes y for a given x for a given set of cubic spline parameters.
    !! Returns a real(dp)

    integer :: i, n
    real(dp) :: xi(:)
      !! input: x array, of length n, which is the original input used to compute spline parameters
    real(dp) :: x
      !! input: x value to compute spline at
    real(dp) :: splineA(:)
      !! input: splineA array, of length n, computed by naturalCubicSpline
    real(dp) :: splineB(:)
      !! input: splineB array, of length n, computed by naturalCubicSpline
    real(dp) :: splineC(:)
      !! input: splineC array, of length n, computed by naturalCubicSpline
    real(dp) :: splineD(:)
      !! input: splineD array, of length n, computed by naturalCubicSpline
    real(dp) :: output, dx

    n = size(xi,1)

    do i=n-1,1,-1
      if(x.ge.xi(i)) then
        exit
      end if
    end do

    dx = x - xi(i)
    output = splineA(i) + dx*(splineB(i) + dx*(splineC(i) + splineD(i)*dx))
    computeNaturalCubicSpline = output
  end function computeNaturalCubicSpline

end module cubicSpline
