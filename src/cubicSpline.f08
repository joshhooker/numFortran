!*********************************************************
!>
! Spline Functions:
!
!  * Natural Cubic Spline

module cubicSpline
  use matrixSolver
  implicit none

  private :: dp
  integer, parameter :: dp = kind(0.d0)

contains

  subroutine naturalCubicSpline(xi, yi, splineA, splineB, splineC, splineD)
    implicit none
    integer :: n
    real(dp) :: xi(:), yi(:), splineA(:), splineB(:), splineC(:), splineD(:)
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
    h(1:n) = abs(xi(2:n+1)-xi(1:n))
    a(2:n) = h(1:n-1)
    b(2:n) = 2.d0*(h(1:n-1)+h(2:n))
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
    implicit none
    integer :: i, n
    real(dp) :: xi(:), x, splineA(:), splineB(:), splineC(:), splineD(:), output, dx

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

end module
