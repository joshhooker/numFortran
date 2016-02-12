module matrixSolver
  implicit none

  private
  integer, parameter :: dp = kind(0.d0)

  public :: tridiagonalSolver

contains

  subroutine tridiagonalSolver(a,b,c,f,u)
    implicit none
    integer :: i, n
    real(dp) :: a(:), b(:), c(:), f(:), u(:), bp
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