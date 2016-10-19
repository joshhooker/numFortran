module nfSpecialFunctions
  use nfConstants
  implicit none

  public :: legendrePoly

  interface legendrePoly
    module procedure legendrePoly_i, legendrePoly_r, legendrePoly_d
  end interface

contains

  recursive function legendrePolyFunc(n,x) result(res)
    integer :: i, n
    real(dp), allocatable :: poly(:)
    real(dp) :: x, res

    allocate(poly(n+1))
    poly(1) = 1.d0
    poly(2) = x
    if(n.eq.0) then
      res = 1
    elseif(n.eq.1) then
      res = x
    else
      do i=2,n
        res = 0.d0
        res = (2*(i-1)+1)*x*poly(i)
        res = res - (i-1)*poly(i-1)
        res = res/(i)
        poly(i+1) = res
      end do
      res = poly(n+1)
    end if
    deallocate(poly)
  end function

  real(dp) function legendrePoly_i(n,x)
    integer :: n, x
    legendrePoly_i = legendrePolyFunc(n,dble(x))
  end function

  real(dp) function legendrePoly_r(n,x)
    integer :: n
    real :: x
    legendrePoly_r = legendrePolyFunc(n,dble(x))
  end function

  real(dp) function legendrePoly_d(n,x)
    integer :: n
    real(dp) :: x
    legendrePoly_d = legendrePolyFunc(n,dble(x))
  end function

  end module nfSpecialFunctions