!*********************************************************
!>
! Special Functions:
!
!  * Bessel Functions - J and Y
!  * Bessel Function - I and K
!  * Spherical Bessel Functions
!  * Airy Function - Ai and Bi
!  * Gamma Function (Not 100%)
!
! Future:
!  * Legendre Function
!  * Hermite Polynomial
!  * Laguerre Polynomial
!  * Chevyshev Polynomial
!  * Beta Function
!  * Error Function?

module specialFunctions
  use constants
  use integration
  implicit none

  private
  integer, parameter :: dp = kind(0.d0)

  public :: besselJ, besselY, besselI, besselK, besselY_dd
  public :: sphBesselJ, sphBesselY
  public :: airyA, airyB
  public :: gammaFunc

  interface besselJ
    module procedure besselJ_ii, besselJ_ir, besselJ_id, besselJ_ri, &
        besselJ_di, besselJ_rr, besselJ_rd, besselJ_dr, besselJ_dd
  end interface

  interface besselY
    module procedure besselY_ii, besselY_ir, besselY_id, besselY_ri, &
        besselY_di, besselY_rr, besselY_rd, besselY_dr, besselY_dd
  end interface

  interface besselI
    module procedure besselI_ii, besselI_ir, besselI_id, besselI_ri, &
        besselI_di, besselI_rr, besselI_rd, besselI_dr, besselI_dd
  end interface

  interface besselK
    module procedure besselK_ii, besselK_ir, besselK_id, besselK_ri, &
        besselK_di, besselK_rr, besselK_rd, besselK_dr, besselK_dd
  end interface

  interface sphBesselJ
    module procedure sphBesselJ_ii, sphBesselJ_ir, sphBesselJ_id, sphBesselJ_ri, &
        sphBesselJ_di, sphBesselJ_rr, sphBesselJ_rd, sphBesselJ_dr, sphBesselJ_dd
  end interface

  interface sphBesselY
    module procedure sphBesselY_ii, sphBesselY_ir, sphBesselY_id, sphBesselY_ri, &
        sphBesselY_di, sphBesselY_rr, sphBesselY_rd, sphBesselY_dr, sphBesselY_dd
  end interface

  interface airyA
    module procedure airyA_i, airyA_r, airyA_d
  end interface

  interface airyB
    module procedure airyB_i, airyB_r, airyB_d
  end interface


contains

  real(dp) function besselJFunction1(n,c,x)
    integer :: n
    real(dp) :: c(n)
    real(dp) :: x
    besselJFunction1 = cos(c(1)*x-c(2)*sin(x))
  end function besselJFunction1

  real(dp) function besselJFunction2(n,c,x)
    integer :: n
    real(dp) :: c(n)
    real(dp) :: x
    real(dp) :: func
    func = exp(-c(2)*sinh(-log(x))+c(1)*log(x))
    besselJFunction2 = func*sin(c(1)*pi_)/x
  end function besselJFunction2

  real(dp) function besselJ_ii(n,x)
    !! Bessel function of the first kind \( J_{n}(x) \)
    integer :: n
      !! input: integer for Bessel function of the form /( J_{n}(x) /)
    integer :: x
      !! input: integer value at where Bessel function to be evaluated at
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselJFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselJFunction2,0.d0,1.d0,2,consts)
    besselJ_ii = intResults/pi_
  end function besselJ_ii

  real(dp) function besselJ_ir(n,x)
    !! Bessel function of the first kind \( J_{n}(x) \)
    integer :: n
      !! input: integer for Bessel function of the form /( J_{n}(x) /)
    real :: x
      !! input: value at where Bessel function to be evaluated at
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselJFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselJFunction2,0.d0,1.d0,2,consts)
    besselJ_ir = intResults/pi_
  end function besselJ_ir

  real(dp) function besselJ_id(n,x)
    !! Bessel function of the first kind \( J_{n}(x) \)
    integer :: n
      !! input: integer for Bessel function of the form /( J_{n}(x) /)
    real(dp) :: x
      !! input: value at where Bessel function to be evaluated at
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselJFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselJFunction2,0.d0,1.d0,2,consts)
    besselJ_id = intResults/pi_
  end function besselJ_id

  real(dp) function besselJ_ri(n,x)
    !! Bessel function of the first kind \( J_{n}(x) \)
    real :: n
      !! input: integer for Bessel function of the form /( J_{n}(x) /)
    integer :: x
      !! input: integer value at where Bessel function to be evaluated at
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselJFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselJFunction2,0.d0,1.d0,2,consts)
    besselJ_ri = intResults/pi_
  end function besselJ_ri

  real(dp) function besselJ_di(n,x)
    !! Bessel function of the first kind \( J_{n}(x) \)
    real(dp) :: n
      !! input: integer for Bessel function of the form /( J_{n}(x) /)
    integer :: x
      !! input: integer value at where Bessel function to be evaluated at
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselJFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselJFunction2,0.d0,1.d0,2,consts)
    besselJ_di = intResults/pi_
  end function besselJ_di

  real(dp) function besselJ_rr(n,x)
    !! Bessel function of the first kind \( J_{n}(x) \)
    real :: n
      !! input: real for Bessel function of the form /( J_{n}(x) /)
    real :: x
      !! input: real value at where Bessel function to be evaluated at
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselJFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselJFunction2,0.d0,1.d0,2,consts)
    besselJ_rr = intResults/pi_
  end function besselJ_rr

  real(dp) function besselJ_rd(n,x)
    !! Bessel function of the first kind \( J_{n}(x) \)
    real :: n
      !! input: real for Bessel function of the form /( J_{n}(x) /)
    real(dp) :: x
      !! input: real value at where Bessel function to be evaluated at
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselJFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselJFunction2,0.d0,1.d0,2,consts)
    besselJ_rd = intResults/pi_
  end function besselJ_rd

  real(dp) function besselJ_dr(n,x)
    !! Bessel function of the first kind \( J_{n}(x) \)
    real(dp) :: n
      !! input: real for Bessel function of the form /( J_{n}(x) /)
    real :: x
      !! input: real value at where Bessel function to be evaluated at
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselJFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselJFunction2,0.d0,1.d0,2,consts)
    besselJ_dr = intResults/pi_
  end function besselJ_dr

  real(dp) function besselJ_dd(n,x)
    !! Bessel function of the first kind \( J_{n}(x) \)
    real(dp) :: n
      !! input: double for Bessel function of the form /( J_{n}(x) /)
    real(dp) :: x
      !! input: double value at where Bessel function to be evaluated at
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselJFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselJFunction2,0.d0,1.d0,2,consts)
    besselJ_dd = intResults/pi_
  end function besselJ_dd

  real(dp) function besselYFunction1(n,c,x)
    integer :: n
    real(dp) :: c(n)
    real(dp) :: x
    besselYFunction1 = sin(c(2)*sin(x)-c(1)*x)
  end function besselYFunction1

  real(dp) function besselYFunction2(n,c,x)
    integer :: n
    real(dp) :: c(n)
    real(dp) :: x
    real(dp) :: func
    func = 0.d0
    func = (exp(c(1)*(-log(x)))+(cos(c(1)*pi_))*exp(-c(1)*(-log(x))))
    func = func*exp(-c(2)*sinh((-log(x))))
    besselYFunction2 = func/x
  end function besselYFunction2

  real(dp) function besselY_ii(n,x)
    integer :: n
    integer :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselYFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselYFunction2,0.d0,1.d0,2,consts)
    besselY_ii = intResults/pi_
  end function besselY_ii

  real(dp) function besselY_ir(n,x)
    integer :: n
    real :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselYFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselYFunction2,0.d0,1.d0,2,consts)
    besselY_ir = intResults/pi_
  end function besselY_ir

  real(dp) function besselY_id(n,x)
    integer :: n
    real(dp) :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselYFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselYFunction2,0.d0,1.d0,2,consts)
    besselY_id = intResults/pi_
  end function besselY_id

  real(dp) function besselY_ri(n,x)
    real :: n
    integer :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselYFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselYFunction2,0.d0,1.d0,2,consts)
    besselY_ri = intResults/pi_
  end function besselY_ri

  real(dp) function besselY_di(n,x)
    real(dp) :: n
    integer :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselYFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselYFunction2,0.d0,1.d0,2,consts)
    besselY_di = intResults/pi_
  end function besselY_di

  real(dp) function besselY_rr(n,x)
    real :: n
    real :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselYFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselYFunction2,0.d0,1.d0,2,consts)
    besselY_rr = intResults/pi_
  end function besselY_rr

  real(dp) function besselY_rd(n,x)
    real :: n
    real(dp) :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselYFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselYFunction2,0.d0,1.d0,2,consts)
    besselY_rd = intResults/pi_
  end function besselY_rd

  real(dp) function besselY_dr(n,x)
    real(dp) :: n
    real :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselYFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselYFunction2,0.d0,1.d0,2,consts)
    besselY_dr = intResults/pi_
  end function besselY_dr

  real(dp) function besselY_dd(n,x)
    real(dp) :: n
    real(dp) :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselYFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselYFunction2,0.d0,1.d0,2,consts)
    besselY_dd = intResults/pi_
  end function besselY_dd

  real(dp) function besselIFunction1(n,c,x)
    integer :: n
    real(dp) :: c(n)
    real(dp) :: x
    real(dp) :: func
    func = exp(c(2)*cos(x))*cos(c(1)*x)
    besselIFunction1 = func
  end function besselIFunction1

  real(dp) function besselIFunction2(n,c,x)
    integer :: n
    real(dp) :: c(n)
    real(dp) :: x
    real(dp) :: func
    func = exp(-c(2)*cosh(-log(x))+c(1)*log(x))
    besselIFunction2 = func*sin(c(1)*pi_)/x
  end function besselIFunction2

  real(dp) function besselI_ii(n,x)
    integer :: n
    integer :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselIFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselIFunction2,0.d0,1.d0,2,consts)
    besselI_ii = intResults/pi_
  end function besselI_ii

  real(dp) function besselI_ir(n,x)
    integer :: n
    real :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselIFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselIFunction2,0.d0,1.d0,2,consts)
    besselI_ir = intResults/pi_
  end function besselI_ir

  real(dp) function besselI_id(n,x)
    integer :: n
    real(dp) :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselIFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselIFunction2,0.d0,1.d0,2,consts)
    besselI_id = intResults/pi_
  end function besselI_id

  real(dp) function besselI_ri(n,x)
    real :: n
    integer :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselIFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselIFunction2,0.d0,1.d0,2,consts)
    besselI_ri = intResults/pi_
  end function besselI_ri

  real(dp) function besselI_di(n,x)
    real(dp) :: n
    integer :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselIFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselIFunction2,0.d0,1.d0,2,consts)
    besselI_di = intResults/pi_
  end function besselI_di

  real(dp) function besselI_rr(n,x)
    real :: n
    real :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselIFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselIFunction2,0.d0,1.d0,2,consts)
    besselI_rr = intResults/pi_
  end function besselI_rr

  real(dp) function besselI_rd(n,x)
    real :: n
    real(dp) :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselIFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselIFunction2,0.d0,1.d0,2,consts)
    besselI_rd = intResults/pi_
  end function besselI_rd

  real(dp) function besselI_dr(n,x)
    real(dp) :: n
    real :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselIFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselIFunction2,0.d0,1.d0,2,consts)
    besselI_dr = intResults/pi_
  end function besselI_dr

  real(dp) function besselI_dd(n,x)
    real(dp) :: n
    real(dp) :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselIFunction1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselIFunction2,0.d0,1.d0,2,consts)
    besselI_dd = intResults/pi_
  end function besselI_dd

  real(dp) function besselKFunction(n,c,x)
    integer :: n
    real(dp) :: c(n)
    real(dp) :: x
    real(dp) :: func
    func = exp(-c(2)*cosh(-log(x)))*cosh(-c(1)*log(x))
    besselKFunction = func/x
  end function besselKFunction

  real(dp) function besselK_ii(n,x)
    integer :: n
    integer :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselKFunction,0.d0,1.d0,2,consts)
    besselK_ii = intResults
  end function besselK_ii

  real(dp) function besselK_ir(n,x)
    integer :: n
    real :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselKFunction,0.d0,1.d0,2,consts)
    besselK_ir = intResults
  end function besselK_ir

  real(dp) function besselK_id(n,x)
    integer :: n
    real(dp) :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselKFunction,0.d0,1.d0,2,consts)
    besselK_id = intResults
  end function besselK_id

  real(dp) function besselK_ri(n,x)
    real :: n
    integer :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselKFunction,0.d0,1.d0,2,consts)
    besselK_ri = intResults
  end function besselK_ri

  real(dp) function besselK_di(n,x)
    real(dp) :: n
    integer :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselKFunction,0.d0,1.d0,2,consts)
    besselK_di = intResults
  end function besselK_di

  real(dp) function besselK_rr(n,x)
    real :: n
    real :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselKFunction,0.d0,1.d0,2,consts)
    besselK_rr = intResults
  end function besselK_rr

  real(dp) function besselK_rd(n,x)
    real :: n
    real(dp) :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselKFunction,0.d0,1.d0,2,consts)
    besselK_rd = intResults
  end function besselK_rd

  real(dp) function besselK_dr(n,x)
    real(dp) :: n
    real :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselKFunction,0.d0,1.d0,2,consts)
    besselK_dr = intResults
  end function besselK_dr

  real(dp) function besselK_dd(n,x)
    real(dp) :: n
    real(dp) :: x
    real(dp) :: intResults
    real(dp) :: consts(2)
    consts(1) = dble(n)
    consts(2) = dble(x)
    intResults = gaussLegendre(besselKFunction,0.d0,1.d0,2,consts)
    besselK_dd = intResults
  end function besselK_dd

  real(dp) function sphBesselJ_ii(n,x)
    integer :: n
    integer :: x
    if(x.eq.0) then
      if(n.eq.0) then
        sphBesselJ_ii = 1.d0
      else
        sphBesselJ_ii = 0.d0
      end if
    else
      sphBesselJ_ii = sqrt(pi_/(2.d0*dble(x)))*besselJ(dble(n)+0.5d0,dble(x))
    end if
  end function sphBesselJ_ii

  real(dp) function sphBesselJ_ir(n,x)
    integer :: n
    real :: x
    if(x.eq.0.0) then
      if(n.eq.0) then
        sphBesselJ_ir = 1.d0
      else
        sphBesselJ_ir = 0.d0
      end if
    else
      sphBesselJ_ir = sqrt(pi_/(2.d0*dble(x)))*besselJ(dble(n)+0.5d0,dble(x))
    end if
  end function sphBesselJ_ir

  real(dp) function sphBesselJ_id(n,x)
    integer :: n
    real(dp) :: x
    if(x.eq.0.d0) then
      if(n.eq.0) then
        sphBesselJ_id = 1.d0
      else
        sphBesselJ_id = 0.d0
      end if
    else
      sphBesselJ_id = sqrt(pi_/(2.d0*dble(x)))*besselJ(dble(n)+0.5d0,dble(x))
    end if
  end function sphBesselJ_id

  real(dp) function sphBesselJ_ri(n,x)
    real :: n
    integer :: x
    if(x.eq.0) then
      if(n.eq.0.0) then
        sphBesselJ_ri = 1.d0
      else
        sphBesselJ_ri = 0.d0
      end if
    else
      sphBesselJ_ri = sqrt(pi_/(2.d0*dble(x)))*besselJ(dble(n)+0.5d0,dble(x))
    end if
  end function sphBesselJ_ri

  real(dp) function sphBesselJ_di(n,x)
    real(dp) :: n
    integer :: x
    if(x.eq.0) then
      if(n.eq.0.d0) then
        sphBesselJ_di = 1.d0
      else
        sphBesselJ_di = 0.d0
      end if
    else
      sphBesselJ_di = sqrt(pi_/(2.d0*dble(x)))*besselJ(dble(n)+0.5d0,dble(x))
    end if
  end function sphBesselJ_di

  real(dp) function sphBesselJ_rr(n,x)
    real :: n
    real :: x
    if(x.eq.0.0) then
      if(n.eq.0.0) then
        sphBesselJ_rr = 1.d0
      else
        sphBesselJ_rr = 0.d0
      end if
    else
      sphBesselJ_rr = sqrt(pi_/(2.d0*dble(x)))*besselJ(dble(n)+0.5d0,dble(x))
    end if
  end function sphBesselJ_rr

  real(dp) function sphBesselJ_rd(n,x)
    real :: n
    real(dp) :: x
    if(x.eq.0.d0) then
      if(n.eq.0.0) then
        sphBesselJ_rd = 1.d0
      else
        sphBesselJ_rd = 0.d0
      end if
    else
      sphBesselJ_rd = sqrt(pi_/(2.d0*dble(x)))*besselJ(dble(n)+0.5d0,dble(x))
    end if
  end function sphBesselJ_rd

  real(dp) function sphBesselJ_dr(n,x)
    real(dp) :: n
    real :: x
    if(x.eq.0.0) then
      if(n.eq.0.d0) then
        sphBesselJ_dr = 1.d0
      else
        sphBesselJ_dr = 0.d0
      end if
    else
      sphBesselJ_dr = sqrt(pi_/(2.d0*dble(x)))*besselJ(dble(n)+0.5d0,dble(x))
    end if
  end function sphBesselJ_dr

  real(dp) function sphBesselJ_dd(n,x)
    real(dp) :: n
    real(dp) :: x
    if(x.eq.0.d0) then
      if(n.eq.0.d0) then
        sphBesselJ_dd = 1.d0
      else
        sphBesselJ_dd = 0.d0
      end if
    else
      sphBesselJ_dd = sqrt(pi_/(2.d0*dble(x)))*besselJ(dble(n)+0.5d0,dble(x))
    end if
  end function sphBesselJ_dd

  real(dp) function sphBesselY_ii(n,x)
    integer :: n
    integer :: x
    sphBesselY_ii = sqrt(pi_/(2.d0*dble(x)))*besselY(dble(n)+0.5d0,dble(x))
  end function sphBesselY_ii

  real(dp) function sphBesselY_ir(n,x)
    integer :: n
    real :: x
    sphBesselY_ir = sqrt(pi_/(2.d0*dble(x)))*besselY(dble(n)+0.5d0,dble(x))
  end function sphBesselY_ir

  real(dp) function sphBesselY_id(n,x)
    integer :: n
    real(dp) :: x
    sphBesselY_id = sqrt(pi_/(2.d0*dble(x)))*besselY(dble(n)+0.5d0,dble(x))
  end function sphBesselY_id

  real(dp) function sphBesselY_ri(n,x)
    real :: n
    integer :: x
    sphBesselY_ri = sqrt(pi_/(2.d0*dble(x)))*besselY(dble(n)+0.5d0,dble(x))
  end function sphBesselY_ri

  real(dp) function sphBesselY_di(n,x)
    real(dp) :: n
    integer :: x
    sphBesselY_di = sqrt(pi_/(2.d0*dble(x)))*besselY(dble(n)+0.5d0,dble(x))
  end function sphBesselY_di

  real(dp) function sphBesselY_rr(n,x)
    real :: n
    real :: x
    sphBesselY_rr = sqrt(pi_/(2.d0*dble(x)))*besselY(dble(n)+0.5d0,dble(x))
  end function sphBesselY_rr

  real(dp) function sphBesselY_rd(n,x)
    real :: n
    real(dp) :: x
    sphBesselY_rd = sqrt(pi_/(2.d0*dble(x)))*besselY(dble(n)+0.5d0,dble(x))
  end function sphBesselY_rd

  real(dp) function sphBesselY_dr(n,x)
    real(dp) :: n
    real :: x
    sphBesselY_dr = sqrt(pi_/(2.d0*dble(x)))*besselY(dble(n)+0.5d0,dble(x))
  end function sphBesselY_dr

  real(dp) function sphBesselY_dd(n,x)
    real(dp) :: n
    real(dp) :: x
    sphBesselY_dd = sqrt(pi_/(2.d0*dble(x)))*besselY(dble(n)+0.5d0,dble(x))
  end function sphBesselY_dd

  real(dp) function airyA_i(x)
    integer :: x
    real(dp) :: airyResult
    airyResult = 0.d0
    if(dble(x).lt.0.d0) then
      airyResult = besselJ(1.d0/3.d0,(2.d0/3.d0)*(abs(dble(x))**(3.d0/2.d0)))
      airyResult = airyResult+besselJ(-1.d0/3.d0,(2.d0/3.d0)*(abs(dble(x))**(3.d0/2.d0)))
      airyResult = airyResult*sqrt(abs(dble(x))/9.d0)
    elseif(dble(x).eq.0.d0) then
      airyResult = 1.d0/((3.d0**(2.d0/3.d0))*gamma(2.d0/3.d0))
    else
      airyResult = besselK(1.d0/3.d0,(2.d0/3.d0)*(dble(x)**(3.d0/2.d0)))
      airyResult = airyResult*sqrt(dble(x)/3.d0)/pi_
    end if
    airyA_i = airyResult
  end function airyA_i

  real(dp) function airyA_r(x)
    real :: x
    real(dp) :: airyResult
    airyResult = 0.d0
    if(dble(x).lt.0.d0) then
      airyResult = besselJ(1.d0/3.d0,(2.d0/3.d0)*(abs(dble(x))**(3.d0/2.d0)))
      airyResult = airyResult+besselJ(-1.d0/3.d0,(2.d0/3.d0)*(abs(dble(x))**(3.d0/2.d0)))
      airyResult = airyResult*sqrt(abs(dble(x))/9.d0)
    elseif(dble(x).eq.0.d0) then
      airyResult = 1.d0/((3.d0**(2.d0/3.d0))*gamma(2.d0/3.d0))
    else
      airyResult = besselK(1.d0/3.d0,(2.d0/3.d0)*(dble(x)**(3.d0/2.d0)))
      airyResult = airyResult*sqrt(dble(x)/3.d0)/pi_
    end if
    airyA_r = airyResult
  end function airyA_r

  real(dp) function airyA_d(x)
    real(dp) :: x
    real(dp) :: airyResult
    airyResult = 0.d0
    if(dble(x).lt.0.d0) then
      airyResult = besselJ(1.d0/3.d0,(2.d0/3.d0)*(abs(dble(x))**(3.d0/2.d0)))
      airyResult = airyResult+besselJ(-1.d0/3.d0,(2.d0/3.d0)*(abs(dble(x))**(3.d0/2.d0)))
      airyResult = airyResult*sqrt(abs(dble(x))/9.d0)
    elseif(dble(x).eq.0.d0) then
      airyResult = 1.d0/((3.d0**(2.d0/3.d0))*gamma(2.d0/3.d0))
    else
      airyResult = besselK(1.d0/3.d0,(2.d0/3.d0)*(dble(x)**(3.d0/2.d0)))
      airyResult = airyResult*sqrt(dble(x)/3.d0)/pi_
    end if
    airyA_d = airyResult
  end function airyA_d

  real(dp) function airyB_i(x)
    integer :: x
    real(dp) :: airyResult
    airyResult = 0.d0
    if(dble(x).lt.0.d0) then
      airyResult = besselJ(-1.d0/3.d0,(2.d0/3.d0)*(abs(dble(x))**(3.d0/2.d0)))
      airyResult = airyResult-besselJ(1.d0/3.d0,(2.d0/3.d0)*(abs(dble(x))**(3.d0/2.d0)))
      airyResult = airyResult*sqrt(abs(dble(x))/3.d0)
    elseif(dble(x).eq.0.d0) then
      airyResult = 1.d0/((3.d0**(1.d0/6.d0))*gamma(2.d0/3.d0))
    else
      airyResult = besselI(1.d0/3.d0,(2.d0/3.d0)*(dble(x)**(3.d0/2.d0)))
      airyResult = airyResult+besselI(-1.d0/3.d0,(2.d0/3.d0)*(dble(x)**(3.d0/2.d0)))
      airyResult = airyResult*sqrt(dble(x)/3.d0)
    end if
    airyB_i = airyResult
  end function airyB_i

  real(dp) function airyB_r(x)
    real :: x
    real(dp) :: airyResult
    airyResult = 0.d0
    if(dble(x).lt.0.d0) then
      airyResult = besselJ(-1.d0/3.d0,(2.d0/3.d0)*(abs(dble(x))**(3.d0/2.d0)))
      airyResult = airyResult-besselJ(1.d0/3.d0,(2.d0/3.d0)*(abs(dble(x))**(3.d0/2.d0)))
      airyResult = airyResult*sqrt(abs(dble(x))/3.d0)
    elseif(dble(x).eq.0.d0) then
      airyResult = 1.d0/((3.d0**(1.d0/6.d0))*gamma(2.d0/3.d0))
    else
      airyResult = besselI(1.d0/3.d0,(2.d0/3.d0)*(dble(x)**(3.d0/2.d0)))
      airyResult = airyResult+besselI(-1.d0/3.d0,(2.d0/3.d0)*(dble(x)**(3.d0/2.d0)))
      airyResult = airyResult*sqrt(dble(x)/3.d0)
    end if
    airyB_r = airyResult
  end function airyB_r

  real(dp) function airyB_d(x)
    real(dp) :: x
    real(dp) :: airyResult
    airyResult = 0.d0
    if(dble(x).lt.0.d0) then
      airyResult = besselJ(-1.d0/3.d0,(2.d0/3.d0)*(abs(dble(x))**(3.d0/2.d0)))
      airyResult = airyResult-besselJ(1.d0/3.d0,(2.d0/3.d0)*(abs(dble(x))**(3.d0/2.d0)))
      airyResult = airyResult*sqrt(abs(dble(x))/3.d0)
    elseif(dble(x).eq.0.d0) then
      airyResult = 1.d0/((3.d0**(1.d0/6.d0))*gamma(2.d0/3.d0))
    else
      airyResult = besselI(1.d0/3.d0,(2.d0/3.d0)*(dble(x)**(3.d0/2.d0)))
      airyResult = airyResult+besselI(-1.d0/3.d0,(2.d0/3.d0)*(dble(x)**(3.d0/2.d0)))
      airyResult = airyResult*sqrt(dble(x)/3.d0)
    end if
    airyB_d = airyResult
  end function airyB_d

  real(dp) function gammaFunction(n,c,x)
    integer :: n
    real(dp) :: c(n)
    real(dp) :: x
    real(dp) :: func
    func = (-log(x))**(c(1)-1.d0)
    gammaFunction = func
  end function gammaFunction

  real(dp) function gammaFunc(x)
    real(dp) :: x
    real(dp) :: intResults
    real(dp) :: consts(1)
    consts(1) = x
    intResults = gaussLegendre(gammaFunction,0.d0,0.01d0,1,consts)
    intResults = intResults+gaussLegendre(gammaFunction,0.01d0,0.02d0,1,consts)
    intResults = intResults+gaussLegendre(gammaFunction,0.02d0,0.03d0,1,consts)
    intResults = intResults+gaussLegendre(gammaFunction,0.03d0,0.04d0,1,consts)
    intResults = intResults+gaussLegendre(gammaFunction,0.04d0,0.05d0,1,consts)
    intResults = intResults+gaussLegendre(gammaFunction,0.05d0,0.06d0,1,consts)
    intResults = intResults+gaussLegendre(gammaFunction,0.06d0,0.07d0,1,consts)
    intResults = intResults+gaussLegendre(gammaFunction,0.07d0,0.08d0,1,consts)
    intResults = intResults+gaussLegendre(gammaFunction,0.08d0,0.09d0,1,consts)
    intResults = intResults+gaussLegendre(gammaFunction,0.09d0,0.1d0,1,consts)
    intResults = intResults+gaussLegendre(gammaFunction,0.1d0,0.2d0,1,consts)
    intResults = intResults+gaussLegendre(gammaFunction,0.2d0,0.3d0,1,consts)
    intResults = intResults+gaussLegendre(gammaFunction,0.3d0,0.4d0,1,consts)
    intResults = intResults+gaussLegendre(gammaFunction,0.4d0,0.5d0,1,consts)
    intResults = intResults+gaussLegendre(gammaFunction,0.5d0,0.6d0,1,consts)
    intResults = intResults+gaussLegendre(gammaFunction,0.6d0,0.7d0,1,consts)
    intResults = intResults+gaussLegendre(gammaFunction,0.7d0,0.8d0,1,consts)
    intResults = intResults+gaussLegendre(gammaFunction,0.8d0,0.9d0,1,consts)
    intResults = intResults+gaussLegendre(gammaFunction,0.9d0,1.d0,1,consts)
    gammaFunc = intResults
  end function gammaFunc


end module specialFunctions
