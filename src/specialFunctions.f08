!*********************************************************
!>
! Special Functions:
!
!  * Bessel Functions - J and Y
!  * Bessel Function - I and K
!  * Gamma Function (Not 100%)
!
! Future:
!  * Spherical Bessel Functions
!  * Airy Function - Ai and Bi
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

  public :: besselJ, besselY, besselI, besselK
  public :: gammaFunc

  interface besselJ
    module procedure besselJ_ii, besselJ_ir, besselJ_id, besselJ_ri, &
        besselJ_di, besselJ_rr, besselJ_dd
  end interface

  interface besselY
    module procedure besselY_ii, besselY_ir, besselY_id, besselY_ri, &
        besselY_di, besselY_rr, besselY_dd
  end interface

  interface besselI
    module procedure besselI_ii, besselI_ir, besselI_id, besselI_ri, &
        besselI_di, besselI_rr, besselI_dd
  end interface

  interface besselK
    module procedure besselK_ii, besselK_ir, besselK_id, besselK_ri, &
        besselK_di, besselK_rr, besselK_dd
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
    func = (exp(c(1)*(-log(x)))+((-1.d0)**c(1))*exp(-c(1)*(-log(x))))
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

  ! real(dp) function airyAFunction(n,c,x)
  !   integer :: n
  !   real(dp) :: c(n)
  !   real(dp) :: x
  !   real(dp) :: func
  !   func = cos(-(log(x)*log(x)*log(x))/3.d0-c(1)*log(x))
  !   airyAFunction = func
  !   !func = cos(x*x*x/3.d0+c(1)*x)
  !   !airyAFunction = func
  ! end function airyAFunction

  ! real(dp) function airyA(x)
  !   real(dp) :: x
  !   real(dp) :: intResults
  !   real(dp) :: consts(1)
  !   consts(1) = x
  !   intResults = gaussLegendre(airyAFunction,0.d0,1.d0,1,consts)
  !   airyA = intResults/pi_
  ! end function airyA

  ! real(dp) function airyBFunction(n,c,x)
  !   integer :: n
  !   real(dp) :: c(n)
  !   real(dp) :: x
  !   real(dp) :: func
  !   func = exp(-x*x*x/3.d0+c(1)*x)+sin(x*x*x/3.d0+c(1)*x)
  !   airyBFunction = func
  ! end function airyBFunction

  ! real(dp) function airyB(x)
  !   real(dp) :: x
  !   real(dp) :: intResults
  !   real(dp) :: consts(1)
  !   consts(1) = x
  !   intResults = gaussLegendre(airyBFunction,0.d0,1.d0,1,consts)
  !   airyB = intResults/pi_
  ! end function airyB

  real(dp) function gammaFunction(n,c,x)
    integer :: n
    real(dp) :: c(n)
    real(dp) :: x
    real(dp) :: func
    func = log(1.d0/x)**(c(1)-1.d0)
    gammaFunction = func
  end function gammaFunction

  real(dp) function gammaFunc(x)
    real(dp) :: x
    real(dp) :: intResults
    real(dp) :: consts(1)
    consts(1) = x
    intResults = gaussLegendre(gammaFunction,0.d0,1.d0,1,consts)
    gammaFunc = intResults
  end function gammaFunc


end module specialFunctions
