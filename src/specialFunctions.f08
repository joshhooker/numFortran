!*********************************************************
!>
! Special Functions:
!
!  * Bessel Functions - J and Y
!  * Bessel Function - I and K
!  * Spherical Bessel Functions
!  * Airy Function - Ai and Bi
!  * Error Function and Complementary Error Function
!  * Beta Function
!  * Incomplete Beta Function
!
! Not 100% Right:
!  * Gamma Function
!
! Future:
!  * Legendre Polynomial
!  * Hermite Polynomial
!  * Laguerre Polynomial
!  * Chevyshev Polynomial

module specialFunctions
  use constants
  use integration
  implicit none

  private
  integer, parameter :: dp = kind(0.d0)

  public :: besselJ, besselY, besselI, besselK, besselY_dd
  public :: sphBesselJ, sphBesselY
  public :: airyA, airyB
  public :: errf, errfc
  public :: beta
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

  interface errf
    module procedure errf_i, errf_r, errf_d
  end interface

  interface errfc
    module procedure errfc_i, errfc_r, errfc_d
  end interface

  interface beta
    module procedure beta_ii, beta_ir, beta_id, beta_ri, beta_di, beta_rr, &
        beta_rd, beta_dr, beta_dd
  end interface

  interface incBeta
    module procedure incBeta_iii, incBeta_iir, incBeta_iid, incBeta_iri, incBeta_irr, &
        incBeta_ird, incBeta_idi, incBeta_idr, incBeta_idd, incBeta_rii, incBeta_rir, &
        incBeta_rid, incBeta_rri, incBeta_rrr, incBeta_rrd, incBeta_rdi, incBeta_rdr, &
        incBeta_rdd, incBeta_dii, incBeta_dir, incBeta_did, incBeta_dri, incBeta_drr, &
        incBeta_drd, incBeta_ddi, incBeta_ddr, incBeta_ddd
  end interface

contains

  !***********************************!
  ! Bessel Function of the First Kind !
  !***********************************!

  real(dp) function besselJ_IntFunc1(n,c,x)
    integer :: n
    real(dp) :: c(n),x
    besselJ_IntFunc1 = cos(c(1)*x-c(2)*sin(x))
  end function

  real(dp) function besselJ_IntFunc2(n,c,x)
    integer :: n
    real(dp) :: c(n), x, func
    func = exp(-c(2)*sinh(-log(x))+c(1)*log(x))
    besselJ_IntFunc2 = func*sin(c(1)*pi_)/x
  end function

  real(dp) function besselJFunc(n,x)
    real(dp) :: n, x, intResults, consts(2)
    consts(1) = n
    consts(2) = x
    intResults = gaussLegendre(besselJ_IntFunc1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselJ_IntFunc2,0.d0,1.d0,2,consts)
    besselJFunc = intResults/pi_
  end function

  real(dp) function besselJ_ii(n,x)
    integer :: n, x
    besselJ_ii = besselJFunc(dble(n),dble(x))
  end function

  real(dp) function besselJ_ir(n,x)
    integer :: n
    real :: x
    besselJ_ir = besselJFunc(dble(n),dble(x))
  end function

  real(dp) function besselJ_id(n,x)
    integer :: n
    real(dp) :: x
    besselJ_id = besselJFunc(dble(n),dble(x))
  end function

  real(dp) function besselJ_ri(n,x)
    real :: n
    integer :: x
    besselJ_ri = besselJFunc(dble(n),dble(x))
  end function

  real(dp) function besselJ_di(n,x)
    real(dp) :: n
    integer :: x
    besselJ_di = besselJFunc(dble(n),dble(x))
  end function

  real(dp) function besselJ_rr(n,x)
    real :: n, x
    besselJ_rr = besselJFunc(dble(n),dble(x))
  end function

  real(dp) function besselJ_rd(n,x)
    real :: n
    real(dp) :: x
    besselJ_rd = besselJFunc(dble(n),dble(x))
  end function

  real(dp) function besselJ_dr(n,x)
    real(dp) :: n
    real :: x
    besselJ_dr = besselJFunc(dble(n),dble(x))
  end function

  real(dp) function besselJ_dd(n,x)
    real(dp) :: n, x
    besselJ_dd = besselJFunc(dble(n),dble(x))
  end function

  !************************************!
  ! Bessel Function of the Second Kind !
  !************************************!

  real(dp) function besselY_IntFunc1(n,c,x)
    integer :: n
    real(dp) :: c(n), x
    besselY_IntFunc1 = sin(c(2)*sin(x)-c(1)*x)
  end function

  real(dp) function besselY_IntFunc2(n,c,x)
    integer :: n
    real(dp) :: c(n), x, func
    func = (exp(c(1)*(-log(x)))+(cos(c(1)*pi_))*exp(-c(1)*(-log(x))))
    func = func*exp(-c(2)*sinh((-log(x))))
    besselY_IntFunc2 = func/x
  end function

  real(dp) function besselYFunc(n,x)
    real(dp) :: n, x, intResults, consts(2)
    consts(1) = n
    consts(2) = x
    intResults = gaussLegendre(besselY_IntFunc1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselY_IntFunc2,0.d0,1.d0,2,consts)
    besselYFunc = intResults/pi_
  end function

  real(dp) function besselY_ii(n,x)
    integer :: n, x
    besselY_ii = besselYFunc(dble(n),dble(x))
  end function

  real(dp) function besselY_ir(n,x)
    integer :: n
    real :: x
    besselY_ir = besselYFunc(dble(n),dble(x))
  end function

  real(dp) function besselY_id(n,x)
    integer :: n
    real(dp) :: x
    besselY_id = besselYFunc(dble(n),dble(x))
  end function

  real(dp) function besselY_ri(n,x)
    real :: n
    integer :: x
    besselY_ri = besselYFunc(dble(n),dble(x))
  end function

  real(dp) function besselY_di(n,x)
    real(dp) :: n
    integer :: x
    besselY_di = besselYFunc(dble(n),dble(x))
  end function

  real(dp) function besselY_rr(n,x)
    real :: n, x
    besselY_rr = besselYFunc(dble(n),dble(x))
  end function

  real(dp) function besselY_rd(n,x)
    real :: n
    real(dp) :: x
    besselY_rd = besselYFunc(dble(n),dble(x))
  end function

  real(dp) function besselY_dr(n,x)
    real(dp) :: n
    real :: x
    besselY_dr = besselYFunc(dble(n),dble(x))
  end function

  real(dp) function besselY_dd(n,x)
    real(dp) :: n, x
    besselY_dd = besselYFunc(dble(n),dble(x))
  end function

  !********************************************!
  ! Modified Bessel Function of the First Kind !
  !********************************************!

  real(dp) function besselI_IntFunc1(n,c,x)
    integer :: n
    real(dp) :: c(n), x
    besselI_IntFunc1 = exp(c(2)*cos(x))*cos(c(1)*x)
  end function

  real(dp) function besselI_IntFunc2(n,c,x)
    integer :: n
    real(dp) :: c(n), x, func
    func = exp(-c(2)*cosh(-log(x))+c(1)*log(x))
    besselI_IntFunc2 = func*sin(c(1)*pi_)/x
  end function

  real(dp) function besselIFunc(n,x)
    real(dp) :: n, x, intResults, consts(2)
    consts(1) = n
    consts(2) = x
    intResults = gaussLegendre(besselI_IntFunc1,0.d0,pi_,2,consts)
    intResults = intResults- &
        gaussLegendre(besselI_IntFunc2,0.d0,1.d0,2,consts)
    besselIFunc = intResults/pi_
  end function

  real(dp) function besselI_ii(n,x)
    integer :: n, x
    besselI_ii = besselIFunc(dble(n),dble(x))
  end function

  real(dp) function besselI_ir(n,x)
    integer :: n
    real :: x
    besselI_ir = besselIFunc(dble(n),dble(x))
  end function

  real(dp) function besselI_id(n,x)
    integer :: n
    real(dp) :: x
    besselI_id = besselIFunc(dble(n),dble(x))
  end function

  real(dp) function besselI_ri(n,x)
    real :: n
    integer :: x
    besselI_ri = besselIFunc(dble(n),dble(x))
  end function

  real(dp) function besselI_di(n,x)
    real(dp) :: n
    integer :: x
    besselI_di = besselIFunc(dble(n),dble(x))
  end function

  real(dp) function besselI_rr(n,x)
    real :: n, x
    besselI_rr = besselIFunc(dble(n),dble(x))
  end function

  real(dp) function besselI_rd(n,x)
    real :: n
    real(dp) :: x
    besselI_rd = besselIFunc(dble(n),dble(x))
  end function

  real(dp) function besselI_dr(n,x)
    real(dp) :: n
    real :: x
    besselI_dr = besselIFunc(dble(n),dble(x))
  end function

  real(dp) function besselI_dd(n,x)
    real(dp) :: n, x
    besselI_dd = besselIFunc(dble(n),dble(x))
  end function

  !*********************************************!
  ! Modified Bessel Function of the Second Kind !
  !*********************************************!

  real(dp) function besselK_IntFunc(n,c,x)
    integer :: n
    real(dp) :: c(n),x, func
    besselK_IntFunc = (exp(-c(2)*cosh(-log(x)))*cosh(-c(1)*log(x)))/x
  end function besselK_IntFunc

  real(dp) function besselKFunc(n,x)
    real(dp) :: n, x, intResults, consts(2)
    consts(1) = n
    consts(2) = x
    besselKFunc = gaussLegendre(besselK_IntFunc,0.d0,1.d0,2,consts)
  end function

  real(dp) function besselK_ii(n,x)
    integer :: n, x
    besselK_ii = besselKFunc(dble(n),dble(x))
  end function

  real(dp) function besselK_ir(n,x)
    integer :: n
    real :: x
    besselK_ir = besselKFunc(dble(n),dble(x))
  end function

  real(dp) function besselK_id(n,x)
    integer :: n
    real(dp) :: x
    besselK_id = besselKFunc(dble(n),dble(x))
  end function

  real(dp) function besselK_ri(n,x)
    real :: n
    integer :: x
    besselK_ri = besselKFunc(dble(n),dble(x))
  end function

  real(dp) function besselK_di(n,x)
    real(dp) :: n
    integer :: x
    besselK_di = besselKFunc(dble(n),dble(x))
  end function

  real(dp) function besselK_rr(n,x)
    real :: n, x
    besselK_rr = besselKFunc(dble(n),dble(x))
  end function

  real(dp) function besselK_rd(n,x)
    real :: n
    real(dp) :: x
    besselK_rd = besselKFunc(dble(n),dble(x))
  end function

  real(dp) function besselK_dr(n,x)
    real(dp) :: n
    real :: x
    besselK_dr = besselKFunc(dble(n),dble(x))
  end function

  real(dp) function besselK_dd(n,x)
    real(dp) :: n, x
    besselK_dd = besselKFunc(dble(n),dble(x))
  end function

  !*********************************************!
  ! Spherical Bessel Function of the First Kind !
  !*********************************************!

  real(dp) function sphBesselJFunc(n,x)
    real(dp) :: n, x
    if(x.eq.0d0) then
      if(n.eq.0.d0) then
        sphBesselJFunc = 1.d0
      else
        sphBesselJFunc = 0.d0
      end if
    else
      sphBesselJFunc = sqrt(pi_/(2.d0*x))*besselJ(n+0.5d0,x)
    end if
  end function

  real(dp) function sphBesselJ_ii(n,x)
    integer :: n, x
    sphBesselJ_ii = sphBesselJFunc(dble(n),dble(x))
  end function

  real(dp) function sphBesselJ_ir(n,x)
    integer :: n
    real :: x
    sphBesselJ_ir = sphBesselJFunc(dble(n),dble(x))
  end function

  real(dp) function sphBesselJ_id(n,x)
    integer :: n
    real(dp) :: x
    sphBesselJ_id = sphBesselJFunc(dble(n),dble(x))
  end function

  real(dp) function sphBesselJ_ri(n,x)
    real :: n
    integer :: x
    sphBesselJ_ri = sphBesselJFunc(dble(n),dble(x))
  end function

  real(dp) function sphBesselJ_di(n,x)
    real(dp) :: n
    integer :: x
    sphBesselJ_di = sphBesselJFunc(dble(n),dble(x))
  end function

  real(dp) function sphBesselJ_rr(n,x)
    real :: n, x
    sphBesselJ_rr = sphBesselJFunc(dble(n),dble(x))
  end function

  real(dp) function sphBesselJ_rd(n,x)
    real :: n
    real(dp) :: x
    sphBesselJ_rd = sphBesselJFunc(dble(n),dble(x))
  end function

  real(dp) function sphBesselJ_dr(n,x)
    real(dp) :: n
    real :: x
    sphBesselJ_dr = sphBesselJFunc(dble(n),dble(x))
  end function

  real(dp) function sphBesselJ_dd(n,x)
    real(dp) :: n, x
    sphBesselJ_dd = sphBesselJFunc(dble(n),dble(x))
  end function

  !**********************************************!
  ! Spherical Bessel Function of the Second Kind !
  !**********************************************!

  real(dp) function sphBesselYFunc(n,x)
    real(dp) :: n, x
    sphBesselYFunc = sqrt(pi_/(2.d0*x))*besselY(n+0.5d0,x)
  end function

  real(dp) function sphBesselY_ii(n,x)
    integer :: n, x
    sphBesselY_ii = sphBesselYFunc(dble(n),dble(x))
  end function

  real(dp) function sphBesselY_ir(n,x)
    integer :: n
    real :: x
    sphBesselY_ir = sphBesselYFunc(dble(n),dble(x))
  end function

  real(dp) function sphBesselY_id(n,x)
    integer :: n
    real(dp) :: x
    sphBesselY_id = sphBesselYFunc(dble(n),dble(x))
  end function

  real(dp) function sphBesselY_ri(n,x)
    real :: n
    integer :: x
    sphBesselY_ri = sphBesselYFunc(dble(n),dble(x))
  end function

  real(dp) function sphBesselY_di(n,x)
    real(dp) :: n
    integer :: x
    sphBesselY_di = sphBesselYFunc(dble(n),dble(x))
  end function

  real(dp) function sphBesselY_rr(n,x)
    real :: n, x
    sphBesselY_rr = sphBesselYFunc(dble(n),dble(x))
  end function

  real(dp) function sphBesselY_rd(n,x)
    real :: n
    real(dp) :: x
    sphBesselY_rd = sphBesselYFunc(dble(n),dble(x))
  end function

  real(dp) function sphBesselY_dr(n,x)
    real(dp) :: n
    real :: x
    sphBesselY_dr = sphBesselYFunc(dble(n),dble(x))
  end function

  real(dp) function sphBesselY_dd(n,x)
    real(dp) :: n, x
    sphBesselY_dd = sphBesselYFunc(dble(n),dble(x))
  end function

  !*********************************!
  ! Airy Function of the First Kind !
  !*********************************!

  real(dp) function airyAFunc(x)
    real(dp) :: x, airyResult
    if(x.lt.0.d0) then
      airyResult = besselJ(1.d0/3.d0,(2.d0/3.d0)*(abs(x)**(3.d0/2.d0)))
      airyResult = airyResult+besselJ(-1.d0/3.d0,(2.d0/3.d0)*(abs(x)**(3.d0/2.d0)))
      airyResult = airyResult*sqrt(abs(x)/9.d0)
    elseif(x.eq.0.d0) then
      airyResult = 1.d0/((3.d0**(2.d0/3.d0))*gamma(2.d0/3.d0))
    else
      airyResult = besselK(1.d0/3.d0,(2.d0/3.d0)*(x**(3.d0/2.d0)))
      airyResult = airyResult*sqrt(x/3.d0)/pi_
    end if
    airyAFunc = airyResult
  end function

  real(dp) function airyA_i(x)
    integer :: x
    airyA_i = airyAFunc(dble(x))
  end function

  real(dp) function airyA_r(x)
    real :: x
    airyA_r = airyAFunc(dble(x))
  end function

  real(dp) function airyA_d(x)
    real(dp) :: x
    airyA_d = airyAFunc(dble(x))
  end function

  !**********************************!
  ! Airy Function of the Second Kind !
  !************************&*********!

  real(dp) function airyBFunc(x)
    real(dp) :: x, airyResult
    if(x.lt.0.d0) then
      airyResult = besselJ(-1.d0/3.d0,(2.d0/3.d0)*(abs(x)**(3.d0/2.d0)))
      airyResult = airyResult-besselJ(1.d0/3.d0,(2.d0/3.d0)*(abs(x)**(3.d0/2.d0)))
      airyResult = airyResult*sqrt(abs(x)/3.d0)
    elseif(x.eq.0.d0) then
      airyResult = 1.d0/((3.d0**(1.d0/6.d0))*gamma(2.d0/3.d0))
    else
      airyResult = besselI(1.d0/3.d0,(2.d0/3.d0)*(x**(3.d0/2.d0)))
      airyResult = airyResult+besselI(-1.d0/3.d0,(2.d0/3.d0)*(x**(3.d0/2.d0)))
      airyResult = airyResult*sqrt(x/3.d0)
    end if
    airyBFunc = airyResult
  end function

  real(dp) function airyB_i(x)
    integer :: x
    airyB_i = airyBFunc(dble(x))
  end function

  real(dp) function airyB_r(x)
    real :: x
    airyB_r = airyBFunc(dble(x))
  end function

  real(dp) function airyB_d(x)
    real(dp) :: x
    airyB_d = airyBFunc(dble(x))
  end function

  !****************!
  ! Error Function !
  !****************!

  real(dp) function errf_IntFunc(n,c,x)
    integer :: n
    real(dp) :: c(n), x
    errf_IntFunc = exp(-x*x)
  end function

  real(dp) function errfFunc(x)
    real(dp) :: x, consts(0)
    errfFunc = gaussLegendre(errf_IntFunc,0.d0,x,0,consts)*2.d0/sqrt(pi_)
  end function

  real(dp) function errf_i(x)
    integer :: x
    errf_i = errfFunc(dble(x))
  end function

  real(dp) function errf_r(x)
    real :: x
    errf_r = errfFunc(dble(x))
  end function

  real(dp) function errf_d(x)
    real(dp) :: x
    errf_d = errfFunc(dble(x))
  end function

  !******************************!
  ! Complementary Error Function !
  !******************************!

  real(dp) function errfcFunc(x)
    real(dp) :: x
    errfcFunc = 1.d0-errf(x)
  end function

  real(dp) function errfc_i(x)
    integer :: x
    errfc_i = errfcFunc(dble(x))
  end function

  real(dp) function errfc_r(x)
    real :: x
    errfc_r = errfcFunc(dble(x))
  end function

  real(dp) function errfc_d(x)
    real(dp) :: x
    errfc_d = errfcFunc(dble(x))
  end function

  !***************!
  ! Beta Function !
  !***************!

  real(dp) function beta_IntFunc(n,c,x)
    integer :: n
    real(dp) :: c(n), x, func
    func = x**(c(1)-1.d0)
    func = func*((1-x)**(c(2)-1.d0))
    beta_IntFunc = func
  end function

  real(dp) function betaFunc(x,y)
    real(dp) :: x, y, consts(2)
    consts(1) = x
    consts(2) = y
    betaFunc = gaussLegendre(beta_IntFunc,0.d0,1.d0,2,consts)
  end function

  real(dp) function beta_ii(x,y)
    integer :: x, y
    beta_ii = betaFunc(dble(x),dble(y))
  end function

  real(dp) function beta_ir(x,y)
    integer :: x
    real :: y
    beta_ir = betaFunc(dble(x),dble(y))
  end function

  real(dp) function beta_id(x,y)
    integer :: x
    real(dp) :: y
    beta_id = betaFunc(dble(x),dble(y))
  end function

  real(dp) function beta_ri(x,y)
    real :: x
    integer :: y
    beta_ri = betaFunc(dble(x),dble(y))
  end function

  real(dp) function beta_di(x,y)
    real(dp) :: x
    integer :: y
    beta_di = betaFunc(dble(x),dble(y))
  end function

  real(dp) function beta_rr(x,y)
    real :: x, y
    beta_rr = betaFunc(dble(x),dble(y))
  end function

  real(dp) function beta_rd(x,y)
    real :: x
    real(dp) :: y
    beta_rd = betaFunc(dble(x),dble(y))
  end function

  real(dp) function beta_dr(x,y)
    real(dp) :: x
    real :: y
    beta_dr = betaFunc(dble(x),dble(y))
  end function

  real(dp) function beta_dd(x,y)
    real(dp) :: x, y
    beta_dd = betaFunc(dble(x),dble(y))
  end function

  !**************************!
  ! Incomplete Beta Function !
  !**************************!

  real(dp) function incBetaFunc(x,a,b)
    real(dp) :: x, a, b, consts(2)
    consts(1) = a
    consts(2) = b
    incBetaFunc = gaussLegendre(beta_IntFunc,0.d0,x,2,consts)
  end function

  real(dp) function incBeta_iii(x,a,b)
    integer :: x, a, b
    incBeta_iii = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_iir(x,a,b)
    integer :: x, a
    real :: b
    incBeta_iir = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_iid(x,a,b)
    integer :: x, a
    real(dp) :: b
    incBeta_iid = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_iri(x,a,b)
    integer :: x, b
    real :: a
    incBeta_iri = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_irr(x,a,b)
    integer :: x
    real :: a, b
    incBeta_irr = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_ird(x,a,b)
    integer :: x
    real :: a
    real(dp) :: b
    incBeta_ird = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_idi(x,a,b)
    integer :: x, b
    real(dp) :: a
    incBeta_idi = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_idr(x,a,b)
    integer :: x
    real(dp) :: a
    real :: b
    incBeta_idr = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_idd(x,a,b)
    integer :: x
    real(dp) :: a, b
    incBeta_idd = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_rii(x,a,b)
    real :: x
    integer :: a, b
    incBeta_rii = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_rir(x,a,b)
    real :: x, b
    integer :: a
    incBeta_rir = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_rid(x,a,b)
    real :: x
    integer :: a
    real(dp) :: b
    incBeta_rid = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_rri(x,a,b)
    real :: x, a
    integer :: b
    incBeta_rri = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_rrr(x,a,b)
    real :: x, a, b
    incBeta_rrr = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_rrd(x,a,b)
    real :: x, a
    real(dp) :: b
    incBeta_rrd = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_rdi(x,a,b)
    real :: x
    real(dp) :: a
    integer :: b
    incBeta_rdi = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_rdr(x,a,b)
    real :: x, b
    real(dp) :: a
    incBeta_rdr = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_rdd(x,a,b)
    real :: x
    real(dp) :: a, b
    incBeta_rdd = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_dii(x,a,b)
    real(dp) :: x
    integer :: a, b
    incBeta_dii = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_dir(x,a,b)
    real(dp) :: x
    integer :: a
    real :: b
    incBeta_dir = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_did(x,a,b)
    real(dp) :: x, b
    integer :: a
    incBeta_did = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_dri(x,a,b)
    real(dp) :: x
    real :: a
    integer :: b
    incBeta_dri = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_drr(x,a,b)
    real(dp) :: x
    real :: a, b
    incBeta_drr = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_drd(x,a,b)
    real(dp) :: x, b
    real :: a
    incBeta_drd = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_ddi(x,a,b)
    real(dp) :: x, a
    integer :: b
    incBeta_ddi = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_ddr(x,a,b)
    real(dp) :: x, a
    real :: b
    incBeta_ddr = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  real(dp) function incBeta_ddd(x,a,b)
    real(dp) :: x, a, b
    incBeta_ddd = incBetaFunc(dble(x),dble(a),dble(b))
  end function

  !****************!
  ! Gamma Function !
  !****************!

  real(dp) function gammaFunction(n,c,x)
    integer :: n
    real(dp) :: c(n), x, func
    gammaFunction = (-log(x))**(c(1)-1.d0)
  end function

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
  end function


end module specialFunctions
