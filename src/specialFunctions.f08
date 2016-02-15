!*********************************************************
!>
! Special Functions:
!  * Airy Function - Ai and Bi
!  * Bessel Functions - J and Y
!  * Bessel Function - I and K
!  * Spherical Bessel Functions
!  * Beta Function
!  * Incomplete Beta Function
!  * Chebyshev Polynomial - T and U
!  * Error Function and Complementary Error Function
!  * Factorial
!  * Hermite Polynomial (Physicists' and Probabilists')
!  * Laguerre Polynomial
!  * Associated Laguerre Polynomial
!  * Legendre Polynomial
!  * Sine Integrals
!  * Hyperbolic Sine Integrals
!  * Cosine Integrals
!  * Hyperbolic Cosine Integrals
!  * Hypergeometric Series: 0F1, 1F0, 1F1, 2F1, 3F2
!  * Dilogarithm (3F2)
!  * Add complex inputs to all functions that can have complex input
!
! Not 100% Right:
!  * Gamma Function
!
! Future:
!  * Hahn polynomials (3F2)
!  * Clausen Functions
!  * Coulomb Wave Functions (1F1)
!  * Coupling Coefficients (Wigner 3-j, 6-j, 9-j)
!  * Dawson Function
!  * Debye Functions
!  * Elliptic Integrals
!  * Fermi-Dirac Function
!  * Gegenbauer Functions
!  * Lambert W Functions
!  * Synchrotron Functions
!  * Zeta Functions


module specialFunctions
  use constants
  use integration
  use odeSolver
  implicit none

  private
  integer, parameter :: dp = kind(0.d0)

  public :: pochhammerF, pochhammerR
  public :: besselJ, besselY, besselI, besselK, besselY_dd
  public :: sphBesselJ, sphBesselY
  public :: airyA, airyB
  public :: errf, errfc
  public :: beta, incBeta
  public :: gammaFunc
  public :: legendrePoly
  public :: hermitePoly, hermitePolyProb
  public :: laguerrePoly, assocLaguerrePoly
  public :: chebyPolyT, chebyPolyU
  public :: sini, hypSini
  public :: cosi, hypCosi
  public :: hypGeo0F1, hypGeo0F1Deriv
  public :: hypGeo1F0, hypGeo1F0Deriv
  public :: hypGeo1F1, hypGeo1F1Deriv
  public :: hypGeo2F1, hypGeo2F1Deriv
  public :: hypGeo3F2, hypGeo3F2Deriv
  public :: diLog

  interface pochhammerF
    module procedure pochhammerF_i, pochhammerF_r, pochhammerF_d, pochhammerF_Cr, pochhammerF_Cd
  end interface

  interface pochhammerR
    module procedure pochhammerR_i, pochhammerR_r, pochhammerR_d, pochhammerR_Cr, pochhammerR_Cd
  end interface

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

  interface legendrePoly
    module procedure legendrePoly_i, legendrePoly_r, legendrePoly_d
  end interface

  interface hermitePoly
    module procedure hermitePoly_i, hermitePoly_r, hermitePoly_d
  end interface

  interface hermitePolyProb
    module procedure hermitePolyProb_i, hermitePolyProb_r, hermitePolyProb_d
  end interface

  interface laguerrePoly
    module procedure laguerrePoly_i, laguerrePoly_r, laguerrePoly_d
  end interface

  interface assocLaguerrePoly
    module procedure assocLaguerrePoly_ii, assocLaguerrePoly_ir, assocLaguerrePoly_id, &
      assocLaguerrePoly_ri, assocLaguerrePoly_rr, assocLaguerrePoly_rd, assocLaguerrePoly_di, &
      assocLaguerrePoly_dr, assocLaguerrePoly_dd
  end interface

  interface chebyPolyT
    module procedure chebyPolyT_i, chebyPolyT_r, chebyPolyT_d
  end interface

  interface chebyPolyU
    module procedure chebyPolyU_i, chebyPolyU_r, chebyPolyU_d
  end interface

  interface sini
    module procedure sini_i, sini_r, sini_d
  end interface

  interface hypSini
    module procedure hypSini_i, hypSini_r, hypSini_d
  end interface

  interface cosi
    module procedure cosi_i, cosi_r, cosi_d
  end interface

  interface hypCosi
    module procedure hypCosi_i, hypCosi_r, hypCosi_d
  end interface

  interface hypGeo0F1
    module procedure hypGeo0F1_ii, hypGeo0F1_ir, hypGeo0F1_id, hypGeo0F1_ri, hypGeo0F1_rr, hypGeo0F1_rd, hypGeo0F1_di, &
      hypGeo0F1_dr, hypGeo0F1_dd, hypGeo0F1_iCr, hypGeo0F1_rCr, hypGeo0F1_dCr, hypGeo0F1_iCd, hypGeo0F1_rCd, hypGeo0F1_dCd, &
      hypGeo0F1_CrCr, hypGeo0F1_CrCd, hypGeo0F1_CdCd
  end interface

  interface hypGeo0F1Deriv
    module procedure hypGeo0F1Deriv_ii, hypGeo0F1Deriv_ir, hypGeo0F1Deriv_id, hypGeo0F1Deriv_ri, hypGeo0F1Deriv_rr, &
      hypGeo0F1Deriv_rd, hypGeo0F1Deriv_di, hypGeo0F1Deriv_dr, hypGeo0F1Deriv_dd, hypGeo0F1Deriv_iCr, hypGeo0F1Deriv_rCr, &
      hypGeo0F1Deriv_dCr, hypGeo0F1Deriv_iCd, hypGeo0F1Deriv_rCd, hypGeo0F1Deriv_dCd, hypGeo0F1Deriv_CrCr, hypGeo0F1Deriv_CrCd, &
      hypGeo0F1Deriv_CdCd
  end interface

  interface hypGeo1F0
    module procedure hypGeo1F0_ii, hypGeo1F0_ir, hypGeo1F0_id, hypGeo1F0_ri, hypGeo1F0_rr, hypGeo1F0_rd, hypGeo1F0_di, &
      hypGeo1F0_dr, hypGeo1F0_dd, hypGeo1F0_iCr, hypGeo1F0_rCr, hypGeo1F0_dCr, hypGeo1F0_iCd, hypGeo1F0_rCd, hypGeo1F0_dCd, &
      hypGeo1F0_CrCr, hypGeo1F0_CrCd, hypGeo1F0_CdCd
  end interface

  interface hypGeo1F0Deriv
    module procedure hypGeo1F0Deriv_ii, hypGeo1F0Deriv_ir, hypGeo1F0Deriv_id, hypGeo1F0Deriv_ri, hypGeo1F0Deriv_rr, &
      hypGeo1F0Deriv_rd, hypGeo1F0Deriv_di, hypGeo1F0Deriv_dr, hypGeo1F0Deriv_dd, hypGeo1F0Deriv_iCr, hypGeo1F0Deriv_rCr, &
      hypGeo1F0Deriv_dCr, hypGeo1F0Deriv_iCd, hypGeo1F0Deriv_rCd, hypGeo1F0Deriv_dCd, hypGeo1F0Deriv_CrCr, hypGeo1F0Deriv_CrCd, &
      hypGeo1F0Deriv_CdCd
  end interface

  interface hypGeo1F1
    module procedure hypGeo1F1_iii, hypGeo1F1_rri, hypGeo1F1_ddi, hypGeo1F1_iir, hypGeo1F1_rrr, hypGeo1F1_ddr, &
      hypGeo1F1_iid, hypGeo1F1_rrd, hypGeo1F1_ddd, hypGeo1F1_iiCr, hypGeo1F1_rrCr, hypGeo1F1_ddCr, hypGeo1F1_iiCd, &
      hypGeo1F1_rrCd, hypGeo1F1_ddCd, hypGeo1F1_CrCrCr, hypGeo1F1_CrCrCd, hypGeo1F1_CdCdCd
  end interface

  interface hypGeo1F1Deriv
    module procedure hypGeo1F1Deriv_iii, hypGeo1F1Deriv_rri, hypGeo1F1Deriv_ddi, hypGeo1F1Deriv_iir, hypGeo1F1Deriv_rrr, &
      hypGeo1F1Deriv_ddr, hypGeo1F1Deriv_iid, hypGeo1F1Deriv_rrd, hypGeo1F1Deriv_ddd, hypGeo1F1Deriv_iiCr, hypGeo1F1Deriv_rrCr, &
      hypGeo1F1Deriv_ddCr, hypGeo1F1Deriv_iiCd, hypGeo1F1Deriv_rrCd, hypGeo1F1Deriv_ddCd, hypGeo1F1Deriv_CrCrCr, &
      hypGeo1F1Deriv_CrCrCd, hypGeo1F1Deriv_CdCdCd
  end interface

  interface hypGeo2F1
    module procedure hypGeo2F1_iiii, hypGeo2F1_iiir, hypGeo2F1_iiid, hypGeo2F1_rrri, &
        hypGeo2F1_rrrr, hypGeo2F1_rrrd, hypGeo2F1_dddi, hypGeo2F1_dddr, hypGeo2F1_dddd, &
        hypGeo2F1_iiiCr, hypGeo2F1_rrrCr, hypGeo2F1_dddCr, hypGeo2F1_iiiCd, hypGeo2F1_rrrCd, &
        hypGeo2F1_dddCd, hypGeo2F1_CrCrCrCr, hypGeo2F1_CrCrCrCd, hypGeo2F1_CdCdCdCd
  end interface

  interface hypGeo2F1Deriv
    module procedure hypGeo2F1Deriv_iiii, hypGeo2F1Deriv_iiir, hypGeo2F1Deriv_iiid, hypGeo2F1Deriv_rrri, &
      hypGeo2F1Deriv_rrrr, hypGeo2F1Deriv_rrrd, hypGeo2F1Deriv_dddi, hypGeo2F1Deriv_dddr, hypGeo2F1Deriv_dddd, &
      hypGeo2F1Deriv_iiiCr, hypGeo2F1Deriv_rrrCr, hypGeo2F1Deriv_dddCr, hypGeo2F1Deriv_iiiCd, hypGeo2F1Deriv_rrrCd, &
      hypGeo2F1Deriv_dddCd, hypGeo2F1Deriv_CrCrCrCr, hypGeo2F1Deriv_CrCrCrCd, hypGeo2F1Deriv_CdCdCdCd
  end interface

  interface hypGeo3F2
    module procedure hypGeo3F2_iiiiii, hypGeo3F2_iiiiir, hypGeo3F2_iiiiid, hypGeo3F2_rrrrri, hypGeo3F2_rrrrrr, hypGeo3F2_rrrrrd, &
      hypGeo3F2_dddddi, hypGeo3F2_dddddr, hypGeo3F2_dddddd, hypGeo3F2_iiiiiCr, hypGeo3F2_iiiiiCd, hypGeo3F2_rrrrrCr, &
      hypGeo3F2_rrrrrCd, hypGeo3F2_dddddCr, hypGeo3F2_dddddCd, hypGeo3F2_CrCrCrCrCrCr, hypGeo3F2_CrCrCrCrCrCd, &
      hypGeo3F2_CdCdCdCdCdCd
  end interface

  interface hypGeo3F2Deriv
    module procedure hypGeo3F2Deriv_iiiiii, hypGeo3F2Deriv_iiiiir, hypGeo3F2Deriv_iiiiid, hypGeo3F2Deriv_rrrrri, &
      hypGeo3F2Deriv_rrrrrr, hypGeo3F2Deriv_rrrrrd, hypGeo3F2Deriv_dddddi, hypGeo3F2Deriv_dddddr, hypGeo3F2Deriv_dddddd, &
      hypGeo3F2Deriv_iiiiiCr, hypGeo3F2Deriv_iiiiiCd, hypGeo3F2Deriv_rrrrrCr, hypGeo3F2Deriv_rrrrrCd, hypGeo3F2Deriv_dddddCr, &
      hypGeo3F2Deriv_dddddCd, hypGeo3F2Deriv_CrCrCrCrCrCr, hypGeo3F2Deriv_CrCrCrCrCrCd, hypGeo3F2Deriv_CdCdCdCdCdCd
  end interface

  interface diLog
    module procedure diLog_i, diLog_r, diLog_d, diLog_Cr, diLog_Cd
  end interface

contains

  !********************!
  ! Pochhammer Falling !
  !********************!

  complex(dp) function pochhammerFFunc(x,n)
    integer :: i, n
    complex(dp) :: x, result
    if(n.eq.0) then
      pochhammerFFunc = 1.d0
    else
      result = x
      do i=1,n+1
        result = result*(x-dble(i))
      end do
      pochhammerFFunc = result
    end if
  end function

  real(dp) function pochhammerF_i(x,n)
    integer :: x, n
    pochhammerF_i = real(pochhammerFFunc(cmplx(real(x),0.d0,dp),n))
  end function

  real(dp) function pochhammerF_r(x,n)
    integer :: n
    real :: x
    pochhammerF_r = real(pochhammerFFunc(cmplx(real(x),0.d0,dp),n))
  end function

  real(dp) function pochhammerF_d(x,n)
    integer :: n
    real(dp) :: x
    pochhammerF_d = real(pochhammerFFunc(cmplx(real(x),0.d0,dp),n))
  end function

  complex(dp) function pochhammerF_Cr(x,n)
    integer :: n
    complex :: x
    pochhammerF_Cr = real(pochhammerFFunc(cmplx(real(x),aimag(x),dp),n))
  end function

  complex(dp) function pochhammerF_Cd(x,n)
    integer :: n
    complex(dp) :: x
    pochhammerF_Cd = real(pochhammerFFunc(x,n))
  end function

  !*******************!
  ! Pochhammer Rising !
  !*******************!

  complex(dp) function pochhammerRFunc(x,n)
    integer :: i, n
    complex(dp) :: x, result
    if(n.eq.0) then
      pochhammerRFunc = cmplx(1.d0,0.d0,dp)
    else
      result = x
      do i=1,n-1
        result = result*(x+i)
      end do
      pochhammerRFunc = result
    end if
  end function

  real(dp) function pochhammerR_i(x,n)
    integer :: x, n
    pochhammerR_i = real(pochhammerRFunc(cmplx(x,0.d0,dp),n))
  end function

  real(dp) function pochhammerR_r(x,n)
    integer :: n
    real :: x
    pochhammerR_r = real(pochhammerRFunc(cmplx(x,0.d0,dp),n))
  end function

  real(dp) function pochhammerR_d(x,n)
    integer :: n
    real(dp) :: x
    pochhammerR_d = real(pochhammerRFunc(cmplx(x,0.d0,dp),n))
  end function

  complex(dp) function pochhammerR_Cr(x,n)
    integer :: n
    complex :: x
    pochhammerR_Cr = pochhammerRFunc(cmplx(real(x),aimag(x),dp),n)
  end function

  complex(dp) function pochhammerR_Cd(x,n)
    integer :: n
    complex(dp) :: x
    pochhammerR_Cd = pochhammerRFunc(x,n)
  end function

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
    real(dp) :: c(n),x
    besselK_IntFunc = (exp(-c(2)*cosh(-log(x)))*cosh(-c(1)*log(x)))/x
  end function besselK_IntFunc

  real(dp) function besselKFunc(n,x)
    real(dp) :: n, x, consts(2)
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
    func = func*((1.d0-x)**(c(2)-1.d0))
    beta_IntFunc = func
  end function

  real(dp) function betaFunc(x,y)
    real(dp) :: x, y, consts(2), intResults
    consts(1) = x
    consts(2) = y
    !betaFunc = gaussLegendre(beta_IntFunc,0.d0,1.d0,2,consts)
    betaFunc = gamma(x)*gamma(y)/gamma(x+y)
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
    real(dp) :: c(n), x
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

  !**********************!
  ! Legendre Polynomials !
  !**********************!

  recursive function legendrePolyFunc(n,x) result(res)
    integer :: n
    real(dp) :: x, res
    if(n.eq.0) then
      res = 1
    elseif(n.eq.1) then
      res = x
    else
      res = (2*(n-1)+1)*x*legendrePolyFunc(n-1,x)
      res = res - (n-1)*legendrePolyFunc(n-2,x)
      res = res/n
    end if
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

  !***********************************!
  ! Hermite Polynomials (Physicists') !
  !***********************************!

  recursive function hermitePolyFunc(n,x) result(res)
    integer :: n
    real(dp) :: x, res
    if(n.eq.0) then
      res = 1.d0
    elseif(n.eq.1) then
      res = 2.d0*x
    else
      res = 2.d0*x*hermitePolyFunc(n-1,x)
      res = res - 2.d0*(n-1)*hermitePolyFunc(n-2,x)
    end if
  end function

  real(dp) function hermitePoly_i(n,x)
    integer :: n, x
    hermitePoly_i = hermitePolyFunc(n,dble(x))
  end function

  real(dp) function hermitePoly_r(n,x)
    integer :: n
    real :: x
    hermitePoly_r = hermitePolyFunc(n,dble(x))
  end function

  real(dp) function hermitePoly_d(n,x)
    integer :: n
    real(dp) :: x
    hermitePoly_d = hermitePolyFunc(n,dble(x))
  end function

  !*************************************!
  ! Hermite Polynomials (Probabilists') !
  !*************************************!

  recursive function hermitePolyProbFunc(n,x) result(res)
    integer :: n
    real(dp) :: x, res
    if(n.eq.0) then
      res = 1.d0
    elseif(n.eq.1) then
      res = x
    else
      res = x*hermitePolyProbFunc(n-1,x)
      res = res - (n-1)*hermitePolyProbFunc(n-2,x)
    end if
  end function

  real(dp) function hermitePolyProb_i(n,x)
    integer :: n, x
    hermitePolyProb_i = hermitePolyProbFunc(n,dble(x))
  end function

  real(dp) function hermitePolyProb_r(n,x)
    integer :: n
    real :: x
    hermitePolyProb_r = hermitePolyProbFunc(n,dble(x))
  end function

  real(dp) function hermitePolyProb_d(n,x)
    integer :: n
    real(dp) :: x
    hermitePolyProb_d = hermitePolyProbFunc(n,dble(x))
  end function

  !**********************!
  ! Laguerre Polynomials !
  !**********************!

  recursive function laguerrePolyFunc(n,x) result(res)
    integer :: n
    real(dp) :: x, res
    if(n.eq.0) then
      res = 1.d0
    elseif(n.eq.1) then
      res = 1.d0 - x
    else
      res = (2.d0*(n-1)+1.d0-x)*laguerrePolyFunc(n-1,x)
      res = res - (n-1)*laguerrePolyFunc(n-2,x)
      res = res/n
    end if
  end function

  real(dp) function laguerrePoly_i(n,x)
    integer :: n, x
    laguerrePoly_i = laguerrePolyFunc(n,dble(x))
  end function

  real(dp) function laguerrePoly_r(n,x)
    integer :: n
    real :: x
    laguerrePoly_r = laguerrePolyFunc(n,dble(x))
  end function

  real(dp) function laguerrePoly_d(n,x)
    integer :: n
    real(dp) :: x
    laguerrePoly_d = laguerrePolyFunc(n,dble(x))
  end function

  !*********************************!
  ! Associated Laguerre Polynomials !
  !*********************************!

  recursive function assocLaguerrePolyFunc(n,a,x) result(res)
    integer :: n
    real(dp) :: a, x, res
    if(n.eq.0) then
      res = 1.d0
    elseif(n.eq.1) then
      res = 1.d0 + a - x
    else
      res = (2.d0*(n-1)+1.d0+a-x)*assocLaguerrePolyFunc(n-1,a,x)
      res = res - (n-1+a)*assocLaguerrePolyFunc(n-2,a,x)
      res = res/n
    end if
  end function

  real(dp) function assocLaguerrePoly_ii(n,a,x)
    integer :: n, a, x
    assocLaguerrePoly_ii = assocLaguerrePolyFunc(n,dble(a),dble(x))
  end function

  real(dp) function assocLaguerrePoly_ir(n,a,x)
    integer :: n, a
    real :: x
    assocLaguerrePoly_ir = assocLaguerrePolyFunc(n,dble(a),dble(x))
  end function

  real(dp) function assocLaguerrePoly_id(n,a,x)
    integer :: n, a
    real(dp) :: x
    assocLaguerrePoly_id = assocLaguerrePolyFunc(n,dble(a),dble(x))
  end function

  real(dp) function assocLaguerrePoly_ri(n,a,x)
    integer :: n, x
    real :: a
    assocLaguerrePoly_ri = assocLaguerrePolyFunc(n,dble(a),dble(x))
  end function

  real(dp) function assocLaguerrePoly_rr(n,a,x)
    integer :: n
    real :: a, x
    assocLaguerrePoly_rr = assocLaguerrePolyFunc(n,dble(a),dble(x))
  end function

  real(dp) function assocLaguerrePoly_rd(n,a,x)
    integer :: n
    real :: a
    real(dp) :: x
    assocLaguerrePoly_rd = assocLaguerrePolyFunc(n,dble(a),dble(x))
  end function

  real(dp) function assocLaguerrePoly_di(n,a,x)
    integer :: n, x
    real(dp) :: a
    assocLaguerrePoly_di = assocLaguerrePolyFunc(n,dble(a),dble(x))
  end function

  real(dp) function assocLaguerrePoly_dr(n,a,x)
    integer :: n
    real(dp) :: a
    real :: x
    assocLaguerrePoly_dr = assocLaguerrePolyFunc(n,dble(a),dble(x))
  end function

  real(dp) function assocLaguerrePoly_dd(n,a,x)
    integer :: n
    real(dp) :: a, x
    assocLaguerrePoly_dd = assocLaguerrePolyFunc(n,dble(a),dble(x))
  end function

  !*****************************************!
  ! Chebyshev Polynomials of the First Kind !
  !*****************************************!

  recursive function chebyPolyTFunc(n,x) result(res)
    integer :: n
    real(dp) :: x, res
    if(n.eq.0) then
      res = 1.d0
    elseif(n.eq.1) then
      res = x
    else
      res = 2.d0*x*chebyPolyTFunc(n-1,x)
      res = res - chebyPolyTFunc(n-2,x)
    end if
  end function

  real(dp) function chebyPolyT_i(n,x)
    integer :: n, x
    chebyPolyT_i = chebyPolyTFunc(n,dble(x))
  end function

  real(dp) function chebyPolyT_r(n,x)
    integer :: n
    real :: x
    chebyPolyT_r = chebyPolyTFunc(n,dble(x))
  end function

  real(dp) function chebyPolyT_d(n,x)
    integer :: n
    real(dp) :: x
    chebyPolyT_d = chebyPolyTFunc(n,dble(x))
  end function

  !******************************************!
  ! Chebyshev Polynomials of the Second Kind !
  !******************************************!

  recursive function chebyPolyUFunc(n,x) result(res)
    integer :: n
    real(dp) :: x, res
    if(n.eq.0) then
      res = 1.d0
    elseif(n.eq.1) then
      res = 2.d0*x
    else
      res = 2.d0*x*chebyPolyUFunc(n-1,x)
      res = res - chebyPolyUFunc(n-2,x)
    end if
  end function

  real(dp) function chebyPolyU_i(n,x)
    integer :: n, x
    chebyPolyU_i = chebyPolyUFunc(n,dble(x))
  end function

  real(dp) function chebyPolyU_r(n,x)
    integer :: n
    real :: x
    chebyPolyU_r = chebyPolyUFunc(n,dble(x))
  end function

  real(dp) function chebyPolyU_d(n,x)
    integer :: n
    real(dp) :: x
    chebyPolyU_d = chebyPolyUFunc(n,dble(x))
  end function

  !****************!
  ! Sine Integrals !
  !****************!

  real(dp) function sini_IntFunc(n,c,x)
    integer :: n
    real(dp) :: c(n),x
    sini_IntFunc = sin(x)/x
  end function

  real(dp) function siniFunc(x)
    real(dp) :: x, consts(0)
    siniFunc = gaussLegendre(sini_IntFunc,0.d0,x,0,consts)
  end function

  real(dp) function sini_i(x)
    integer :: x
    sini_i = siniFunc(dble(x))
  end function

  real(dp) function sini_r(x)
    real :: x
    sini_r = siniFunc(dble(x))
  end function

  real(dp) function sini_d(x)
    real(dp) :: x
    sini_d = siniFunc(dble(x))
  end function

  !***************************!
  ! Hyperbolic Sine Integrals !
  !***************************!

  real(dp) function hypSini_IntFunc(n,c,x)
    integer :: n
    real(dp) :: c(n),x
    hypSini_IntFunc = sinh(x)/x
  end function

  real(dp) function hypSiniFunc(x)
    real(dp) :: x, consts(0)
    hypSiniFunc = gaussLegendre(hypSini_IntFunc,0.d0,x,0,consts)
  end function

  real(dp) function hypSini_i(x)
    integer :: x
    hypSini_i = hypSiniFunc(dble(x))
  end function

  real(dp) function hypSini_r(x)
    real :: x
    hypSini_r = hypSiniFunc(dble(x))
  end function

  real(dp) function hypSini_d(x)
    real(dp) :: x
    hypSini_d = hypSiniFunc(dble(x))
  end function

  !******************!
  ! Cosine Integrals !
  !******************!

  real(dp) function cosi_IntFunc(n,c,x)
    integer :: n
    real(dp) :: c(n),x
    cosi_IntFunc = (cos(x) - 1.d0)/x
  end function

  real(dp) function cosiFunc(x)
    real(dp) :: x, consts(0)
    cosiFunc = em_ + log(x) + gaussLegendre(cosi_IntFunc,0.d0,x,0,consts)
  end function

  real(dp) function cosi_i(x)
    integer :: x
    cosi_i = cosiFunc(dble(x))
  end function

  real(dp) function cosi_r(x)
    real :: x
    cosi_r = cosiFunc(dble(x))
  end function

  real(dp) function cosi_d(x)
    real(dp) :: x
    cosi_d = cosiFunc(dble(x))
  end function

  !*****************************!
  ! Hyperbolic Cosine Integrals !
  !*****************************!

  real(dp) function hypCosi_IntFunc(n,c,x)
    integer :: n
    real(dp) :: c(n),x
    hypCosi_IntFunc = (cosh(x) - 1.d0)/x
  end function

  real(dp) function hypCosiFunc(x)
    real(dp) :: x, consts(0)
    hypCosiFunc = em_ + log(x) + gaussLegendre(hypCosi_IntFunc,0.d0,x,0,consts)
  end function

  real(dp) function hypCosi_i(x)
    integer :: x
    hypCosi_i = hypCosiFunc(dble(x))
  end function

  real(dp) function hypCosi_r(x)
    real :: x
    hypCosi_r = hypCosiFunc(dble(x))
  end function

  real(dp) function hypCosi_d(x)
    real(dp) :: x
    hypCosi_d = hypCosiFunc(dble(x))
  end function

  !*****************************!
  ! Hypergeometric Function 0F1 !
  !*****************************!

  function hypGeo0F1_odeFunc(nF,nC,c,x,y)
    integer :: nF, nC
    real(dp) :: hypGeo0F1_odeFunc(nF), x, y(nF)
    complex(dp) :: c(nC), zs, fT, fpT, f, fp, dz
    dz = c(2)-c(1)
    zs = c(1) + x*dz
    fT = cmplx(y(1),y(2),dp); fpT = cmplx(y(3),y(4),dp)
    f = dz*fpT
    fp = (fT-c(3)*fpT)*(dz/zs)
    hypGeo0F1_odeFunc(1) = real(f); hypGeo0F1_odeFunc(2) = aimag(f)
    hypGeo0F1_odeFunc(3) = real(fp); hypGeo0F1_odeFunc(4) = aimag(fp)
  end function

  complex(dp) function hypGeo0F1_series(a,z)
    !! Returns Hypergeometric Function 0F1
    !! Should be only used for |z|<1.0
    integer :: i
    complex(dp) :: a, z, result, oldResult, indvResult
    result = cmplx(0.d0,0.d0,dp)
    oldResult = result
    do i=0, 100000
      indvResult = (z**i)/(pochhammerR(a,i)*gamma(dble(i)+1.d0))
      if (indvResult /= indvResult) exit
      result = result + indvResult
      if (abs(result-oldResult).le.1.d-16) exit
      oldResult = result
    end do
  hypGeo0F1_series = result
  end function

  complex(dp) function hypGeo0F1Deriv_series(a,z)
    !! Returns Derivative of the Hypergeometric Function 0F1
    !! Should be only used for |z|<1.0
    complex(dp) :: a, z
    hypGeo0F1Deriv_series = hypGeo1F0_series(a+1.d0,z)/a
  end function

  complex(dp) function hypGeo0F1Func(a,z)
    complex(dp) :: a, z
    real(dp) :: y0(4), odeResult(4)
    complex(dp) :: z0, series(2), consts(3)
    if(abs(z).lt.1.d0) then
      hypGeo0F1Func = hypGeo0F1_series(a,z)
      return
    else if(real(z).gt.0.d0 .and. real(z).lt.1.d0) then
      z0 = cmplx(0.5d0,0.d0,dp)
    else if(real(z).gt.-1.d0 .and. real(z).lt.0.d0) then
      z0 = cmplx(-0.5d0,0.d0,dp)
    else if(real(z).gt.1.d0) then
      z0 = cmplx(0.d0,0.5d0,dp)
    else if(real(z).lt.-1.d0) then
      z0 = cmplx(0.d0,-0.5d0,dp)
    end if
    series(1) = hypGeo0F1_series(a,z0)
    series(2) = hypGeo0F1Deriv_series(a,z0)
    consts(1) = z0; consts(2) = z; consts(3) = a;
    y0(1) = real(series(1)); y0(2) = aimag(series(1))
    y0(3) = real(series(2)); y0(4) = aimag(series(2))
    odeResult = rk1AdaptStepCmplxC(hypGeo0F1_odeFunc,4,0.d0,1.d0,y0,3,consts)
    hypGeo0F1Func = cmplx(odeResult(1),odeResult(2),dp)
  end function

  complex(dp) function hypGeo0F1DerivFunc(a,z)
    complex(dp) :: a, z
    hypGeo0F1DerivFunc = hypGeo0F1Func(a+1.d0,z)/a
  end function

  complex(dp) function hypGeo0F1_ii(a,z)
    integer :: a, z
    hypGeo0F1_ii = hypGeo0F1Func(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo0F1_ri(a,z)
    integer :: z
    real :: a
    hypGeo0F1_ri = hypGeo0F1Func(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo0F1_di(a,z)
    integer :: z
    real(dp) :: a
    hypGeo0F1_di = hypGeo0F1Func(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo0F1_ir(a,z)
    integer :: a
    real :: z
    hypGeo0F1_ir = hypGeo0F1Func(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo0F1_rr(a,z)
    real :: a, z
    hypGeo0F1_rr = hypGeo0F1Func(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo0F1_dr(a,z)
    real :: z
    real(dp) :: a
    hypGeo0F1_dr = hypGeo0F1Func(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo0F1_id(a,z)
    integer :: a
    real(dp) :: z
    hypGeo0F1_id = hypGeo0F1Func(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo0F1_rd(a,z)
    real :: a
    real(dp) :: z
    hypGeo0F1_rd = hypGeo0F1Func(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo0F1_dd(a,z)
    real(dp) :: a, z
    hypGeo0F1_dd = hypGeo0F1Func(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo0F1_iCr(a,z)
    integer :: a
    complex :: z
    hypGeo0F1_iCr = hypGeo0F1Func(cmplx(a,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo0F1_rCr(a,z)
    real :: a
    complex :: z
    hypGeo0F1_rCr = hypGeo0F1Func(cmplx(a,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo0F1_dCr(a,z)
    real(dp) :: a
    complex :: z
    hypGeo0F1_dCr = hypGeo0F1Func(cmplx(a,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo0F1_iCd(a,z)
    integer :: a
    complex(dp) :: z
    hypGeo0F1_iCd = hypGeo0F1Func(cmplx(a,0.d0,dp),z)
  end function

  complex(dp) function hypGeo0F1_rCd(a,z)
    real :: a
    complex(dp) :: z
    hypGeo0F1_rCd = hypGeo0F1Func(cmplx(a,0.d0,dp),z)
  end function

  complex(dp) function hypGeo0F1_dCd(a,z)
    real(dp) :: a
    complex(dp) :: z
    hypGeo0F1_dCd = hypGeo0F1Func(cmplx(a,0.d0,dp),z)
  end function

  complex(dp) function hypGeo0F1_CrCr(a,z)
    complex :: a, z
    hypGeo0F1_CrCr = hypGeo0F1Func(cmplx(real(a),aimag(a),dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo0F1_CrCd(a,z)
    complex :: a
    complex(dp) :: z
    hypGeo0F1_CrCd = hypGeo0F1Func(cmplx(real(a),aimag(a),dp),z)
  end function

  complex(dp) function hypGeo0F1_CdCd(a,z)
    complex(dp) :: a, z
    hypGeo0F1_CdCd = hypGeo0F1Func(a,z)
  end function

  complex(dp) function hypGeo0F1Deriv_ii(a,z)
    integer :: a, z
    hypGeo0F1Deriv_ii = hypGeo0F1DerivFunc(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo0F1Deriv_ri(a,z)
    integer :: z
    real :: a
    hypGeo0F1Deriv_ri = hypGeo0F1DerivFunc(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo0F1Deriv_di(a,z)
    integer :: z
    real(dp) :: a
    hypGeo0F1Deriv_di = hypGeo0F1DerivFunc(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo0F1Deriv_ir(a,z)
    integer :: a
    real :: z
    hypGeo0F1Deriv_ir = hypGeo0F1DerivFunc(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo0F1Deriv_rr(a,z)
    real :: a, z
    hypGeo0F1Deriv_rr = hypGeo0F1DerivFunc(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo0F1Deriv_dr(a,z)
    real :: z
    real(dp) :: a
    hypGeo0F1Deriv_dr = hypGeo0F1DerivFunc(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo0F1Deriv_id(a,z)
    integer :: a
    real(dp) :: z
    hypGeo0F1Deriv_id = hypGeo0F1DerivFunc(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo0F1Deriv_rd(a,z)
    real :: a
    real(dp) :: z
    hypGeo0F1Deriv_rd = hypGeo0F1DerivFunc(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo0F1Deriv_dd(a,z)
    real(dp) :: a, z
    hypGeo0F1Deriv_dd = hypGeo0F1DerivFunc(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo0F1Deriv_iCr(a,z)
    integer :: a
    complex :: z
    hypGeo0F1Deriv_iCr = hypGeo0F1DerivFunc(cmplx(a,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo0F1Deriv_rCr(a,z)
    real :: a
    complex :: z
    hypGeo0F1Deriv_rCr = hypGeo0F1DerivFunc(cmplx(a,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo0F1Deriv_dCr(a,z)
    real(dp) :: a
    complex :: z
    hypGeo0F1Deriv_dCr = hypGeo0F1DerivFunc(cmplx(a,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo0F1Deriv_iCd(a,z)
    integer :: a
    complex(dp) :: z
    hypGeo0F1Deriv_iCd = hypGeo0F1DerivFunc(cmplx(a,0.d0,dp),z)
  end function

  complex(dp) function hypGeo0F1Deriv_rCd(a,z)
    real :: a
    complex(dp) :: z
    hypGeo0F1Deriv_rCd = hypGeo0F1DerivFunc(cmplx(a,0.d0,dp),z)
  end function

  complex(dp) function hypGeo0F1Deriv_dCd(a,z)
    real(dp) :: a
    complex(dp) :: z
    hypGeo0F1Deriv_dCd = hypGeo0F1DerivFunc(cmplx(a,0.d0,dp),z)
  end function

  complex(dp) function hypGeo0F1Deriv_CrCr(a,z)
    complex :: a, z
    hypGeo0F1Deriv_CrCr = hypGeo0F1DerivFunc(cmplx(real(a),aimag(a),dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo0F1Deriv_CrCd(a,z)
    complex :: a
    complex(dp) :: z
    hypGeo0F1Deriv_CrCd = hypGeo0F1DerivFunc(cmplx(real(a),aimag(a),dp),z)
  end function

  complex(dp) function hypGeo0F1Deriv_CdCd(a,z)
    complex(dp) :: a, z
    hypGeo0F1Deriv_CdCd = hypGeo0F1DerivFunc(a,z)
  end function

  !*****************************!
  ! Hypergeometric Function 1F0 !
  !*****************************!

  function hypGeo1F0_odeFunc(nF,nC,c,x,y)
    integer :: nF, nC
    real(dp) :: hypGeo1F0_odeFunc(nF), x, y(nF)
    complex(dp) :: c(nC), zs, fT, f, dz
    dz = c(2)-c(1)
    zs = c(1) + x*dz
    fT = cmplx(y(1),y(2),dp)
    f = dz*fT*c(3)/(1.d0-zs)
    hypGeo1F0_odeFunc(1) = real(f)
    hypGeo1F0_odeFunc(2) = aimag(f)
  end function

  complex(dp) function hypGeo1F0_series(a,z)
    !! Returns Hypergeometric Function 1F0
    !! Should be only used for |z|<1.0
    integer :: i
    complex(dp) :: a, z, result, oldResult, indvResult
    result = cmplx(0.d0,0.d0,dp)
    oldResult = result
    do i=0, 100000
      indvResult = pochhammerR(a,i)*(z**i)/gamma(dble(i)+1.d0)
      if (indvResult /= indvResult) exit
      result = result + indvResult
      if (abs(result-oldResult).le.1.d-16) exit
      oldResult = result
    end do
  hypGeo1F0_series = result
  end function

  complex(dp) function hypGeo1F0Deriv_series(a,z)
    !! Returns Derivative of the Hypergeometric Function 1F0
    !! Should be only used for |z|<1.0
    complex(dp) :: a, z
    hypGeo1F0Deriv_series = a*hypGeo1F0_series(a+1.d0,z)
  end function

  complex(dp) function hypGeo1F0Func(a,z)
    complex(dp) :: a, z
    real(dp) :: y0(2), odeResult(2)
    complex(dp) :: z0, series, consts(3)
    z0 = cmplx(0.d0,0.d0,dp)
    if(abs(z).lt.1.d0) then
      hypGeo1F0Func = hypGeo1F0_series(a,z)
      return
    else if(real(z).gt.0.d0 .and. real(z).lt.1.d0) then
      z0 = cmplx(0.5d0,0.d0,dp)
    else if(real(z).gt.-1.d0 .and. real(z).lt.0.d0) then
      z0 = cmplx(-0.5d0,0.d0,dp)
    else if(real(z).gt.1.d0) then
      z0 = cmplx(0.d0,0.5d0,dp)
    else if(real(z).lt.-1.d0) then
      z0 = cmplx(0.d0,-0.5d0,dp)
    end if
    series = hypGeo1F0_series(a,z0)
    consts(1) = z0; consts(2) = z; consts(3) = a;
    y0(1) = real(series); y0(2) = aimag(series)
    odeResult = rk1AdaptStepCmplxC(hypGeo1F0_odeFunc,2,0.d0,1.d0,y0,3,consts)
    hypGeo1F0Func = cmplx(odeResult(1),odeResult(2),dp)
  end function

  complex(dp) function hypGeo1F0DerivFunc(a,z)
    complex(dp) :: a, z
    hypGeo1F0DerivFunc = a*hypGeo1F0Func(a+1.d0,z)
  end function

  complex(dp) function hypGeo1F0_ii(a,z)
    integer :: a, z
    hypGeo1F0_ii = hypGeo1F0Func(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F0_ri(a,z)
    integer :: z
    real :: a
    hypGeo1F0_ri = hypGeo1F0Func(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F0_di(a,z)
    integer :: z
    real(dp) :: a
    hypGeo1F0_di = hypGeo1F0Func(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F0_ir(a,z)
    integer :: a
    real :: z
    hypGeo1F0_ir = hypGeo1F0Func(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F0_rr(a,z)
    real :: a, z
    hypGeo1F0_rr = hypGeo1F0Func(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F0_dr(a,z)
    real :: z
    real(dp) :: a
    hypGeo1F0_dr = hypGeo1F0Func(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F0_id(a,z)
    integer :: a
    real(dp) :: z
    hypGeo1F0_id = hypGeo1F0Func(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F0_rd(a,z)
    real :: a
    real(dp) :: z
    hypGeo1F0_rd = hypGeo1F0Func(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F0_dd(a,z)
    real(dp) :: a, z
    hypGeo1F0_dd = hypGeo1F0Func(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F0_iCr(a,z)
    integer :: a
    complex :: z
    hypGeo1F0_iCr = hypGeo1F0Func(cmplx(a,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo1F0_rCr(a,z)
    real :: a
    complex :: z
    hypGeo1F0_rCr = hypGeo1F0Func(cmplx(a,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo1F0_dCr(a,z)
    real(dp) :: a
    complex :: z
    hypGeo1F0_dCr = hypGeo1F0Func(cmplx(a,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo1F0_iCd(a,z)
    integer :: a
    complex(dp) :: z
    hypGeo1F0_iCd = hypGeo1F0Func(cmplx(a,0.d0,dp),z)
  end function

  complex(dp) function hypGeo1F0_rCd(a,z)
    real :: a
    complex(dp) :: z
    hypGeo1F0_rCd = hypGeo1F0Func(cmplx(a,0.d0,dp),z)
  end function

  complex(dp) function hypGeo1F0_dCd(a,z)
    real(dp) :: a
    complex(dp) :: z
    hypGeo1F0_dCd = hypGeo1F0Func(cmplx(a,0.d0,dp),z)
  end function

  complex(dp) function hypGeo1F0_CrCr(a,z)
    complex :: a, z
    hypGeo1F0_CrCr = hypGeo1F0Func(cmplx(real(a),aimag(a),dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo1F0_CrCd(a,z)
    complex :: a
    complex(dp) :: z
    hypGeo1F0_CrCd = hypGeo1F0Func(cmplx(real(a),aimag(a),dp),z)
  end function

  complex(dp) function hypGeo1F0_CdCd(a,z)
    complex(dp) :: a, z
    hypGeo1F0_CdCd = hypGeo1F0Func(a,z)
  end function

  complex(dp) function hypGeo1F0Deriv_ii(a,z)
    integer :: a, z
    hypGeo1F0Deriv_ii = hypGeo1F0DerivFunc(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F0Deriv_ri(a,z)
    integer :: z
    real :: a
    hypGeo1F0Deriv_ri = hypGeo1F0DerivFunc(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F0Deriv_di(a,z)
    integer :: z
    real(dp) :: a
    hypGeo1F0Deriv_di = hypGeo1F0DerivFunc(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F0Deriv_ir(a,z)
    integer :: a
    real :: z
    hypGeo1F0Deriv_ir = hypGeo1F0DerivFunc(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F0Deriv_rr(a,z)
    real :: a, z
    hypGeo1F0Deriv_rr = hypGeo1F0DerivFunc(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F0Deriv_dr(a,z)
    real :: z
    real(dp) :: a
    hypGeo1F0Deriv_dr = hypGeo1F0DerivFunc(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F0Deriv_id(a,z)
    integer :: a
    real(dp) :: z
    hypGeo1F0Deriv_id = hypGeo1F0DerivFunc(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F0Deriv_rd(a,z)
    real :: a
    real(dp) :: z
    hypGeo1F0Deriv_rd = hypGeo1F0DerivFunc(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F0Deriv_dd(a,z)
    real(dp) :: a, z
    hypGeo1F0Deriv_dd = hypGeo1F0DerivFunc(cmplx(a,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F0Deriv_iCr(a,z)
    integer :: a
    complex :: z
    hypGeo1F0Deriv_iCr = hypGeo1F0DerivFunc(cmplx(a,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo1F0Deriv_rCr(a,z)
    real :: a
    complex :: z
    hypGeo1F0Deriv_rCr = hypGeo1F0DerivFunc(cmplx(a,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo1F0Deriv_dCr(a,z)
    real(dp) :: a
    complex :: z
    hypGeo1F0Deriv_dCr = hypGeo1F0DerivFunc(cmplx(a,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo1F0Deriv_iCd(a,z)
    integer :: a
    complex(dp) :: z
    hypGeo1F0Deriv_iCd = hypGeo1F0DerivFunc(cmplx(a,0.d0,dp),z)
  end function

  complex(dp) function hypGeo1F0Deriv_rCd(a,z)
    real :: a
    complex(dp) :: z
    hypGeo1F0Deriv_rCd = hypGeo1F0DerivFunc(cmplx(a,0.d0,dp),z)
  end function

  complex(dp) function hypGeo1F0Deriv_dCd(a,z)
    real(dp) :: a
    complex(dp) :: z
    hypGeo1F0Deriv_dCd = hypGeo1F0DerivFunc(cmplx(a,0.d0,dp),z)
  end function

  complex(dp) function hypGeo1F0Deriv_CrCr(a,z)
    complex :: a, z
    hypGeo1F0Deriv_CrCr = hypGeo1F0DerivFunc(cmplx(real(a),aimag(a),dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo1F0Deriv_CrCd(a,z)
    complex :: a
    complex(dp) :: z
    hypGeo1F0Deriv_CrCd = hypGeo1F0DerivFunc(cmplx(real(a),aimag(a),dp),z)
  end function

  complex(dp) function hypGeo1F0Deriv_CdCd(a,z)
    complex(dp) :: a, z
    hypGeo1F0Deriv_CdCd = hypGeo1F0DerivFunc(a,z)
  end function

  !*****************************!
  ! Hypergeometric Function 1F1 !
  !*****************************!

  function hypGeo1F1_odeFunc(nF,nC,c,x,y)
    integer :: nF, nC
    real(dp) :: hypGeo1F1_odeFunc(nF), x, y(nF)
    complex(dp) :: c(nC), zs, fT, fpT, f, fp, dz
    dz = c(2)-c(1)
    zs = c(1) + x*dz
    fT = cmplx(y(1),y(2),dp); fpT = cmplx(y(3),y(4),dp)
    f = dz*fpT
    fp = (c(3)*fT-(c(4)-zs)*fpT)*dz/zs
    hypGeo1F1_odeFunc(1) = real(f); hypGeo1F1_odeFunc(2) = aimag(f)
    hypGeo1F1_odeFunc(3) = real(fp); hypGeo1F1_odeFunc(4) = aimag(fp)
  end function

  complex(dp) function hypGeo1F1_series(a,b,z)
    !! Returns Hypergeometric Function 1F1
    !! Should be only used for |z|<1.0
    integer :: i
    complex(dp) :: a, b, z, result, oldResult, indvResult
    result = cmplx(0.d0,0.d0,dp)
    oldResult = result
    do i=0, 100000
      indvResult = (pochhammerR(a,i)/pochhammerR(b,i))*(z**i)/gamma(dble(i)+1.d0)
      if (indvResult /= indvResult) exit
      result = result + indvResult
      if (abs(result-oldResult).le.1.d-16) exit
      oldResult = result
    end do
  hypGeo1F1_series = result
  end function

  complex(dp) function hypGeo1F1Deriv_series(a,b,z)
    !! Returns Derivative of the Hypergeometric Function 1F1
    !! Should be only used for |z|<1.0
    complex(dp) :: a, b, z
    hypGeo1F1Deriv_series = (a/b)*hypGeo1F1_series(a+1.d0,b+1.d0,z)
  end function

  complex(dp) function hypGeo1F1Func(a,b,z)
    complex(dp) :: a, b, z
    real(dp) :: y0(4), odeResult(4)
    complex(dp) :: z0, series(2), consts(4)
    if(abs(z).lt.1.d0) then
      hypGeo1F1Func = hypGeo1F1_series(a,b,z)
      return
    else if(real(z).gt.0.d0 .and. real(z).lt.1.d0) then
      z0 = cmplx(0.5d0,0.d0,dp)
    else if(real(z).gt.-1.d0 .and. real(z).lt.0.d0) then
      z0 = cmplx(-0.5d0,0.d0,dp)
    else if(real(z).gt.1.d0) then
      z0 = cmplx(0.d0,0.5d0,dp)
    else if(real(z).lt.-1.d0) then
      z0 = cmplx(0.d0,-0.5d0,dp)
    end if
    series(1) = hypGeo1F1_series(a,b,z0)
    series(2) = hypGeo1F1Deriv_series(a,b,z0)
    consts(1) = z0; consts(2) = z
    consts(3) = a; consts(4) = b;
    y0(1) = real(series(1)); y0(2) = aimag(series(1))
    y0(3) = real(series(2)); y0(4) = aimag(series(2))
    odeResult = rk1AdaptStepCmplxC(hypGeo1F1_odeFunc,4,0.d0,1.d0,y0,4,consts)
    hypGeo1F1Func = cmplx(odeResult(1),odeResult(2),dp)
  end function

  complex(dp) function hypGeo1F1DerivFunc(a,b,z)
    complex(dp) :: a, b, z
    hypGeo1F1DerivFunc = (a/b)*hypGeo1F1Func(a+1.d0,b+1.d0,z)
  end function

  complex(dp) function hypGeo1F1_iii(a,b,z)
    integer :: a, b, z
    hypGeo1F1_iii = hypGeo1F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F1_rri(a,b,z)
    integer :: z
    real :: a, b
    hypGeo1F1_rri = hypGeo1F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F1_ddi(a,b,z)
    integer :: z
    real(dp) :: a, b
    hypGeo1F1_ddi = hypGeo1F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F1_iir(a,b,z)
    integer :: a, b
    real :: z
    hypGeo1F1_iir = hypGeo1F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F1_rrr(a,b,z)
    real :: a, b, z
    hypGeo1F1_rrr = hypGeo1F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F1_ddr(a,b,z)
    real :: z
    real(dp) :: a, b
    hypGeo1F1_ddr = hypGeo1F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F1_iid(a,b,z)
    integer :: a, b
    real(dp) :: z
    hypGeo1F1_iid = hypGeo1F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F1_rrd(a,b,z)
    real :: a, b
    real(dp) :: z
    hypGeo1F1_rrd = hypGeo1F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F1_ddd(a,b,z)
    real(dp) :: a, b, z
    hypGeo1F1_ddd = hypGeo1F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F1_iiCr(a,b,z)
    integer :: a, b
    complex :: z
    hypGeo1F1_iiCr = hypGeo1F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo1F1_rrCr(a,b,z)
    real :: a, b
    complex :: z
    hypGeo1F1_rrCr = hypGeo1F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo1F1_ddCr(a,b,z)
    real(dp) :: a, b
    complex :: z
    hypGeo1F1_ddCr = hypGeo1F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo1F1_iiCd(a,b,z)
    integer :: a, b
    complex(dp) :: z
    hypGeo1F1_iiCd = hypGeo1F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),z)
  end function

  complex(dp) function hypGeo1F1_rrCd(a,b,z)
    real :: a, b
    complex(dp) :: z
    hypGeo1F1_rrCd = hypGeo1F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),z)
  end function

  complex(dp) function hypGeo1F1_ddCd(a,b,z)
    real(dp) :: a, b
    complex(dp) :: z
    hypGeo1F1_ddCd = hypGeo1F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),z)
  end function

  complex(dp) function hypGeo1F1_CrCrCr(a,b,z)
    complex :: a, b, z
    hypGeo1F1_CrCrCr = hypGeo1F1Func(cmplx(real(a),aimag(a),dp),cmplx(real(b),aimag(b),dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo1F1_CrCrCd(a,b,z)
    complex :: a, b
    complex(dp) :: z
    hypGeo1F1_CrCrCd = hypGeo1F1Func(cmplx(real(a),aimag(a),dp),cmplx(real(b),aimag(b),dp),z)
  end function

  complex(dp) function hypGeo1F1_CdCdCd(a,b,z)
    complex(dp) :: a, b, z
    hypGeo1F1_CdCdCd = hypGeo1F1Func(a,b,z)
  end function

  complex(dp) function hypGeo1F1Deriv_iii(a,b,z)
    integer :: a, b, z
    hypGeo1F1Deriv_iii = hypGeo1F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F1Deriv_rri(a,b,z)
    integer :: z
    real :: a, b
    hypGeo1F1Deriv_rri = hypGeo1F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F1Deriv_ddi(a,b,z)
    integer :: z
    real(dp) :: a, b
    hypGeo1F1Deriv_ddi = hypGeo1F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F1Deriv_iir(a,b,z)
    integer :: a, b
    real :: z
    hypGeo1F1Deriv_iir = hypGeo1F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F1Deriv_rrr(a,b,z)
    real :: a, b, z
    hypGeo1F1Deriv_rrr = hypGeo1F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F1Deriv_ddr(a,b,z)
    real :: z
    real(dp) :: a, b
    hypGeo1F1Deriv_ddr = hypGeo1F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F1Deriv_iid(a,b,z)
    integer :: a, b
    real(dp) :: z
    hypGeo1F1Deriv_iid = hypGeo1F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F1Deriv_rrd(a,b,z)
    real :: a, b
    real(dp) :: z
    hypGeo1F1Deriv_rrd = hypGeo1F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F1Deriv_ddd(a,b,z)
    real(dp) :: a, b, z
    hypGeo1F1Deriv_ddd = hypGeo1F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo1F1Deriv_iiCr(a,b,z)
    integer :: a, b
    complex :: z
    hypGeo1F1Deriv_iiCr = hypGeo1F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo1F1Deriv_rrCr(a,b,z)
    real :: a, b
    complex :: z
    hypGeo1F1Deriv_rrCr = hypGeo1F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo1F1Deriv_ddCr(a,b,z)
    real(dp) :: a, b
    complex :: z
    hypGeo1F1Deriv_ddCr = hypGeo1F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo1F1Deriv_iiCd(a,b,z)
    integer :: a, b
    complex(dp) :: z
    hypGeo1F1Deriv_iiCd = hypGeo1F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),z)
  end function

  complex(dp) function hypGeo1F1Deriv_rrCd(a,b,z)
    real :: a, b
    complex(dp) :: z
    hypGeo1F1Deriv_rrCd = hypGeo1F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),z)
  end function

  complex(dp) function hypGeo1F1Deriv_ddCd(a,b,z)
    real(dp) :: a, b
    complex(dp) :: z
    hypGeo1F1Deriv_ddCd = hypGeo1F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),z)
  end function

  complex(dp) function hypGeo1F1Deriv_CrCrCr(a,b,z)
    complex :: a, b, z
    hypGeo1F1Deriv_CrCrCr = hypGeo1F1DerivFunc(cmplx(real(a),aimag(a),dp),cmplx(real(b),aimag(b),dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo1F1Deriv_CrCrCd(a,b,z)
    complex :: a, b
    complex(dp) :: z
    hypGeo1F1Deriv_CrCrCd = hypGeo1F1DerivFunc(cmplx(real(a),aimag(a),dp),cmplx(real(b),aimag(b),dp),z)
  end function

  complex(dp) function hypGeo1F1Deriv_CdCdCd(a,b,z)
    complex(dp) :: a, b, z
    hypGeo1F1Deriv_CdCdCd = hypGeo1F1DerivFunc(a,b,z)
  end function


  !*****************************!
  ! Hypergeometric Function 2F1 !
  !*****************************!

  function hypGeo2F1_odeFunc(nF,nC,c,x,y)
    integer :: nF, nC
    real(dp) :: hypGeo2F1_odeFunc(nF), x, y(nF)
    complex(dp) :: c(nC), zs, fT, fpT, f, fp, dz
    dz = c(2)-c(1)
    zs = c(1) + x*dz
    fT = cmplx(y(1),y(2),dp); fpT = cmplx(y(3),y(4),dp)
    f = dz*fpT
    fp = (c(3)*c(4)*fT-(c(5)-(c(3)+c(4)+1.d0)*zs)*fpT)*dz/(zs*(1.d0-zs))
    hypGeo2F1_odeFunc(1) = real(f); hypGeo2F1_odeFunc(2) = aimag(f)
    hypGeo2F1_odeFunc(3) = real(fp); hypGeo2F1_odeFunc(4) = aimag(fp)
  end function

  complex(dp) function hypGeo2F1_series(a,b,c,z)
    !! Returns Hypergeometric Function 2F1
    !! Should be only used for |z|<1.0
    integer :: i
    complex(dp) :: a, b, c, z, result, oldResult, indvResult
    result = cmplx(0.d0,0.d0,dp)
    oldResult = result
    do i=0, 100000
      indvResult = ((pochhammerR(a,i)*pochhammerR(b,i))/pochhammerR(c,i))*(z**i)/gamma(dble(i)+1.d0)
      if (indvResult /= indvResult) exit
      result = result + indvResult
      if (abs(result-oldResult).le.1.d-16) exit
      oldResult = result
    end do
  hypGeo2F1_series = result
  end function

  complex(dp) function hypGeo2F1Deriv_series(a,b,c,z)
    !! Returns Derivative of the Hypergeometric Function 2F1
    !! Should be only used for |z|<1.0
    complex(dp) :: a, b, c, z
    hypGeo2F1Deriv_series = (a*b/c)*hypGeo2F1_series(a+1.d0,b+1.d0,c+1.d0,z)
  end function

  complex(dp) function hypGeo2F1Func(a,b,c,z)
    complex(dp) :: a, b, c, z
    real(dp) :: y0(4), odeResult(4)
    complex(dp) :: z0, series(2), consts(5)
    if(abs(z).lt.1.d0) then
      hypGeo2F1Func = hypGeo2F1_series(a,b,c,z)
      return
    else if(real(z).gt.0.d0 .and. real(z).lt.1.d0) then
      z0 = cmplx(0.5d0,0.d0,dp)
    else if(real(z).gt.-1.d0 .and. real(z).lt.0.d0) then
      z0 = cmplx(-0.5d0,0.d0,dp)
    else if(real(z).gt.1.d0) then
      z0 = cmplx(0.d0,0.5d0,dp)
    else if(real(z).lt.-1.d0) then
      z0 = cmplx(0.d0,-0.5d0,dp)
    end if
    series(1) = hypGeo2F1_series(a,b,c,z0)
    series(2) = hypGeo2F1Deriv_series(a,b,c,z0)
    consts(1) = z0; consts(2) = z
    consts(3) = a; consts(4) = b; consts(5) = c
    y0(1) = real(series(1)); y0(2) = aimag(series(1))
    y0(3) = real(series(2)); y0(4) = aimag(series(2))
    odeResult = rk1AdaptStepCmplxC(hypGeo2F1_odeFunc,4,0.d0,1.d0,y0,5,consts)
    hypGeo2F1Func = cmplx(odeResult(1),odeResult(2),dp)
  end function

  complex(dp) function hypGeo2F1DerivFunc(a,b,c,z)
    complex(dp) :: a, b, c, z
    hypGeo2F1DerivFunc = (a*b/c)*hypGeo2F1Func(a+1.d0,b+1.d0,c+1.d0,z)
  end function

  complex(dp) function hypGeo2F1_iiii(a,b,c,z)
    integer :: a, b, c, z
    hypGeo2F1_iiii = hypGeo2F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo2F1_rrri(a,b,c,z)
    real :: a, b, c
    integer :: z
    hypGeo2F1_rrri = hypGeo2F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo2F1_dddi(a,b,c,z)
    real(dp) :: a, b, c
    integer :: z
    hypGeo2F1_dddi = hypGeo2F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo2F1_iiir(a,b,c,z)
    integer :: a, b, c
    real :: z
    hypGeo2F1_iiir = hypGeo2F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo2F1_rrrr(a,b,c,z)
    real :: a, b, c, z
    hypGeo2F1_rrrr = hypGeo2F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo2F1_dddr(a,b,c,z)
    real(dp) :: a, b, c
    real :: z
    hypGeo2F1_dddr = hypGeo2F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo2F1_iiid(a,b,c,z)
    integer :: a, b, c
    real(dp) :: z
    hypGeo2F1_iiid = hypGeo2F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo2F1_rrrd(a,b,c,z)
    real :: a, b, c
    real(dp) :: z
    hypGeo2F1_rrrd = hypGeo2F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo2F1_dddd(a,b,c,z)
    real(dp) :: a, b, c, z
    hypGeo2F1_dddd = hypGeo2F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo2F1_iiiCr(a,b,c,z)
    integer :: a, b, c
    complex :: z
    hypGeo2F1_iiiCr = hypGeo2F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo2F1_rrrCr(a,b,c,z)
    real :: a, b, c
    complex :: z
    hypGeo2F1_rrrCr = hypGeo2F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo2F1_dddCr(a,b,c,z)
    real(dp) :: a, b, c
    complex :: z
    hypGeo2F1_dddCr = hypGeo2F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo2F1_iiiCd(a,b,c,z)
    integer :: a, b, c
    complex(dp) :: z
    hypGeo2F1_iiiCd = hypGeo2F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),z)
  end function

  complex(dp) function hypGeo2F1_rrrCd(a,b,c,z)
    real :: a, b, c
    complex(dp) :: z
    hypGeo2F1_rrrCd = hypGeo2F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),z)
  end function

  complex(dp) function hypGeo2F1_dddCd(a,b,c,z)
    real(dp) :: a, b, c
    complex(dp) :: z
    hypGeo2F1_dddCd = hypGeo2F1Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),z)
  end function

  complex(dp) function hypGeo2F1_CrCrCrCr(a,b,c,z)
    complex :: a, b, c, z
    hypGeo2F1_CrCrCrCr = hypGeo2F1Func(cmplx(real(a),aimag(a),dp),cmplx(real(b),aimag(b),dp),cmplx(real(c),aimag(c),dp),&
      cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo2F1_CrCrCrCd(a,b,c,z)
    complex :: a, b, c
    complex(dp) :: z
    hypGeo2F1_CrCrCrCd = hypGeo2F1Func(cmplx(real(a),aimag(a),dp),cmplx(real(b),aimag(b),dp),cmplx(real(c),aimag(c),dp),z)
  end function

  complex(dp) function hypGeo2F1_CdCdCdCd(a,b,c,z)
    complex(dp) :: a, b, c, z
    hypGeo2F1_CdCdCdCd = hypGeo2F1Func(a,b,c,z)
  end function

  complex(dp) function hypGeo2F1Deriv_iiii(a,b,c,z)
    integer :: a, b, c, z
    hypGeo2F1Deriv_iiii = hypGeo2F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo2F1Deriv_rrri(a,b,c,z)
    real :: a, b, c
    integer :: z
    hypGeo2F1Deriv_rrri = hypGeo2F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo2F1Deriv_dddi(a,b,c,z)
    real(dp) :: a, b, c
    integer :: z
    hypGeo2F1Deriv_dddi = hypGeo2F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo2F1Deriv_iiir(a,b,c,z)
    integer :: a, b, c
    real :: z
    hypGeo2F1Deriv_iiir = hypGeo2F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo2F1Deriv_rrrr(a,b,c,z)
    real :: a, b, c, z
    hypGeo2F1Deriv_rrrr = hypGeo2F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo2F1Deriv_dddr(a,b,c,z)
    real(dp) :: a, b, c
    real :: z
    hypGeo2F1Deriv_dddr = hypGeo2F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo2F1Deriv_iiid(a,b,c,z)
    integer :: a, b, c
    real(dp) :: z
    hypGeo2F1Deriv_iiid = hypGeo2F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo2F1Deriv_rrrd(a,b,c,z)
    real :: a, b, c
    real(dp) :: z
    hypGeo2F1Deriv_rrrd = hypGeo2F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo2F1Deriv_dddd(a,b,c,z)
    real(dp) :: a, b, c, z
    hypGeo2F1Deriv_dddd = hypGeo2F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo2F1Deriv_iiiCr(a,b,c,z)
    integer :: a, b, c
    complex :: z
    hypGeo2F1Deriv_iiiCr = hypGeo2F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo2F1Deriv_rrrCr(a,b,c,z)
    real :: a, b, c
    complex :: z
    hypGeo2F1Deriv_rrrCr = hypGeo2F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo2F1Deriv_dddCr(a,b,c,z)
    real(dp) :: a, b, c
    complex :: z
    hypGeo2F1Deriv_dddCr = hypGeo2F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo2F1Deriv_iiiCd(a,b,c,z)
    integer :: a, b, c
    complex(dp) :: z
    hypGeo2F1Deriv_iiiCd = hypGeo2F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),z)
  end function

  complex(dp) function hypGeo2F1Deriv_rrrCd(a,b,c,z)
    real :: a, b, c
    complex(dp) :: z
    hypGeo2F1Deriv_rrrCd = hypGeo2F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),z)
  end function

  complex(dp) function hypGeo2F1Deriv_dddCd(a,b,c,z)
    real(dp) :: a, b, c
    complex(dp) :: z
    hypGeo2F1Deriv_dddCd = hypGeo2F1DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp),z)
  end function

  complex(dp) function hypGeo2F1Deriv_CrCrCrCr(a,b,c,z)
    complex :: a, b, c, z
    hypGeo2F1Deriv_CrCrCrCr = hypGeo2F1DerivFunc(cmplx(real(a),aimag(a),dp),cmplx(real(b),aimag(b),dp),cmplx(real(c),aimag(c),dp),&
      cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo2F1Deriv_CrCrCrCd(a,b,c,z)
    complex :: a, b, c
    complex(dp) :: z
    hypGeo2F1Deriv_CrCrCrCd = hypGeo2F1DerivFunc(cmplx(real(a),aimag(a),dp),cmplx(real(b),aimag(b),dp),cmplx(real(c),aimag(c),dp),z)
  end function

  complex(dp) function hypGeo2F1Deriv_CdCdCdCd(a,b,c,z)
    complex(dp) :: a, b, c, z
    hypGeo2F1Deriv_CdCdCdCd = hypGeo2F1DerivFunc(a,b,c,z)
  end function

  !*****************************!
  ! Hypergeometric Function 3F2 !
  !*****************************!

  function hypGeo3F2_odeFunc(nF,nC,c,x,y)
    integer :: nF, nC
    real(dp) :: hypGeo3F2_odeFunc(nF), x, y(nF)
    complex(dp) :: c(nC), zs, fT, fpT, fppT, f, fp, fpp, dz
    complex(dp) :: factfT, factfpT, factfpT2, factfppT
    dz = c(2)-c(1)
    zs = c(1) + x*dz
    factfT = -c(3)*c(4)*c(5)
    factfpT = c(4)*c(3)+c(5)*c(3)+c(3)+c(4)+c(4)*c(5)+c(5)+1.d0
    factfpT2 = c(6)*c(7) - factfpT*zs
    factfppT = (-(c(3)+c(4)+c(5)+3.d0)*zs+c(6)+c(7)+1.d0)*zs
    fT = cmplx(y(1),y(2),dp); fpT = cmplx(y(3),y(4),dp); fppT = cmplx(y(5),y(6),dp)
    f = dz*fpT
    fp = dz*fppT
    fpp = -dz*(factfppT*fppT+factfpT2*fpT+factfT*fT)/(zs*zs*(1.d0-zs))
    hypGeo3F2_odeFunc(1) = real(f); hypGeo3F2_odeFunc(2) = aimag(f)
    hypGeo3F2_odeFunc(3) = real(fp); hypGeo3F2_odeFunc(4) = aimag(fp)
    hypGeo3F2_odeFunc(5) = real(fpp); hypGeo3F2_odeFunc(6) = aimag(fpp)
  end function

  complex(dp) function hypGeo3F2_series(a,b,c,d,e,z)
    !! Returns Hypergeometric Function 3F2
    !! Should be only used for |z|<1.0
    integer :: i
    complex(dp) :: a, b, c, d, e, z, result, oldResult, indvResult
    result = cmplx(0.d0,0.d0,dp)
    oldResult = result
    do i=0, 100000
      indvResult = ((pochhammerR(a,i)*pochhammerR(b,i)*pochhammerR(c,i))/(pochhammerR(d,i)*pochhammerR(e,i)))* &
        (z**i)/gamma(dble(i)+1.d0)
      if (indvResult /= indvResult) exit
      result = result + indvResult
      if (abs(result-oldResult).le.1.d-16) exit
      oldResult = result
    end do
    hypGeo3F2_series = result
  end function

  complex(dp) function hypGeo3F2Deriv_series(a,b,c,d,e,z)
    !! Returns Hypergeometric Function 3F2
    !! Should be only used for |z|<1.0
    complex(dp) :: a, b, c, d, e, z
    hypGeo3F2Deriv_series = (a*b*c/(d*e))*hypGeo3F2_series(a+1.d0,b+1.d0,c+1.d0,d+1.d0,e+1.d0,z)
  end function

  complex(dp) function hypGeo3F2Func(a,b,c,d,e,z)
    complex(dp) :: a, b, c, d, e, z
    real(dp) :: y0(6), odeResult(6)
    complex(dp) :: z0, series(3), consts(7)
    if(abs(z).lt.1.d0) then
      hypGeo3F2Func = hypGeo3F2_series(a,b,c,d,e,z)
      return
    else if(real(z).gt.0.d0 .and. real(z).lt.1.d0) then
      z0 = cmplx(0.5d0,0.d0,dp)
    else if(real(z).gt.-1.d0 .and. real(z).lt.0.d0) then
      z0 = cmplx(-0.5d0,0.d0,dp)
    else if(real(z).gt.1.d0) then
      z0 = cmplx(0.d0,0.5d0,dp)
    else if(real(z).lt.-1.d0) then
      z0 = cmplx(0.d0,-0.5d0,dp)
    end if
    series(1) = hypGeo3F2_series(a,b,c,d,e,z0)
    series(2) = hypGeo3F2Deriv_series(a,b,c,d,e,z0)
    series(3) = (a*(a+1.d0)*b*(b+1.d0)*c*(c+1.d0))/(d*(d+1.d0)*e*(e+1.d0))*&
      hypGeo3F2_series(a+2.d0,b+2.d0,c+2.d0,d+2.d0,e+2.d0,z0)
    consts(1) = z0; consts(2) = z
    consts(3) = a; consts(4) = b; consts(5) = c; consts(6) = d; consts(7) = e
    y0(1) = real(series(1)); y0(2) = aimag(series(1))
    y0(3) = real(series(2)); y0(4) = aimag(series(2))
    y0(5) = real(series(3)); y0(6) = aimag(series(3))
    odeResult = rk1AdaptStepCmplxC(hypGeo3F2_odeFunc,6,0.d0,1.d0,y0,7,consts)
    hypGeo3F2Func = cmplx(odeResult(1),odeResult(2),dp)
  end function

  complex(dp) function hypGeo3F2DerivFunc(a,b,c,d,e,z)
    complex(dp) :: a, b, c, d, e, z
    hypGeo3F2DerivFunc = (a*b*c/(d*e))*hypGeo3F2Func(a+1.d0,b+1.d0,c+1.d0,d+1.d0,e+1.d0,z)
  end function

  complex(dp) function hypGeo3F2_iiiiii(a,b,c,d,e,z)
    integer :: a, b, c, d, e, z
    hypGeo3F2_iiiiii = hypGeo3F2Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo3F2_iiiiir(a,b,c,d,e,z)
    integer :: a, b, c, d, e
    real :: z
    hypGeo3F2_iiiiir = hypGeo3F2Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo3F2_iiiiid(a,b,c,d,e,z)
    integer :: a, b, c, d, e
    real(dp) :: z
    hypGeo3F2_iiiiid = hypGeo3F2Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo3F2_rrrrri(a,b,c,d,e,z)
    real :: a, b, c, d, e
    integer :: z
    hypGeo3F2_rrrrri = hypGeo3F2Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo3F2_rrrrrr(a,b,c,d,e,z)
    real :: a, b, c, d, e, z
    hypGeo3F2_rrrrrr = hypGeo3F2Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo3F2_rrrrrd(a,b,c,d,e,z)
    real :: a, b, c, d, e
    real(dp) :: z
    hypGeo3F2_rrrrrd = hypGeo3F2Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo3F2_dddddi(a,b,c,d,e,z)
    real(dp) :: a, b, c, d, e
    integer :: z
    hypGeo3F2_dddddi = hypGeo3F2Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo3F2_dddddr(a,b,c,d,e,z)
    real(dp) :: a, b, c, d, e
    real :: z
    hypGeo3F2_dddddr = hypGeo3F2Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo3F2_dddddd(a,b,c,d,e,z)
    real(dp) :: a, b, c, d, e, z
    hypGeo3F2_dddddd = hypGeo3F2Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo3F2_iiiiiCr(a,b,c,d,e,z)
    integer :: a, b, c, d, e
    complex :: z
    hypGeo3F2_iiiiiCr = hypGeo3F2Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo3F2_iiiiiCd(a,b,c,d,e,z)
    integer :: a, b, c, d, e
    complex(dp) :: z
    hypGeo3F2_iiiiiCd = hypGeo3F2Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),z)
  end function

  complex(dp) function hypGeo3F2_rrrrrCr(a,b,c,d,e,z)
    real :: a, b, c, d, e
    complex :: z
    hypGeo3F2_rrrrrCr = hypGeo3F2Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo3F2_rrrrrCd(a,b,c,d,e,z)
    real :: a, b, c, d, e
    complex(dp) :: z
    hypGeo3F2_rrrrrCd = hypGeo3F2Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),z)
  end function

  complex(dp) function hypGeo3F2_dddddCr(a,b,c,d,e,z)
    real(dp) :: a, b, c, d, e
    complex :: z
    hypGeo3F2_dddddCr = hypGeo3F2Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo3F2_dddddCd(a,b,c,d,e,z)
    real(dp) :: a, b, c, d, e
    complex(dp) :: z
    hypGeo3F2_dddddCd = hypGeo3F2Func(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),z)
  end function

  complex(dp) function hypGeo3F2_CrCrCrCrCrCr(a,b,c,d,e,z)
    complex :: a, b, c, d, e, z
    hypGeo3F2_CrCrCrCrCrCr = hypGeo3F2Func(cmplx(real(a),aimag(a),dp),cmplx(real(b),aimag(b),dp),&
      cmplx(real(c),aimag(c),dp), cmplx(real(d),aimag(d),dp),cmplx(real(e),aimag(e),dp), &
      cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo3F2_CrCrCrCrCrCd(a,b,c,d,e,z)
    complex :: a, b, c, d, e
    complex(dp) :: z
    hypGeo3F2_CrCrCrCrCrCd = hypGeo3F2Func(cmplx(real(a),aimag(a),dp),cmplx(real(b),aimag(b),dp),&
      cmplx(real(c),aimag(c),dp), cmplx(real(d),aimag(d),dp),cmplx(real(e),aimag(e),dp),z)
  end function

  complex(dp) function hypGeo3F2_CdCdCdCdCdCd(a,b,c,d,e,z)
    complex(dp) :: a, b, c, d, e, z
    hypGeo3F2_CdCdCdCdCdCd = hypGeo3F2Func(a,b,c,d,e,z)
  end function

  complex(dp) function hypGeo3F2Deriv_iiiiii(a,b,c,d,e,z)
    integer :: a, b, c, d, e, z
    hypGeo3F2Deriv_iiiiii = hypGeo3F2DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo3F2Deriv_iiiiir(a,b,c,d,e,z)
    integer :: a, b, c, d, e
    real :: z
    hypGeo3F2Deriv_iiiiir = hypGeo3F2DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo3F2Deriv_iiiiid(a,b,c,d,e,z)
    integer :: a, b, c, d, e
    real(dp) :: z
    hypGeo3F2Deriv_iiiiid = hypGeo3F2DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo3F2Deriv_rrrrri(a,b,c,d,e,z)
    real :: a, b, c, d, e
    integer :: z
    hypGeo3F2Deriv_rrrrri = hypGeo3F2DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo3F2Deriv_rrrrrr(a,b,c,d,e,z)
    real :: a, b, c, d, e, z
    hypGeo3F2Deriv_rrrrrr = hypGeo3F2DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo3F2Deriv_rrrrrd(a,b,c,d,e,z)
    real :: a, b, c, d, e
    real(dp) :: z
    hypGeo3F2Deriv_rrrrrd = hypGeo3F2DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo3F2Deriv_dddddi(a,b,c,d,e,z)
    real(dp) :: a, b, c, d, e
    integer :: z
    hypGeo3F2Deriv_dddddi = hypGeo3F2DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo3F2Deriv_dddddr(a,b,c,d,e,z)
    real(dp) :: a, b, c, d, e
    real :: z
    hypGeo3F2Deriv_dddddr = hypGeo3F2DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo3F2Deriv_dddddd(a,b,c,d,e,z)
    real(dp) :: a, b, c, d, e, z
    hypGeo3F2Deriv_dddddd = hypGeo3F2DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(z,0.d0,dp))
  end function

  complex(dp) function hypGeo3F2Deriv_iiiiiCr(a,b,c,d,e,z)
    integer :: a, b, c, d, e
    complex :: z
    hypGeo3F2Deriv_iiiiiCr = hypGeo3F2DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo3F2Deriv_iiiiiCd(a,b,c,d,e,z)
    integer :: a, b, c, d, e
    complex(dp) :: z
    hypGeo3F2Deriv_iiiiiCd = hypGeo3F2DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),z)
  end function

  complex(dp) function hypGeo3F2Deriv_rrrrrCr(a,b,c,d,e,z)
    real :: a, b, c, d, e
    complex :: z
    hypGeo3F2Deriv_rrrrrCr = hypGeo3F2DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo3F2Deriv_rrrrrCd(a,b,c,d,e,z)
    real :: a, b, c, d, e
    complex(dp) :: z
    hypGeo3F2Deriv_rrrrrCd = hypGeo3F2DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),z)
  end function

  complex(dp) function hypGeo3F2Deriv_dddddCr(a,b,c,d,e,z)
    real(dp) :: a, b, c, d, e
    complex :: z
    hypGeo3F2Deriv_dddddCr = hypGeo3F2DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo3F2Deriv_dddddCd(a,b,c,d,e,z)
    real(dp) :: a, b, c, d, e
    complex(dp) :: z
    hypGeo3F2Deriv_dddddCd = hypGeo3F2DerivFunc(cmplx(a,0.d0,dp),cmplx(b,0.d0,dp),cmplx(c,0.d0,dp), &
      cmplx(d,0.d0,dp),cmplx(e,0.d0,dp),z)
  end function

  complex(dp) function hypGeo3F2Deriv_CrCrCrCrCrCr(a,b,c,d,e,z)
    complex :: a, b, c, d, e, z
    hypGeo3F2Deriv_CrCrCrCrCrCr = hypGeo3F2DerivFunc(cmplx(real(a),aimag(a),dp),cmplx(real(b),aimag(b),dp),&
      cmplx(real(c),aimag(c),dp), cmplx(real(d),aimag(d),dp),cmplx(real(e),aimag(e),dp), &
      cmplx(real(z),aimag(z),dp))
  end function

  complex(dp) function hypGeo3F2Deriv_CrCrCrCrCrCd(a,b,c,d,e,z)
    complex :: a, b, c, d, e
    complex(dp) :: z
    hypGeo3F2Deriv_CrCrCrCrCrCd = hypGeo3F2DerivFunc(cmplx(real(a),aimag(a),dp),cmplx(real(b),aimag(b),dp),&
      cmplx(real(c),aimag(c),dp), cmplx(real(d),aimag(d),dp),cmplx(real(e),aimag(e),dp),z)
  end function

  complex(dp) function hypGeo3F2Deriv_CdCdCdCdCdCd(a,b,c,d,e,z)
    complex(dp) :: a, b, c, d, e, z
    hypGeo3F2Deriv_CdCdCdCdCdCd = hypGeo3F2DerivFunc(a,b,c,d,e,z)
  end function

  !*************!
  ! Dilogarithm !
  !*************!

  complex(dp) function diLog_series(x)
    integer :: i
    complex(dp) :: x, result, oldResult, indvResult
    result = 0.d0
    oldResult = result
    do i=1, 1000000
      indvResult = (x**i)/(i**2)
      if (indvResult /= indvResult) exit
      result = result + indvResult
      if (abs(result-oldResult).le.1.d-16) exit
      oldResult = result
    end do
    diLog_series = result
  end function

  complex(dp) function diLogFunc(x)
    complex(dp) :: x
    if(abs(x).lt.1.d0) then
      diLogFunc = diLog_series(x)
      return
    else
      diLogFunc = x*hypGeo3F2(1.d0,1.d0,1.d0,2.d0,2.d0,x)
      return
    end if
  end function

  complex(dp) function diLog_i(x)
    integer :: x
    diLog_i = diLogFunc(cmplx(x,0.d0,dp))
  end function

  complex(dp) function diLog_r(x)
    real :: x
    diLog_r = diLogFunc(cmplx(x,0.d0,dp))
  end function

  complex(dp) function diLog_d(x)
    real(dp) :: x
    diLog_d = diLogFunc(cmplx(x,0.d0,dp))
  end function

  complex(dp) function diLog_Cr(x)
    complex :: x
    diLog_Cr = diLogFunc(cmplx(real(x),aimag(x),dp))
  end function

  complex(dp) function diLog_Cd(x)
    complex(dp) :: x
    diLog_Cd = diLogFunc(x)
  end function


end module
