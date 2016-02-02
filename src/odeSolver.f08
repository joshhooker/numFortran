!*********************************************************
!>
! Functions to solve ODES:
!
!  * 1st order solvers:
!    * rk1FixedStep
!    * rk1AdaptStep

module odeSolver
  implicit none

  private
  integer, parameter :: dp = kind(0.d0)

  public :: rk1FixedStep, rk1AdaptStep
  public :: rk1AdaptStepCmplxC
contains

  function rk1FixedStep(f,nF,a,b,h,y0,nC,consts)
    !! 4th order Runga Kutta ODE Solver with fixed step h and initial
    !! value \(y(a) = y_0\). If the fixed step does not evenly divide
    !! into \(b-a\), a smaller step is placed at the end so the final
    !! \(x\) position is at \(b\)

    integer :: i, nF, nC
    real(dp) :: xn, k1(nF), k2(nF), k3(nF), k4(nF), a, b, h
    real(dp) :: y0(nF), yi(nF), yn(nF), consts(nC)
    real(dp) :: rk1FixedStep(nF)
    interface
      function f(nF,nC,c,x,y)
        integer, parameter :: dp = kind(0.d0)
        integer :: nF, nC
        real(dp) :: c(nC), x, y(nF), f(nF)
      end function f
    end interface
    xn = a
    yn = y0
    do while(xn.lt.b)
      if(xn+h.gt.b) then
        h = b-xn
      end if
      k1 = h*f(nF,nC,consts,xn,yn)
      k2 = h*f(nF,nC,consts,xn+h/2.d0,yn+k1/2.d0)
      k3 = h*f(nF,nC,consts,xn+h/2.d0,yn+k2/2.d0)
      k4 = h*f(nF,nC,consts,xn+h,yn+k3)
      xn = xn+h
      yn = yn+(1.d0/6.d0)*(k1+2.d0*k2+2.d0*k3+k4)
    end do
    rk1FixedStep = yn
  end function

  function rk1AdaptStep(f,nF,a,b,y0,nC,consts)
    !! Runga Kutta ODE Solver with adaptive h size to satisfy a fixed
    !! error and initial value \(y(a) = y_0\). Calculates \(y\) for two
    !! different \(c_i\)'s in \(y = y + \sum_i c_i*k_i\) and compares
    !! them to determine if the step size is too big or small.

    integer :: i, j, nF, nC
    real(dp) :: xn, k1(nF), k2(nF), k3(nF), k4(nF), k5(nF), k6(nF)
    real(dp) :: a, b, h
    real(dp) :: y0(nF), yi(nF), yn(nF), yns(nF), dyn, consts(nC)
    real(dp) :: xni, yni(nF), kci(nF), kcip(nF), ynT(nF), ynsT(nF)
    real(dp) :: ai(6),bi(6,6),ci(6),cip(6)
    real(dp) :: rk1AdaptStep(nF)
    logical :: goodStep
    interface
      function f(nF,nC,c,x,y)
        integer, parameter :: dp = kind(0.d0)
        integer :: nF, nC
        real(dp) :: c(nC), x, y(nF), f(nF)
      end function f
    end interface
    call rk1AdaptStepSetupConstants(ai,bi,ci,cip)
    xn = a
    yn = y0
    h = 0.001d0
    do while(xn<b)
      goodStep = .false.
      if(xn+h>b) then
        h = b-xn
      end if
      do while(goodStep .eqv. .false.)
        k1 = h*f(nF,nC,consts,xn,yn)
        k2 = h*f(nF,nC,consts,xn+ai(2)*h,yn+bi(2,1)*k1)
        k3 = h*f(nF,nC,consts,xn+ai(3)*h,yn+bi(3,1)*k1+bi(3,2)*k2)
        k4 = h*f(nF,nC,consts,xn+ai(4)*h,yn+bi(4,1)*k1+bi(4,2)*k2+bi(4,3)*k3)
        k5 = h*f(nF,nC,consts,xn+ai(5)*h,yn+bi(5,1)*k1+bi(5,2)*k2+bi(5,3)*k3+bi(5,4)*k4)
        k6 = h*f(nF,nC,consts,xn+ai(6)*h,yn+bi(6,1)*k1+bi(6,2)*k2+bi(6,3)*k3+bi(6,4)*k4+bi(6,5)*k5)
        kci = k1*ci(1)+k2*ci(2)+k3*ci(3)+k4*ci(4)+k5*ci(5)+k6*ci(6)
        kcip = k1*cip(1)+k2*cip(2)+k3*cip(3)+k4*cip(4)+k5*cip(5)+k6*cip(6)
        ynT = yn+kci
        ynsT = yn+kcip
        dyn = abs(ynT(1)-ynsT(1))
        if(dyn .gt. 1.d-10) then
          h = h*((1.d-12/dyn)**(0.2d0))
        else
          goodStep = .true.
        end if
      end do
      xn = xn+h
      yns = yn+kcip
      yn = yn+kci
      dyn = abs(yn(1)-yns(1))
      h = h*((1.d-12/dyn)**(0.2d0))
    end do
    rk1AdaptStep = yn
  end function

  function rk1AdaptStepCmplxC(f,nF,a,b,y0,nC,consts)
    !! Runga Kutta ODE Solver with adaptive h size to satisfy a fixed
    !! error and initial value \(y(a) = y_0\). Calculates \(y\) for two
    !! different \(c_i\)'s in \(y = y + \sum_i c_i*k_i\) and compares
    !! them to determine if the step size is too big or small.

    integer :: i, j, nF, nC
    real(dp) :: xn, k1(nF), k2(nF), k3(nF), k4(nF), k5(nF), k6(nF)
    real(dp) :: a, b, h
    real(dp) :: y0(nF), yi(nF), yn(nF), yns(nF), dyn
    real(dp) :: xni, yni(nF), kci(nF), kcip(nF)
    real(dp) :: ynT(nF), ynsT(nF)
    real(dp) :: ai(6),bi(6,6),ci(6),cip(6)
    real(dp) :: rk1AdaptStepCmplxC(nF)
    complex(dp) :: consts(nC)
    logical :: goodStep
    interface
      function f(nF,nC,c,x,y)
        integer, parameter :: dp = kind(0.d0)
        integer :: nF, nC
        real(dp) :: x, y(nF), f(nF)
        complex(dp) :: c(nC)
      end function f
    end interface
    call rk1AdaptStepSetupConstants(ai,bi,ci,cip)
    xn = a
    yn = y0
    h = 0.001d0
    do while(xn.le.b)
      if(xn.eq.b) exit
      goodStep = .false.
      do while(goodStep .eqv. .false.)
        if(xn+h>b) then
          h = b-xn
        end if
        k1 = h*f(nF,nC,consts,xn,yn)
        k2 = h*f(nF,nC,consts,xn+ai(2)*h,yn+bi(2,1)*k1)
        k3 = h*f(nF,nC,consts,xn+ai(3)*h,yn+bi(3,1)*k1+bi(3,2)*k2)
        k4 = h*f(nF,nC,consts,xn+ai(4)*h,yn+bi(4,1)*k1+bi(4,2)*k2+bi(4,3)*k3)
        k5 = h*f(nF,nC,consts,xn+ai(5)*h,yn+bi(5,1)*k1+bi(5,2)*k2+bi(5,3)*k3+bi(5,4)*k4)
        k6 = h*f(nF,nC,consts,xn+ai(6)*h,yn+bi(6,1)*k1+bi(6,2)*k2+bi(6,3)*k3+bi(6,4)*k4+bi(6,5)*k5)
        kci = k1*ci(1)+k2*ci(2)+k3*ci(3)+k4*ci(4)+k5*ci(5)+k6*ci(6)
        kcip = k1*cip(1)+k2*cip(2)+k3*cip(3)+k4*cip(4)+k5*cip(5)+k6*cip(6)
        ynT = yn+kci
        ynsT = yn+kcip
        dyn = abs(ynT(1)-ynsT(1))
        if(dyn .gt. 1.d-10) then
          h = h*((1.d-12/dyn)**(0.2d0))
        else
          goodStep = .true.
        end if
      end do
      xn = xn+h
      yns = yn+kcip
      yn = yn+kci
      dyn = abs(yn(1)-yns(1))
      h = h*((1.d-12/dyn)**(0.2d0))
    end do
    rk1AdaptStepCmplxC = yn
  end function

  subroutine rk1AdaptStepSetupConstants(ai,bi,ci,cip)
    !! Assigns values of constants for RK adaptive routine

    real(dp) :: ai(6), bi(6,6), ci(6), cip(6)
    ai(1:6) = 0.d0
    bi(1:6,1:6) = 0.d0
    ci(1:6) = 0.d0
    cip(1:6) = 0.d0
    ai(1) = 0.d0
    ai(2) = 1.d0/5.d0
    ai(3) = 3.d0/10.d0
    ai(4) = 3.d0/5.d0
    ai(5) = 1.d0
    ai(6) = 7.d0/8.d0
    bi(2,1) = 1.d0/5.d0
    bi(3,1) = 3.d0/40.d0
    bi(3,2) = 9.d0/40.d0
    bi(4,1) = 3.d0/10.d0
    bi(4,2) = -9.d0/10.d0
    bi(4,3) = 6.d0/5.d0
    bi(5,1) = -11.d0/54.d0
    bi(5,2) = 5.d0/2.d0
    bi(5,3) = -70.d0/27.d0
    bi(5,4) = 35.d0/27.d0
    bi(6,1) = 1631.d0/55296.d0
    bi(6,2) = 175.d0/512.d0
    bi(6,3) = 575.d0/13824.d0
    bi(6,4) = 44275.d0/110592.d0
    bi(6,5) = 253.d0/4096.d0
    ci(1) = 37.d0/378.d0
    ci(2) = 0.d0
    ci(3) = 250.d0/621.d0
    ci(4) = 125.d0/594.d0
    ci(5) = 0.d0
    ci(6) = 512.d0/1771.d0
    cip(1) = 2825.d0/27648.d0
    cip(2) = 0.d0
    cip(3) = 18575.d0/48384.d0
    cip(4) = 13525.d0/55296.d0
    cip(5) = 277.d0/14336.d0
    cip(6) = 1.d0/4.d0
  end subroutine rk1AdaptStepSetupConstants

end module odeSolver
