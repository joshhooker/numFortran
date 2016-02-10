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

  real(dp), dimension(6), parameter :: rk1ai = [0.d0, 1.d0/5.d0, 3.d0/1.d1, 3.d0/5.d0, 1.d0, 7.d0/8.d0]
  real(dp), dimension(6,6), parameter :: rk1bi = reshape((/0.d0, 1.d0/5.d0, 3.d0/4.d1, 3.d0/1.d1, -11.d0/54.d0, 1631.d0/55296.d0, &
      0.d0, 0.d0, 9.d0/4.d1, -9.d0/1.d1, 5.d0/2.d0, 175.d0/512.d0, 0.d0, 0.d0, 0.d0, 6.d0/5.d0, -7.d1/27.d0, 575.d0/13824.d0, &
      0.d0, 0.d0, 0.d0, 0.d0, 35.d0/27.d0, 44275.d0/110592.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 253.d0/4096.d0, 0.d0, 0.d0, 0.d0, &
      0.d0, 0.d0, 0.d0, 0.d0/),(/6,6/))
  real(dp), dimension(6), parameter :: rk1ci = [37.d0/378.d0, 0.d0, 25.d1/621.d0, 125.d0/594.d0, 0.d0, 512.d0/1771.d0]
  real(dp), dimension(6), parameter :: rk1cip = [2825.d0/27648.d0, 0.d0, 18575.d0/48384.d0, 13525.d0/55296.d0, &
      277.d0/14336.d0, 1.d0/4.d0]

  public :: rk1FixedDriver, rk1FixedStep
  public :: rk1AdaptDriver, rk1AdaptStep
  public :: rk1AdaptStepCmplxC

contains

  subroutine rk1FixedDriver(f,nF,x,y,h,nC,consts)
    integer :: nF, nC
    real(dp) :: x, y(nF), k1(nF), k2(nF), k3(nF), k4(nF), h, consts(nC)
    interface
      function f(nF,nC,c,x,y)
        integer, parameter :: dp = kind(0.d0)
        integer :: nF, nC
        real(dp) :: c(nC), x, y(nF), f(nF)
      end function f
    end interface
    k1 = h*f(nF,nC,consts,x,y)
    k2 = h*f(nF,nC,consts,x+h/2.d0,y+k1/2.d0)
    k3 = h*f(nF,nC,consts,x+h/2.d0,y+k2/2.d0)
    k4 = h*f(nF,nC,consts,x+h,y+k3)
    y = y+(1.d0/6.d0)*(k1+2.d0*k2+2.d0*k3+k4)
  end subroutine

  function rk1FixedStep(f,nF,a,b,h,y0,nC,consts)
    !! 4th order Runga Kutta ODE Solver with fixed step h
    integer :: i, nF, nC
    real(dp) :: xn, a, b, h
    real(dp) :: y0(nF), yi(nF), yn(nF), consts(nC), rk1FixedStep(nF)
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
      call rk1FixedDriver(f,nF,xn,yn,h,nC,consts)
      xn = xn+h
    end do
    rk1FixedStep = yn
  end function

  subroutine rk1AdaptDriver(f,nF,x,y,h,nC,consts,err)
    integer nF, nC
    real(dp) :: x, y(nF), ys(nF), ynT(nF), ynsT(nF), dyn, consts(nC), h, err
    real(dp) :: k1(nF), k2(nF), k3(nF), k4(nF), k5(nF), k6(nF), kci(nF), kcip(nF)
    logical :: goodStep
    interface
      function f(nF,nC,c,x,y)
        integer, parameter :: dp = kind(0.d0)
        integer :: nF, nC
        real(dp) :: c(nC), x, y(nF), f(nF)
      end function f
    end interface
    goodStep = .false.
    do while(goodStep .eqv. .false.)
      k1 = h*f(nF,nC,consts,x,y)
      k2 = h*f(nF,nC,consts,x+rk1ai(2)*h,y+rk1bi(2,1)*k1)
      k3 = h*f(nF,nC,consts,x+rk1ai(3)*h,y+rk1bi(3,1)*k1+rk1bi(3,2)*k2)
      k4 = h*f(nF,nC,consts,x+rk1ai(4)*h,y+rk1bi(4,1)*k1+rk1bi(4,2)*k2+rk1bi(4,3)*k3)
      k5 = h*f(nF,nC,consts,x+rk1ai(5)*h,y+rk1bi(5,1)*k1+rk1bi(5,2)*k2+rk1bi(5,3)*k3+rk1bi(5,4)*k4)
      k6 = h*f(nF,nC,consts,x+rk1ai(6)*h,y+rk1bi(6,1)*k1+rk1bi(6,2)*k2+rk1bi(6,3)*k3+rk1bi(6,4)*k4+rk1bi(6,5)*k5)
      kci = k1*rk1ci(1)+k2*rk1ci(2)+k3*rk1ci(3)+k4*rk1ci(4)+k5*rk1ci(5)+k6*rk1ci(6)
      kcip = k1*rk1cip(1)+k2*rk1cip(2)+k3*rk1cip(3)+k4*rk1cip(4)+k5*rk1cip(5)+k6*rk1cip(6)
      ynT = y+kci
      ynsT = y+kcip
      dyn = abs(ynT(1)-ynsT(1))
      if(dyn .gt. 1.d-8) then
        h = h*((err/dyn)**(0.2d0))
      else
        x = x+h
        y = ynT
        h = h*((err/dyn)**(0.2d0))
        goodStep = .true.
      end if
    end do
  end subroutine

  function rk1AdaptStep(f,nF,a,b,y0,nC,consts)
    !! Runga Kutta ODE Solver with adaptive h size
    integer :: i, j, nF, nC
    real(dp) :: xn, a, b, h, yn(nF), y0(nF), consts(nC), rk1AdaptStep(nF)
    interface
      function f(nF,nC,c,x,y)
        integer, parameter :: dp = kind(0.d0)
        integer :: nF, nC
        real(dp) :: c(nC), x, y(nF), f(nF)
      end function f
    end interface
    xn = a
    yn = y0
    h = 0.001d0
    call rk1AdaptDriver(f,nF,xn,yn,h,nC,consts,1.d-12)
    do while(xn<b)
      if(xn+h>b) then
        h = b-xn
      end if
      call rk1AdaptDriver(f,nF,xn,yn,h,nC,consts,1.d-12)
    end do
    rk1AdaptStep = yn
  end function

  function rk1AdaptStepCmplxC(f,nF,a,b,y0,nC,consts)
    !! Runga Kutta ODE Solver with adaptive h size with Complex constant arguments
    integer :: i, j, nF, nC
    real(dp) :: xn, k1(nF), k2(nF), k3(nF), k4(nF), k5(nF), k6(nF)
    real(dp) :: a, b, h, y0(nF), yi(nF), yn(nF), yns(nF), dyn, ynT(nF), ynsT(nF)
    real(dp) :: xni, yni(nF), kci(nF), kcip(nF), rk1AdaptStepCmplxC(nF)
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
        k2 = h*f(nF,nC,consts,xn+rk1ai(2)*h,yn+rk1bi(2,1)*k1)
        k3 = h*f(nF,nC,consts,xn+rk1ai(3)*h,yn+rk1bi(3,1)*k1+rk1bi(3,2)*k2)
        k4 = h*f(nF,nC,consts,xn+rk1ai(4)*h,yn+rk1bi(4,1)*k1+rk1bi(4,2)*k2+rk1bi(4,3)*k3)
        k5 = h*f(nF,nC,consts,xn+rk1ai(5)*h,yn+rk1bi(5,1)*k1+rk1bi(5,2)*k2+rk1bi(5,3)*k3+rk1bi(5,4)*k4)
        k6 = h*f(nF,nC,consts,xn+rk1ai(6)*h,yn+rk1bi(6,1)*k1+rk1bi(6,2)*k2+rk1bi(6,3)*k3+rk1bi(6,4)*k4+rk1bi(6,5)*k5)
        kci = k1*rk1ci(1)+k2*rk1ci(2)+k3*rk1ci(3)+k4*rk1ci(4)+k5*rk1ci(5)+k6*rk1ci(6)
        kcip = k1*rk1cip(1)+k2*rk1cip(2)+k3*rk1cip(3)+k4*rk1cip(4)+k5*rk1cip(5)+k6*rk1cip(6)
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

end module
