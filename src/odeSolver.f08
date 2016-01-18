!*********************************************************
!>
! Functions to solve ODES:
!
!  * 1st order solvers:
!    * rk1FixedStep
!    * rk1AdaptStep

module odeSolver
  use parameters
  implicit none

  private :: rk1AdaptStepSetupConstants
contains

  real(dp) function odeFunction(x,y)
    !! Function in the form \[y' = f(x,y)\] to be solved
    real(dp) :: x
      !! input: x value in f(x,y)
    real(dp) :: y
      !! input: y value in f(x,y)
    odeFunction = y*y+1.d0
  end function odeFunction

  real(dp) function rk1FixedStep(a,b,h,y0)
    !! date: January 8, 2016
    !! version: v0.1
    !!
    !! 4th order Runga Kutta ODE Solver with fixed step h and initial
    !! value \(y(a) = y_0\). If the fixed step does not evenly divide
    !! into \(b-a\), a smaller step is placed at the end so the final
    !! \(x\) position is at \(b\)

    integer :: i
    real(dp) :: xn, yn, k1, k2, k3, k4
    real(dp) :: a
      !! input: starting x value
    real(dp) :: b
      !! input: final x value
    real(dp) :: h
      !! input: step size
    real(dp) :: y0
      !! input: inital value \(y(a) = y0\)

    xn = a
    yn = y0
    do while (xn.lt.b)
      if(xn+h.gt.b) then
        h = b-xn
      end if
      k1 = h*odeFunction(xn,yn)
      k2 = h*odeFunction(xn+h/2.d0,yn+k1/2.d0)
      k3 = h*odeFunction(xn+h/2.d0,yn+k2/2.d0)
      k4 = h*odeFunction(xn+h,yn+k3)
      xn = xn+h
      yn = yn+(1.d0/6.d0)*(k1+2.d0*k2+2.d0*k3+k4)
    end do
    rk1FixedStep = yn
  end function rk1FixedStep

  real(dp) function rk1AdaptStep(a,b,y0)
    !! date: January 8, 2016
    !! version: v0.1
    !!
    !! Runga Kutta ODE Solver with adaptive h size to satisfy a fixed
    !! error and initial value \(y(a) = y_0\). Calculates \(y\) for two
    !! different \(c_i\)'s in \(y = y + \sum_i c_i*k_i\) and compares
    !! them to determine if the step size is too big or small.

    integer :: i,j
    real(dp) :: a
      !! input: starting x value
    real(dp) :: b
      !! input: final x value
    real(dp) :: y0
      !! input: inital value \(y(a) = y0\)
    real(dp) :: xn,yn,yns,dyn
    real(dp) :: xni,yni,kci,kcip,k(6)
    real(dp) :: ai(6),bi(6,6),ci(6),cip(6)
    real(dp) :: h
    call rk1AdaptStepSetupConstants(ai,bi,ci,cip)
    xn = a
    yn = y0
    h = 0.1d0
    do while(xn<b)
      if(xn+h>b) then
        h = b-xn
      end if
      kci = 0.d0
      kcip = 0.d0
      do i=1,6
        xni = xn+ai(i)*h
        yni = yn
        if (i.gt.1) then
          do j=1,i-1
            yni = yni+bi(i,j)*k(j)
          end do
        end if
        k(i) = h*odeFunction(xni,yni)
        kci = kci+k(i)*ci(i)
        kcip = kcip+k(i)*cip(i)
      end do
      xn = xn+h
      yns = yn+kcip
      yn = yn+kci
      dyn = abs(yn-yns)
      h = h*((1.d-12/dyn)**(0.2d0))
    end do
    rk1AdaptStep = yn
  end function rk1AdaptStep

  subroutine rk1AdaptStepSetupConstants(ai,bi,ci,cip)
    !! date: January 8, 2016
    !! version: v0.1
    !!
    !! Assigns values of constants for RK adaptive routine

    real(dp) :: ai(6)
      !! output: \(a_i\) array, of size 6, used for steps in \(x\)
    real(dp) :: bi(6,6)
      !! output: \(b_i\) array, of size 6x6, used for steps in \(y\)
    real(dp) :: ci(6)
      !! output: \(c_i\) array, of size 6, used to add \(c_i k\) for \(y\)
    real(dp) :: cip(6)
      !! output: \(c_{ip}\) array, of size 6, used to add \(c_{ip} k\) for \(y_p\)


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
