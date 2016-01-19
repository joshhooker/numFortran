!*********************************************************
!>
! Integration Functions:
!
!  * Gauss-Legendre Integration:

module integration
  implicit none

  private :: dp
  integer, parameter :: dp = kind(0.d0)

  real(dp) :: glx128(64), glw128(64)
  real(dp) :: glx512(256), glw512(256)
  real(dp) :: glx1024(512), glw1024(512)

  private :: glx128, glw128
  private :: glx512, glw512
  private :: glx1024, glw1024

contains

  subroutine readGLParameters()
    integer :: i

    open(unit=1,file='data/gl128.dat',status='old')
    do i=1,64
      read(1,*) glx128(i), glw128(i)
    end do
    close(1)

    open(unit=1,file='data/gl512.dat',status='old')
    do i=1,256
      read(1,*) glx512(i), glw512(i)
    end do
    close(1)

    open(unit=1,file='data/gl1024.dat',status='old')
    do i=1,512
      read(1,*) glx1024(i), glw1024(i)
    end do
    close(1)

  end subroutine readGLParameters

  real(dp) function gaussLegendre(f,a,b,n,consts)
    !! date: January 14, 2016
    !! version: v0.1
    !!
    !! Integrates some function f from a to b
    integer :: i
    integer :: n
      !! input: number of constants
    real(dp) :: consts(n)
      !! input: array of constants of size n
    real(dp) :: a
      !! input: lower bound of integral
    real(dp) :: b
      !! input: upper bound of integral
    interface
      function f(n,c,x)
        integer, parameter :: dp = kind(0.d0)
        integer :: n
        real(dp) :: c(n)
        real(dp) :: x
        real(dp) :: f
      end function f
    end interface
    real(dp) :: t0, dt
    real(dp) :: sumIntegral

    t0 = (a+b)/2.d0
    dt = (b-a)/2.d0
    sumIntegral = 0.d0

    ! Using glx128 and glw128
    do i=64,1,-1
      sumIntegral = sumIntegral+glw128(i)*f(n,consts,t0-dt*glx128(i))
    end do
    do i=1,64
      sumIntegral = sumIntegral+glw128(i)*f(n,consts,t0+dt*glx128(i))
    end do

    ! Using glx512 and glw512
   ! do i=256,1,-1
   !   sumIntegral = sumIntegral+glw512(i)*f(n,consts,t0-dt*glx512(i))
   ! end do
   ! do i=1,256
   !   sumIntegral = sumIntegral+glw512(i)*f(n,consts,t0+dt*glx512(i))
   ! end do

    ! Using glx1024 and glw1024
   ! do i=512,1,-1
   !   sumIntegral = sumIntegral+glw1024(i)*f(n,consts,t0-dt*glx1024(i))
   ! end do
   ! do i=1,512
   !   sumIntegral = sumIntegral+glw1024(i)*f(n,consts,t0+dt*glx1024(i))
   ! end do

    gaussLegendre = sumIntegral*dt
  end function gaussLegendre

end module integration
