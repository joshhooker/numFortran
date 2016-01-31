!*********************************************************
!>
! Integration Functions:
!
!  * Gauss-Legendre Integration:

module integration
  implicit none

  private
  integer, parameter :: dp = kind(0.d0)

  real(dp) :: glx128(64), glw128(64)
  real(dp) :: glx512(256), glw512(256)
  real(dp) :: glx1024(512), glw1024(512)

  public :: readGLParameters
  public :: gaussLegendre, gaussLegendreCmplx

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
  end subroutine

  real(dp) function gaussLegendre(f,a,b,n,consts)
    integer :: i, n
    real(dp) :: consts(n), a, b, t0, dt, sumIntegral
    interface
      function f(n,c,x)
        integer, parameter :: dp = kind(0.d0)
        integer :: n
        real(dp) :: c(n), x, f
      end function f
    end interface

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
    !    sumIntegral = sumIntegral+glw512(i)*f(n,consts,t0+dt*glx512(i))
    ! end do

    ! Using glx1024 and glw1024
    ! do i=512,1,-1
    !   sumIntegral = sumIntegral+glw1024(i)*f(n,consts,t0-dt*glx1024(i))
    ! end do
    ! do i=1,512
    !   sumIntegral = sumIntegral+glw1024(i)*f(n,consts,t0+dt*glx1024(i))
    ! end do
    gaussLegendre = sumIntegral*dt
  end function

  complex(dp) function gaussLegendreCmplx(f,a,b,n,consts)
    integer :: i, n
    real(dp) :: a, b, t0, dt, sumReal, sumCmplx
    complex(dp) :: consts(n)
    interface
      function f(n,c,x)
        integer, parameter :: dp = kind(0.d0)
        integer :: n
        real(dp) :: x
        complex(dp) :: c(n), f
      end function
    end interface

    t0 = (a+b)/2.d0
    dt = (b-a)/2.d0
    sumReal = 0.d0
    sumCmplx = 0.d0

    ! Using glx128 and glw128
    do i=64,1,-1
      sumReal = sumReal+glw128(i)*real(f(n,consts,t0-dt*glx128(i)))
      sumCmplx = sumCmplx+glw128(i)*aimag(f(n,consts,t0-dt*glx128(i)))
    end do
    do i=1,64
      sumReal = sumReal+glw128(i)*real(f(n,consts,t0+dt*glx128(i)))
      sumCmplx = sumCmplx+glw128(i)*aimag(f(n,consts,t0+dt*glx128(i)))
    end do

    ! Using glx512 and glw512
    ! do i=256,1,-1
    !   sumIntegral = sumIntegral+glw512(i)*f(n,consts,t0-dt*glx512(i))
    ! end do
    ! do i=1,256
    !    sumIntegral = sumIntegral+glw512(i)*f(n,consts,t0+dt*glx512(i))
    ! end do

    ! Using glx1024 and glw1024
    ! do i=512,1,-1
    !   sumIntegral = sumIntegral+glw1024(i)*f(n,consts,t0-dt*glx1024(i))
    ! end do
    ! do i=1,512
    !   sumIntegral = sumIntegral+glw1024(i)*f(n,consts,t0+dt*glx1024(i))
    ! end do

    gaussLegendreCmplx = cmplx(sumReal*dt, sumCmplx*dt)
  end function

end module integration
