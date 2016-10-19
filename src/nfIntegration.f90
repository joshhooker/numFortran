!*********************************************************
!>
! Integration Functions:
!
!  * Gauss-Legendre Integration:

module nfIntegration
  use nfConstants
  use nfSpecialFunctions
  implicit none

  real(dp) :: glx128(64), glw128(64)
  real(dp) :: glx512(256), glw512(256)
  real(dp) :: glx1024(512), glw1024(512)
  real(dp) :: g30k61_gx(15), g30k61_gw(15), g30k61_kx(31), g30k61_kw(31)

  type, public :: nfInt
    real(dp), public, allocatable :: glx(:), glw(:)
  contains
    procedure, public :: setGLNumPoints
    procedure, private :: setGLPoints
    procedure, private :: bisectionGLPoints
  end type nfInt

  type, private :: nfIntGLPoints
    real(dp), private, allocatable :: lower(:), upper(:)
  end type nfIntGLPoints

  public :: readGLParameters
  public :: gaussLegendre, gaussLegendreCmplx
  public :: readGKParameters

contains

  subroutine setGLNumPoints(this, numPoints)
    class(nfInt) :: this
    integer :: numPoints

    allocate(this%glx(numPoints))
    allocate(this%glw(numPoints))
    call setGLPoints(this)
  end subroutine setGLNumPoints

  subroutine setGLPoints(this)
    class(nfInt) :: this
    type(nfIntGLPoints) :: points
    integer :: i, numPoints
    real(dp) :: x, dx, value, valueOLD

    numPoints = size(this%glx,1)
    allocate(points%lower(numPoints))
    allocate(points%upper(numPoints))
    dx = 0.000001
    i = 1
    valueOLD = legendrePoly(numPoints, -1)
    do x=-1.d0+dx,1.d0,dx
      value = legendrePoly(numPoints, x)
      if(value*valueOLD.lt.0.d0) then
        points%lower(i) = x-dx
        points%upper(i) = x
        i = i+1
        print *, points%lower(i), points%upper(i)
      end if
      valueOLD = value
    end do
    ! do i=1,numPoints
    !   print *, points%lower(i), points%upper(i)
    ! end do
    !print *, 'hi 2'
    call bisectionGLPoints(this, points)
  end subroutine setGLPoints

  subroutine bisectionGLPoints(this, points)
    class(nfInt) :: this
    type(nfIntGLPoints) :: points
    integer :: i, numPoints
    real(dp) :: x, value, eps

    eps = 1.d-10
    numPoints = size(this%glx,1)

    do i=1,numPoints
      do while(abs(points%upper(i)-points%lower(i)).gt.eps)
        x = (points%upper(i)+points%lower(i))/2.d0
        value = legendrePoly(numPoints, x)
        !print *, points%lower(i), points%upper(i), x, value
        if(value.lt.0.d0) then
          points%upper(i) = x
        else
          points%lower(i) = x
        end if
      end do
      !print *, x
    end do



  end subroutine bisectionGLPoints

  subroutine readGLParameters()
    implicit none
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

  subroutine readGKParameters()
    implicit none
    integer :: i, numLines

    open(unit=1,file='data/g30k61.dat',status='old')
    read(1,*) numLines
    do i=1,numLines
      read(1,*) g30k61_gx(i), g30k61_gw(i)
    end do
    read(1,*) numLines
    do i=1,numLines
      read(1,*) g30k61_kx(i), g30k61_kw(i)
    end do
    close(1)
  end subroutine

  real(dp) function gaussLegendre(f,a,b,n,consts)
    implicit none
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
    implicit none
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
    gaussLegendreCmplx = cmplx(sumReal*dt,sumCmplx*dt,dp)
  end function

end module nfIntegration