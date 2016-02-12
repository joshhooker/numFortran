!*********************************************************
!>
! Collection of statisitcal calculations:
!
! * meanArr - Calculates mean of an array
! * stdDev - Calculates Mean and Std. Deviation
!
! Future:
! * bootstrap, etc.

module statTests
  use randomNumbers
  use sorting
  implicit none

  private
  integer, parameter :: dp = kind(0.d0)

  public :: medianArr
  public :: meanArr
  public :: stdDev

  interface medianArr
    module procedure median_i, median_r, median_d
  end interface

  interface meanArr
    module procedure mean_i, mean_r, mean_d
  end interface

  interface stdDev
    module procedure stdDev_i, stdDev_r, stdDev_d
  end interface

contains

  real(dp) function median_i(array)
    implicit none
    integer :: n, array(:)

    n = size(array,1)
    if(mod(n,2).eq.0) then
      median_i = array(n/2)
    else
      median_i = (array((n+1)/2-1)+array((n+1)/2+1))/2.d0
    end if
  end function

  real(dp) function median_r(array)
    implicit none
    integer :: n
    real :: array(:)

    n = size(array,1)
    if(mod(n,2).eq.0) then
      median_r = array(n/2)
    else
      median_r = (array((n+1)/2-1)+array((n+1)/2+1))/2.d0
    end if
  end function

  real(dp) function median_d(array)
    implicit none
    integer :: n
    real(dp) :: array(:)

    n = size(array,1)
    if(mod(n,2).eq.0) then
      median_d = array(n/2)
    else
      median_d = (array((n+1)/2-1)+array((n+1)/2+1))/2.d0
    end if
  end function

  real(dp) function mean_i(array)
    implicit none
    integer :: i, n, array(:)
    real(dp) :: arraySum

    n = size(array,1)
    arraySum = 0.d0
    do i=1,n
      arraySum = arraySum+dble(array(i))
    end do
    mean_i = arraySum/dble(n)
  end function

  real(dp) function mean_r(array)
    implicit none
    integer :: i, n
    real :: array(:)
    real(dp) :: arraySum

    n = size(array,1)
    arraySum = 0.d0
    do i=1,n
      arraySum = arraySum+dble(array(i))
    end do
    mean_r = arraySum/dble(n)
  end function

  real(dp) function mean_d(array)
    implicit none
    integer :: i, n
    real(dp) :: array(:), arraySum

    n = size(array,1)
    arraySum = 0.d0
    do i=1,n
      arraySum = arraySum+array(i)
    end do
    mean_d = arraySum/dble(n)
  end function

  subroutine stdDev_i(array, mean, sigma)
    implicit none
    integer :: i, n, array(:)
    real(dp) :: arrayStdDev, mean, sigma

    n = size(array,1)
    mean = meanArr(array)
    arrayStdDev = 0.d0
    do i=1,n
      arrayStdDev = arrayStdDev+(dble(array(i))-mean)* &
          (dble(array(i))-mean)
    end do
    arrayStdDev = arrayStdDev/dble(n)
    sigma = sqrt(arrayStdDev)
  end subroutine

  subroutine stdDev_r(array, mean, sigma)
    implicit none
    integer :: i, n
    real :: array(:)
    real(dp) :: arrayStdDev, mean, sigma

    n = size(array,1)
    mean = meanArr(array)
    arrayStdDev = 0.d0
    do i=1,n
      arrayStdDev = arrayStdDev+(dble(array(i))-mean)* &
          (dble(array(i))-mean)
    end do
    arrayStdDev = arrayStdDev/dble(n)
    sigma = sqrt(arrayStdDev)
  end subroutine

  subroutine stdDev_d(array, mean, sigma)
    implicit none
    integer :: i, n
    real(dp) :: array(:), arrayStdDev, mean, sigma

    n = size(array,1)
    mean = meanArr(array)
    arrayStdDev = 0.d0
    do i=1,n
      arrayStdDev = arrayStdDev+(array(i)-mean)* &
          (array(i)-mean)
    end do
    arrayStdDev = arrayStdDev/dble(n)
    sigma = sqrt(arrayStdDev)
  end subroutine

  subroutine bootstrap(array, mean, sigma, nboot)
    implicit none
    integer :: i, n, nboot, boot_i, clock
    real(dp) :: array(:), mean, sigma, bootMean(nboot)
    real(dp), allocatable :: bootArr(:)

    call system_clock(count=clock)
    call randSeedGenerator(clock)
    n = size(array,1)
    allocate(bootArr(n))
    do boot_i=1,nboot
      do i=1,n
        bootArr(i) = array(xorshift128Int(n)+1)
      end do
      call qsort(bootArr)
    end do
    deallocate(bootArr)

  end subroutine

end module