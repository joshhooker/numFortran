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
  implicit none

  private
  integer, parameter :: dp = kind(0.d0)

  public :: meanArr, stdDev

  interface meanArr
    module procedure mean_i, mean_r, mean_d
  end interface

  interface stdDev
    module procedure stdDev_i, stdDev_r, stdDev_d
  end interface

contains

  real(dp) function mean_i(array)
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



end module statTests