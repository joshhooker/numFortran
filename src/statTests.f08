!*********************************************************
!>
! Collection of statisitcal calculations:
!
! * mean - Calculates mean of an array
! * stdDev - Calculates Mean and Std. Deviation

module statTests
  implicit none

  private
  integer, parameter :: dp = kind(0.d0)

  public :: mean
  public :: stdDev, stdDevInt

  interface mean
    module procedure mean_i, mean_r, mean_d
  end interface

  ! interface sDev
  !   module procedure sDev_

contains

  real(dp) function meanFunc(array)
    integer :: i, n
    real(dp) :: array(:), arraySum
    n = size(array,1)
    arraySum = 0.d0
    do i=1,n
      arraySum = arraySum+array(i)
    end do
    meanFunc = arraySum/dble(n)
  end function

  real(dp) function mean_i(array)
    integer :: i, n, array(:)
    real(dp), allocatable :: newArray(:)
    n = size(array,1)
    allocate(newArray(n))
    do i=1,n
      newArray(i) = dble(array(i))
    end do
    mean_i = meanFunc(newArray)
    deallocate(newArray)
  end function

  real(dp) function mean_r(array)
    integer :: i, n
    real :: array(:)
    real(dp), allocatable :: newArray(:)
    n = size(array,1)
    allocate(newArray(n))
    do i=1,n
      newArray(i) = dble(array(i))
    end do
    mean_r = meanFunc(newArray)
    deallocate(newArray)
  end function

  real(dp) function mean_d(array)
    real(dp) :: array(:)
    mean_d = meanFunc(array)
  end function


  subroutine stdDev(array, mean, sigma)
    !! date: January 9, 2016
    !! version: v0.1
    !!
    !! Calculates the mean and standard deviation of
    !! an input array of doubles

    integer :: i, n
    real(dp) :: arraySum, arrayStdDev
    real(dp) :: array(:)
      !! input: array of numbers of size n
    real(dp) :: mean
      !! output: mean of array
    real(dp) :: sigma
      !! output: standard deviation of array

    n = size(array,1)

    arraySum = 0.d0
    do i=1,n
      arraySum = arraySum + array(i)
    end do
    mean = arraySum/dble(n)

    arrayStdDev = 0.d0
    do i=1,n
      arrayStdDev = arrayStdDev+(array(i)-mean)*(array(i)-mean)
    end do
    arrayStdDev = arrayStdDev/dble(n)
    sigma = sqrt(arrayStdDev)

  end subroutine stdDev

  subroutine stdDevInt(array, mean, sigma)
    !! date: January 10, 2016
    !! version: v0.1
    !!
    !! Calculates the mean and standard deviation of
    !! an input array of integers

    integer :: i, n
    integer :: array(:)
      !! input: array of numbers of size n
    real(dp) :: arraySum, arrayStdDev
    real(dp) :: mean
      !! output: mean of array
    real(dp) :: sigma
      !! output: standard deviation of array

    n = size(array,1)

    arraySum = 0.d0
    do i=1,n
      arraySum = arraySum + dble(array(i))
    end do
    mean = arraySum/dble(n)

    arrayStdDev = 0.d0
    do i=1,n
      arrayStdDev = arrayStdDev+(dble(array(i))-mean)*(dble(array(i))-mean)
    end do
    arrayStdDev = arrayStdDev/dble(n)
    sigma = sqrt(arrayStdDev)

  end subroutine stdDevInt



end module statTests