!*********************************************************
!>
! Collection of statisitcal calculations:
!

module statTests
  use parameters
  implicit none

contains

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