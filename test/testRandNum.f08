module testRandNum
  use randomNumbers
  use statTests

  use iso_fortran_env
  implicit none

  private
  integer, parameter :: dp = kind(0.d0)

  public :: randNumTest

contains

  subroutine randNumTest()
    implicit none

    integer :: i

    ! Variables to time functions
    integer(int64) :: t1, t2, clock_rate, clock_max

    ! Variables to test random number generators
    integer, parameter :: randArrayN = 50000000
    integer, parameter :: randChiSqK = 500
    integer :: randArrayInt(randArrayN)
    real(dp) :: randArray(randArrayN)
    real(dp) :: randMean, randStdDev

    integer, parameter :: chiSqDOF = 999
    real(dp) :: chiSq

    write(*,'(a)') '*************************************'
    write(*,'(a,1x,i0,1x,a)') 'Testing random number generators with', randArrayN,'numbers:'

    call randSeedGenerator(1234)

    !xorshift
    write(*,'(2x,a)') 'xorshift:'
    call system_clock(t1, clock_rate, clock_max)
    do i=1,randArrayN
      randArray(i) = xorshift()
    end do
    call system_clock(t2, clock_rate, clock_max)
    write(*,'(4x,a,1x,i0,1x,a,1x,es12.5,1x,a)') 'Time to generate', randArrayN, 'numbers =', dble(t2-t1)/dble(clock_rate), ' sec'
    call stdDev(randArray,randMean,randStdDev)
    write(*,'(3x,4(1x,a,1x,es12.5))') 'mean =', randMean, 'stdDev =', randStdDev, &
              'min =', minval(randArray,1), 'max =', maxval(randArray,1)
    call chiSqTest(randArray, chiSqDOF+1, 0.d0, 1.d0, chiSq)
    write(*,'(4x,a,1x,es12.5,1x,a,1x,i0,1x,a)') 'chi-squared value:', chiSq, 'with', chiSqDOF, 'degrees of freedom'
    write(*,*)

    write(*,'(2x,a)') 'xorshiftInt:'
    call randSeedGenerator(12345679)
    call system_clock(t1, clock_rate, clock_max)
    do i=1,randArrayN
      randArrayInt(i) = xorshiftInt(101)
    end do
    call system_clock(t2, clock_rate, clock_max)
    write(*,'(4x,a,1x,i0,1x,a,1x,es12.5,1x,a)') 'Time for', randArrayN, ' numbers =', dble(t2-t1)/dble(clock_rate), 'sec'
    call stdDev(randArrayInt,randMean,randStdDev)
    write(*,'(3x,2(1x,a,1x,es12.5),2(1x,a,1x,i0))') 'mean =', randMean, 'stdDev =', randStdDev, &
              'min =', minval(randArrayInt,1), 'max =', maxval(randArrayInt,1)
    write(*,*)

    write(*,'(2x,a)') 'xorshift64star:'
    call randSeedGenerator(12345679)
    call system_clock(t1, clock_rate, clock_max)
    do i=1,randArrayN
      randArray(i) = xorshift64star()
    end do
    call system_clock(t2, clock_rate, clock_max)
    write(*,'(4x,a,1x,i0,1x,a,1x,es12.5,1x,a)') 'Time to generate', randArrayN, 'numbers =', dble(t2-t1)/dble(clock_rate), ' sec'
    call stdDev(randArray,randMean,randStdDev)
    write(*,'(3x,4(1x,a,1x,es12.5))') 'mean =', randMean, 'stdDev =', randStdDev, &
              'min =', minval(randArray,1), 'max =', maxval(randArray,1)
    call chiSqTest(randArray, chiSqDOF+1, 0.d0, 1.d0, chiSq)
    write(*,'(4x,a,1x,es12.5,1x,a,1x,i0,1x,a)') 'chi-squared value:', chiSq, 'with', chiSqDOF, 'degrees of freedom'
    write(*,*)

    write(*,'(2x,a)') 'xorshift64starInt:'
    call randSeedGenerator(12345679)
    call system_clock(t1, clock_rate, clock_max)
    do i=1,randArrayN
      randArrayInt(i) = xorshift64starInt(101)
    end do
    call system_clock(t2, clock_rate, clock_max)
    write(*,'(4x,a,1x,i0,1x,a,1x,es12.5,1x,a)') 'Time for', randArrayN, ' numbers =', dble(t2-t1)/dble(clock_rate), 'sec'
    call stdDev(randArrayInt,randMean,randStdDev)
    write(*,'(3x,2(1x,a,1x,es12.5),2(1x,a,1x,i0))') 'mean =', randMean, 'stdDev =', randStdDev, &
              'min =', minval(randArrayInt,1), 'max =', maxval(randArrayInt,1)
    write(*,*)

    write(*,'(2x,a)') 'xorshift128:'
    call system_clock(t1, clock_rate, clock_max)
    do i=1,randArrayN
      randArray(i) = xorshift128()
    end do
    call system_clock(t2, clock_rate, clock_max)
    write(*,'(4x,a,1x,i0,1x,a,1x,es12.5,1x,a)') 'Time to generate', randArrayN, 'numbers =', dble(t2-t1)/dble(clock_rate), ' sec'
    call stdDev(randArray,randMean,randStdDev)
    write(*,'(3x,4(1x,a,1x,es12.5))') 'mean =', randMean, 'stdDev =', randStdDev, &
              'min =', minval(randArray,1), 'max =', maxval(randArray,1)
    call chiSqTest(randArray, chiSqDOF+1, 0.d0, 1.d0, chiSq)
    write(*,'(4x,a,1x,es12.5,1x,a,1x,i0,1x,a)') 'chi-squared value:', chiSq, 'with', chiSqDOF, 'degrees of freedom'
    write(*,*)

    write(*,'(2x,a)') 'xorshift128Int:'
    call randSeedGenerator(12345679)
    call system_clock(t1, clock_rate, clock_max)
    do i=1,randArrayN
      randArrayInt(i) = xorshift128Int(101)
    end do
    call system_clock(t2, clock_rate, clock_max)
    write(*,'(4x,a,1x,i0,1x,a,1x,es12.5,1x,a)') 'Time for', randArrayN, ' numbers =', dble(t2-t1)/dble(clock_rate), 'sec'
    call stdDev(randArrayInt,randMean,randStdDev)
    write(*,'(3x,2(1x,a,1x,es12.5),2(1x,a,1x,i0))') 'mean =', randMean, 'stdDev =', randStdDev, &
              'min =', minval(randArrayInt,1), 'max =', maxval(randArrayInt,1)
    write(*,*)

    write(*,'(2x,a)') 'mt19937:'
    call randSeedGenerator(12345678)
    call system_clock(t1, clock_rate, clock_max)
    do i=1,randArrayN
      randArray(i) = mt19937()
    end do
    call system_clock(t2, clock_rate, clock_max)
    write(*,'(4x,a,1x,i0,1x,a,1x,es12.5,1x,a)') 'Time to generate', randArrayN, 'numbers =', dble(t2-t1)/dble(clock_rate), ' sec'
    call stdDev(randArray,randMean,randStdDev)
    write(*,'(3x,4(1x,a,1x,es12.5))') 'mean =', randMean, 'stdDev =', randStdDev, &
              'min =', minval(randArray,1), 'max =', maxval(randArray,1)
    call chiSqTest(randArray, chiSqDOF+1, 0.d0, 1.d0, chiSq)
    write(*,'(4x,a,1x,es12.5,1x,a,1x,i0,1x,a)') 'chi-squared value:', chiSq, 'with', chiSqDOF, 'degrees of freedom'
    write(*,*)

  end subroutine

end module