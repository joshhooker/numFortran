module testRandNum
  use randomNumbers

  use iso_fortran_env
  implicit none

  public :: randNumTest

contains

  subroutine randNumTest()
    implicit none

    write(*,'(a)') '*************************************'
    write(*,'(a)') 'TESTING RANDOM NUMBER GENERATORS'

  end subroutine

end module