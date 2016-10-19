!*********************************************************
!>
! Sorting Routines

module sorting
  implicit none

  private
  integer, parameter :: dp = kind(0.d0)

  public :: swap
  public :: qsort

  interface swap
    module procedure swap_i, swap_r, swap_d, swap_Cr, swap_Cd
  end interface

  interface qsort
    module procedure qsort_i, qsort_r, qsort_d
  end interface

contains

  subroutine swap_i(a,b,condition)
    implicit none
    integer :: a, b, aDum, bDum
    logical :: ifSwap
    logical, optional :: condition
    if(present(condition)) then
      ifSwap = condition
    else
      ifSwap = .true.
    end if
    if(ifSwap) then
      aDum = a; bDum = b
      a = bDum; b = aDum
    end if
  end subroutine

  subroutine swap_r(a,b,condition)
    implicit none
    real :: a, b, aDum, bDum
    logical :: ifSwap
    logical, optional :: condition
    if(present(condition)) then
      ifSwap = condition
    else
      ifSwap = .true.
    end if
    if(ifSwap) then
      aDum = a; bDum = b
      a = bDum; b = aDum
    end if
  end subroutine

  subroutine swap_d(a,b,condition)
    implicit none
    real(dp) :: a, b, aDum, bDum
    logical :: ifSwap
    logical, optional :: condition
    if(present(condition)) then
      ifSwap = condition
    else
      ifSwap = .true.
    end if
    if(ifSwap) then
      aDum = a; bDum = b
      a = bDum; b = aDum
    end if
  end subroutine

  subroutine swap_Cr(a,b,condition)
    implicit none
    complex :: a, b, aDum, bDum
    logical :: ifSwap
    logical, optional :: condition
    if(present(condition)) then
      ifSwap = condition
    else
      ifSwap = .true.
    end if
    if(ifSwap) then
      aDum = a; bDum = b
      a = bDum; b = aDum
    end if
  end subroutine

  subroutine swap_Cd(a,b,condition)
    implicit none
    complex(dp) :: a, b, aDum, bDum
    logical :: ifSwap
    logical, optional :: condition
    if(present(condition)) then
      ifSwap = condition
    else
      ifSwap = .true.
    end if
    if(ifSwap) then
      aDum = a; bDum = b
      a = bDum; b = aDum
    end if
  end subroutine

  subroutine qsort_i(array)
    implicit none
    integer :: k, n, left, right, marker
    integer :: pivot, array(:)
    n = size(array,1)
    if(n>1) then
      k = (n/2)
      pivot = array(k)
      left = 0
      right = n+1
      do while(left<right)
        right = right-1
        do while(array(right).gt.pivot)
          right = right-1
        end do
        left = left+1
        do while(array(left).lt.pivot)
          left = left+1
        end do
        if(left<right) then
          call swap(array(left),array(right))
        end if
      end do
      if(left==right) then
        marker = left +1
      else
        marker = left
      end if
      call qsort(array(:marker-1))
      call qsort(array(marker:))
    end if
  end subroutine

  subroutine qsort_r(array)
    implicit none
    integer :: k, n, left, right, marker
    real :: pivot, array(:)
    n = size(array,1)
    if(n>1) then
      k = (n/2)
      pivot = array(k)
      left = 0
      right = n+1
      do while(left<right)
        right = right-1
        do while(array(right).gt.pivot)
          right = right-1
        end do
        left = left+1
        do while(array(left).lt.pivot)
          left = left+1
        end do
        if(left<right) then
          call swap(array(left),array(right))
        end if
      end do
      if(left==right) then
        marker = left +1
      else
        marker = left
      end if
      call qsort(array(:marker-1))
      call qsort(array(marker:))
    end if
  end subroutine

  subroutine qsort_d(array)
    implicit none
    integer :: k, n, left, right, marker
    real(dp) :: pivot, array(:)
    n = size(array,1)
    if(n>1) then
      k = (n/2)
      pivot = array(k)
      left = 0
      right = n+1
      do while(left<right)
        right = right-1
        do while(array(right).gt.pivot)
          right = right-1
        end do
        left = left+1
        do while(array(left).lt.pivot)
          left = left+1
        end do
        if(left<right) then
          call swap(array(left),array(right))
        end if
      end do
      if(left==right) then
        marker = left +1
      else
        marker = left
      end if
      call qsort(array(:marker-1))
      call qsort(array(marker:))
    end if
  end subroutine

end module