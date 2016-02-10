!*********************************************************
!>
! Sorting Routines

module sorting
  implicit none

  private
  integer, parameter :: dp = kind(0.d0)

  public :: swap
  public :: sort

  interface swap
    module procedure swap_i, swap_r, swap_d
  end interface

  interface sort
    module procedure sort_i
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

  subroutine sort_i(array)
    implicit none
    integer :: i, j, k, l, m, n, iStack
    integer, parameter :: sizeSub = 15, sizeStack = 50
    integer :: arrValue, array(:), arrStack(sizeStack)
    n = size(array,1)
    iStack = 0
    l = 1
    m = n
    do
      if(m-1.lt.sizeSub) then
        do j=l+1,m
          arrValue = array(j)
          do i=j-1,l,-1
            if(array(i).le.arrValue) exit
            array(i+1) = array(i)
          end do
          array(i+1) = arrValue
        end do
        if(iStack.eq.0) return
        m = arrStack(iStack)
        l = arrStack(iStack-1)
        iStack = iStack-2
      else
        k=(l+m)/2
        call swap(array(k),array(l+1))
        call swap(array(l),array(m),array(l).gt.array(m))
        call swap(array(l+1),array(m),array(l+1).gt.array(m))
        call swap(array(l),array(l+1),array(l).gt.array(l+1))
        i = l+1
        j = m
        arrValue = array(l+1)
        do
          do
            i=i+1
            if(array(i).ge.arrValue) exit
          end do
          do
            j=j-1
            if(array(j).le.arrValue) exit
          end do
          if(j.lt.i) exit
          call swap(array(i),array(j))
        end do
        array(l+1) = array(j)
        array(j) = arrValue
        iStack = iStack+2
        if(iStack.gt.sizeStack) write(*,*) 'Sort Warning: sizeStack too small!'
        if(m-i+1 .ge. j-l) then
          arrStack(iStack) = m
          arrStack(iStack-1) = i
          m = j-1
        else
          arrStack(iStack) = j-1
          arrStack(iStack-1) = l
          l = i
        end if
      end if
    end do
  end subroutine



end module