!*********************************************************
!>
! Different random number generators:
!
!  1. XORSHIFT:
!     * Random float from (0,1):
!        * xorshift
!        * xorshift64star
!        * xorshift128
!     * Integer from [0,a):
!        * xorshiftInt
!  2. Mersenne Twister:
!     * Random float from (0,1):
!        * mt19937

module randomNumbers
  use, intrinsic :: iso_fortran_env
  use parameters
  implicit none

  integer :: x, y, z, w, v

  integer, parameter :: wMT = 32
  integer, parameter :: nMT = 624
  integer, parameter :: mMT = 397
  integer, parameter :: rMT = 31
  integer, parameter :: uMT = 11
  integer, parameter :: sMT = 7
  integer, parameter :: tMT = 15
  integer, parameter :: lMT = 18
  integer :: mti = nMT+1
  integer :: mt(nMT)

  private :: x, y, z, w, v
  private :: wMT,nMT,mMT,rMT
  private :: uMT,sMT,tMT,lMT
  private :: mti,mt

contains

  subroutine randSeedGenerator(seed)

    integer :: seed
      !! input: seed is an integer for random number
      !! generators
    integer :: i

    x = seed
    y = ishft(x,2)
    z = ishft(y,-5)
    w = ishft(z,-2)
    v = ishft(w,3)

    mt(1) = seed
    do i=2,nMT
      mt(i) = 1812433253*(ieor(mt(i-1),ishft(mt(i-1),-30))+i)
    end do
    mti = nMT+1

  end subroutine randSeedGenerator

  real(dp) function xorshift()
    !! date: January 11, 2016
    !! version: v0.2
    !!
    !! Generic xorshift. Unknown period?. Returns a value
    !! from (0,1)

    integer :: t

    t = ieor(x,ishft(x,-7))
    x = y
    y = z
    z = w
    w = v
    v = ieor(ieor(v,ishft(v,6)),ieor(t,ishft(t,13)))
    xorshift = abs(dble(((y+y+1)*v))/dble(huge(y)))

  end function xorshift

  real(dp) function xorshift64star()
    !! date: January 11, 2016
    !! version: v0.2
    !!
    !! xorshift64s takes a xorshift generator and applies
    !! an invertible multiplication to its output. Has a
    !! maximal period of \(2^{64}-1\). Returns a value
    !! from (0,1)

    x = ieor(x,ishft(x,-12))
    x = ieor(x,ishft(x,25))
    x = ieor(x,ishft(x,-27))
    xorshift64star = abs(dble(x)/dble(huge(x)))

  end function xorshift64star

  real(dp) function xorshift128()
    !! date: January 11, 2016
    !! version: v0.2
    !!
    !! xorshift128 has a maximal period of \(2^{128}-1\).
    !! Returns a value from (0,1)

    integer :: t

    t = ieor(x,ishft(x,11))
    x = y
    y = z
    z = w
    w = ieor(w,ieor(ishft(w,-19),ieor(t,ishft(t,-8))))
    xorshift128 = abs(dble(w)/dble(huge(x)))

  end function xorshift128

  integer function xorshiftInt(a)
    !! date: January 10, 2016
    !! version: v0.1
    !!
    !! xorshiftInt returns a random integer from
    !! [0,a-1]. Has an unknown period?

    integer :: a
      !! input: upper bound for random integer

    xorshiftInt = floor(xorshift()*dble(a))

  end function xorshiftInt

  integer function xorshift64starInt(a)
    !! date: January 10, 2016
    !! version: v0.1
    !!
    !! xorshift64starInt returns a random integer
    !! from [0,a-1]. Has a
    !! maximal period of \(2^{64}-1\)

    integer :: a
      !! input: upper bound for random integer

    xorshift64starInt = floor(xorshift64star()*dble(a))

  end function xorshift64starInt

  integer function xorshift128Int(a)
    !! date: January 10, 2016
    !! version: v0.1
    !!
    !! xorshift128Int returns a random integer
    !! from [0,a-1]. Has a maximal period
    !! of \(2^{128}-1\)

    integer :: a
      !! input: upper bound for random integer

    xorshift128Int = floor(xorshift128()*dble(a))

  end function xorshift128Int

  real(dp) function mt19937()

    integer :: i, y

    if(mti.ge.nMT) then !retwist
      do i=1,nMT
        y = iand(mt(i),int(z'80000000',int32)) + iand(mt(mod(i+1,nMT)),int(z'7fffffff',int32))
        mt(i) = ieor(mt(mod(i+mMT,nMT)),ishft(y,-1))
        if(mod(y,2).ne.0) then
          mt(i) = ieor(mt(i),int(z'9908b0df',int32))
        end if
      end do
      mti = 1
    end if

    y = mt(mti)
    y = ieor(y,ishft(y,-uMT))
    y = ieor(y,iand(ishft(y,sMT),int(z'9d2c5680',int32)))
    y = ieor(y,iand(ishft(y,-tMT),int(z'efc60000',int32)))
    y = ieor(y,ishft(y,lMT))
    mti = mti + 1

    mt19937 = abs(dble(y)/dble(huge(y)))

  end function mt19937

end module randomNumbers
