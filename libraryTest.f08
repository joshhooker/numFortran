program libraryTest
  use constants
  use cubicSpline
  use integration
  use odeSolver
  use randomNumbers
  use specialFunctions
  use statTests

  use iso_fortran_env
  implicit none

  integer, parameter :: dp = kind(0.d0)

  integer :: i

  ! Variables to time functions
  integer(int64) :: t1, t2, clock_rate, clock_max

  ! Variables to test cubic spline
  real(dp) :: x, y
  real(dp) :: splineX(10), splineY(10)
  real(dp) :: splineA(10), splineB(10), splineC(10), splineD(10)

  ! Variables to test ode solver
  real(dp) :: h, h0, odeResult

  ! Variables to test random number generators
  integer, parameter :: randArrayN = 5000000
  integer :: randArrayInt(randArrayN)
  real(dp) :: randArray(randArrayN)
  real(dp) :: randMean, randStdDev

  ! Variables to test numerical integration
  integer :: intN
  real(dp) :: intA, intB, intResult
  real(dp), allocatable :: intConsts(:)

  ! Variables to test special functions
  integer :: sfTestsN
  real(dp) :: sfResult, sfSumDiff
  real(dp), allocatable :: testBesselN(:), testBesselX(:), testBessel(:)
  real(dp), allocatable :: testAiryX(:), testAiry(:)

  !*******************************
  ! CONSTANT TESTING
  !*******************************
  write(*,'(a)') '*************************************'
  write(*,'(a)') 'TESTING CONSTANTS:'
  write(*,'(2x,a,1x,es25.18)') 'pi =', pi_
  write(*,'(2x,a,1x,es25.18)') 'e =', e_
  write(*,'(2x,a,1x," ",2f8.5)') 'i =', i_
  write(*,*)
  write(*,*)

  !*******************************
  ! CUBIC SPLINE TESTING
  !*******************************
  write(*,'(a)') '*************************************'
  ! Test cubic spline of function y = cos(sin(x))
  write(*,'(a)') 'TESTING CUBIC SPLINE:'
  write(*,'(a)') 'Function: y = cos(sin(x))'
  write(*,*)
  write(*,'(a)') 'Starting values for the spline:'
  do i=1,10
    splineX(i) = 1.3d0 + (i-1)*0.73d0
    splineY(i) = cos(sin(splineX(i)))
    write(*,'(1x,2(1x,a,1x,es12.5))') 'x =', splineX(i), 'y =', splineY(i)
  end do
  call naturalCubicSpline(splineX, splineY, splineA, splineB, splineC, splineD)
  x = minval(splineX,1)+0.1d0
  write(*,*)
  write(*,'(a)') 'Values calculated from spline:'
  do while (x<= maxval(splineX,1))
    y = computeNaturalCubicSpline(splineX, x, splineA, splineB, splineC, splineD)
    write(*,'(1x,3(1x,a,1x,es12.5))') 'x =', x, 'y =', y, 'diff =', y-cos(sin(x))
    x = x + 0.23d0
  end do
  write(*,*)
  write(*,*)

  !*******************************
  ! ODE SOLVER TESTING
  !*******************************
  write(*,'(a)') '*************************************'
  write(*,'(a)') 'TESTING 1ST ORDER ODE SOLVERS:'
  write(*,'(2x,a)') "Function: y' = y^2 + 1"
  write(*,'(2x,a)') 'From 0 to 1.5 with y(0) = 0'
  write(*,'(2x,a)') 'Solution function: tan(x)'
  write(*,'(2x,a,1x,es12.5)') 'Numerical solution: tan(1.5) =', tan(1.5d0)
  write(*,*)
  write(*,'(4x,a)') '4th order RK Fixed Step:'
  do i=1,10
    h = 0.2/i
    h0 = h
    call system_clock(t1, clock_rate, clock_max)
    odeResult = rk1FixedStep(0.d0,1.5d0,h,0.d0)
    call system_clock(t2, clock_rate, clock_max)
    write(*,'(6x,a,1x,f6.4,3(1x,a,1x,es12.5),1x,a)') 'h =', h0, 'result =', odeResult, &
      'diff =', odeResult-tan(1.5d0), ' time =', dble(t2-t1)/dble(clock_rate), 'sec'
  end do
  write(*,*)
  write(*,'(4x,a)') 'RK Adaptive Step:'
  call system_clock(t1, clock_rate, clock_max)
  odeResult = rk1AdaptStep(0.d0,1.5d0,0.d0)
  call system_clock(t2, clock_rate, clock_max)
  write(*,'(5x,3(1x,a,1x,es12.5),1x,a)') 'result =', odeResult, &
      'diff =', odeResult-tan(1.5d0), ' time =', dble(t2-t1)/dble(clock_rate), 'sec'
  write(*,*)
  write(*,*)

  !*******************************
  ! RANDOM NUMBER TESTING
  !*******************************
  write(*,'(a)') '*************************************'
  write(*,'(a,i0,1x,a)') 'TESTING RANDOM NUMBER GENERATORS (', randArrayN, 'numbers):'

  write(*,'(2x,a)') 'xorshift:'
  call randSeedGenerator(12345679)
  call system_clock(t1, clock_rate, clock_max)
  do i=1,randArrayN
    randArray(i) = xorshift()
  end do
  call system_clock(t2, clock_rate, clock_max)
  write(*,'(4x,a,1x,i0,1x,a,1x,es12.5,1x,a)') 'Time for', randArrayN, 'numbers =', dble(t2-t1)/dble(clock_rate), ' sec'
  call stdDev(randArray,randMean,randStdDev)
  write(*,'(3x,4(1x,a,1x,es12.5))') 'mean =', randMean, 'stdDev =', randStdDev, &
            'min =', minval(randArray,1), 'max =', maxval(randArray,1)
  write(*,*)

  write(*,'(2x,a)') 'xorshift64star:'
  call randSeedGenerator(12345679)
  call system_clock(t1, clock_rate, clock_max)
  do i=1,randArrayN
    randArray(i) = xorshift64star()
  end do
  call system_clock(t2, clock_rate, clock_max)
  write(*,'(4x,a,1x,i0,1x,a,1x,es12.5,1x,a)') 'Time for', randArrayN, 'numbers =', dble(t2-t1)/dble(clock_rate), 'sec'
  call stdDev(randArray,randMean,randStdDev)
  write(*,'(3x,4(1x,a,1x,es12.5))') 'mean =', randMean, 'stdDev =', randStdDev, &
            'min =', minval(randArray,1), 'max =', maxval(randArray,1)
  write(*,*)

  write(*,'(2x,a)') 'xorshift128:'
  call randSeedGenerator(12345679)
  call system_clock(t1, clock_rate, clock_max)
  do i=1,randArrayN
    randArray(i) = xorshift128()
  end do
  call system_clock(t2, clock_rate, clock_max)
  write(*,'(4x,a,1x,i0,1x,a,1x,es12.5,1x,a)') 'Time for', randArrayN, 'numbers =', dble(t2-t1)/dble(clock_rate), 'sec'
  call stdDev(randArray,randMean,randStdDev)
  write(*,'(3x,4(1x,a,1x,es12.5))') 'mean =', randMean, 'stdDev =', randStdDev, &
            'min =', minval(randArray,1), 'max =', maxval(randArray,1)
  write(*,*)

  write(*,'(2x,a)') 'xorshiftInt:'
  call randSeedGenerator(12345679)
  call system_clock(t1, clock_rate, clock_max)
  do i=1,randArrayN
    randArrayInt(i) = xorshiftInt(101)
  end do
  call system_clock(t2, clock_rate, clock_max)
  write(*,'(4x,a,1x,i0,1x,a,1x,es12.5,1x,a)') 'Time for', randArrayN, ' numbers =', dble(t2-t1)/dble(clock_rate), 'sec'
  call stdDevInt(randArrayInt,randMean,randStdDev)
  write(*,'(3x,2(1x,a,1x,es12.5),2(1x,a,1x,i0))') 'mean =', randMean, 'stdDev =', randStdDev, &
            'min =', minval(randArrayInt,1), 'max =', maxval(randArrayInt,1)
  write(*,*)

  write(*,'(2x,a)') 'xorshift64starInt:'
  call randSeedGenerator(12345679)
  call system_clock(t1, clock_rate, clock_max)
  do i=1,randArrayN
    randArrayInt(i) = xorshift64starInt(101)
  end do
  call system_clock(t2, clock_rate, clock_max)
  write(*,'(4x,a,1x,i0,1x,a,1x,es12.5,1x,a)') 'Time for', randArrayN, ' numbers =', dble(t2-t1)/dble(clock_rate), 'sec'
  call stdDevInt(randArrayInt,randMean,randStdDev)
  write(*,'(3x,2(1x,a,1x,es12.5),2(1x,a,1x,i0))') 'mean =', randMean, 'stdDev =', randStdDev, &
            'min =', minval(randArrayInt,1), 'max =', maxval(randArrayInt,1)
  write(*,*)

  write(*,'(2x,a)') 'xorshift128Int:'
  call randSeedGenerator(12345679)
  call system_clock(t1, clock_rate, clock_max)
  do i=1,randArrayN
    randArrayInt(i) = xorshift128Int(101)
  end do
  call system_clock(t2, clock_rate, clock_max)
  write(*,'(4x,a,1x,i0,1x,a,1x,es12.5,1x,a)') 'Time for', randArrayN, ' numbers =', dble(t2-t1)/dble(clock_rate), 'sec'
  call stdDevInt(randArrayInt,randMean,randStdDev)
  write(*,'(3x,2(1x,a,1x,es12.5),2(1x,a,1x,i0))') 'mean =', randMean, 'stdDev =', randStdDev, &
            'min =', minval(randArrayInt,1), 'max =', maxval(randArrayInt,1)
  write(*,*)

  write(*,'(2x,a)') 'mt19937:'
  call randSeedGenerator(12349)
  call system_clock(t1, clock_rate, clock_max)
  do i=1,randArrayN
    randArray(i) = mt19937()
  end do
  call system_clock(t2, clock_rate, clock_max)
  write(*,'(4x,a,1x,i0,1x,a,1x,es12.5,1x,a)') 'Time for', randArrayN, ' numbers =', dble(t2-t1)/dble(clock_rate), 'sec'
  call stdDev(randArray,randMean,randStdDev)
  write(*,'(3x,4(1x,a,1x,es12.5))') 'mean =', randMean, 'stdDev =', randStdDev, &
            'min =', minval(randArray,1), 'max =', maxval(randArray,1)
  write(*,*)
  write(*,*)

  !*******************************
  ! INTEGRAL TESTING
  !*******************************
  write(*,'(a)') '*************************************'
  write(*,'(a)') 'TESTING INTEGRALS:'
  call readGLParameters()

  intA = 0.d0
  intB = 54.d0
  intN = 0
  allocate(intConsts(intN))
  write(*,'(1x,2(1x,a,1x,e10.4))') 'Function = tan(x)*cos(x) from ', intA, 'to', intB
  write(*,'(2x,a)') 'Result = cos(a)-cos(b)'
  write(*,*)

  write(*,'(2x,a)') 'Gauss-Legendre Integration:'
  intResult = gaussLegendre(integralFunction1,intA,intB,intN,intConsts)
  write(*,'(3x,2(1x,a,1x,e10.4))') 'Result =', intResult, ', Diff =', intResult-(cos(intA)-cos(intB))
  deallocate(intConsts)
  write(*,*)

  intA = -1.d0
  intB = 1.d0
  intN = 3
  allocate(intConsts(intN))
  intConsts(1) = 2.d0
  intConsts(2) = 3.d0
  intConsts(3) = 2.d0
  write(*,'(2x,a)') 'Function = 2x^3 - 3x^2*cos(2x) from -1 to 1'
  write(*,'(2x,a)') 'Result ~ -0.115505630597095'
  write(*,*)

  write(*,'(2x,a)') 'Gauss-Legendre Integration:'
  intResult = gaussLegendre(integralFunction2,intA,intB,intN,intConsts)
  write(*,'(3x,2(1x,a,1x,e10.4))') 'Result =', intResult, ', Diff =', &
      intResult-(-0.115505630597095)
  deallocate(intConsts)
  write(*,*)
  write(*,*)

  !*******************************
  ! SPECIAL FUNCTION TESTING
  !*******************************
  write(*,'(a)') '*************************************'
  write(*,'(a)') 'TESTING SPECIAL FUNCTIONS:'

  write(*,'(2x,a)') 'TESTING BESSEL FUNCTION OF THE FIRST KIND:'
  open(unit=1,file='tests/test_jn.dat',status='old')
  read(1,*) sfTestsN
  allocate(testBesselN(sfTestsN), testBesselX(sfTestsN), testBessel(sfTestsN))
  do i=1,sfTestsN
    read(1,*) testBesselN(i), testBesselX(i), testBessel(i)
  end do
  close(1)
  sfSumDiff = 0.d0
  do i=1, sfTestsN
    sfSumDiff = sfSumDiff+abs(besselJ(testBesselN(i),testBesselX(i))-testBessel(i))
  end do
  write(*,'(3x,2(1x,a,1x,es16.8))') 'Total Error =', sfSumDiff, 'Average Error =', sfSumDiff/dble(sfTestsN)
  deallocate(testBesselN, testBesselX, testBessel)
  write(*,*)

  write(*,'(2x,a)') 'TESTING BESSEL FUNCTION OF THE SECOND KIND:'
  open(unit=1,file='tests/test_yn.dat',status='old')
  read(1,*) sfTestsN
  allocate(testBesselN(sfTestsN), testBesselX(sfTestsN), testBessel(sfTestsN))
  do i=1,sfTestsN
    read(1,*) testBesselN(i), testBesselX(i), testBessel(i)
  end do
  close(1)
  sfSumDiff = 0.d0
  do i=1, sfTestsN
    sfSumDiff = sfSumDiff+abs(besselY(testBesselN(i),testBesselX(i))-testBessel(i))
    write(*,*) testBesselN(i), testBesselX(i), besselY(testBesselN(i),testBesselX(i))
  end do
  write(*,'(3x,2(1x,a,1x,es16.8))') 'Total Error =', sfSumDiff, 'Average Error =', sfSumDiff/dble(sfTestsN)
  deallocate(testBesselN, testBesselX, testBessel)
  write(*,*)

  write(*,'(2x,a)') 'TESTING MODIFIED BESSEL FUNCTION OF THE FIRST KIND:'
  open(unit=1,file='tests/test_in.dat',status='old')
  read(1,*) sfTestsN
  allocate(testBesselN(sfTestsN), testBesselX(sfTestsN), testBessel(sfTestsN))
  do i=1,sfTestsN
    read(1,*) testBesselN(i), testBesselX(i), testBessel(i)
  end do
  close(1)
  sfSumDiff = 0.d0
  do i=1, sfTestsN
    sfSumDiff = sfSumDiff+abs(besselI(testBesselN(i),testBesselX(i))-testBessel(i))
  end do
  write(*,'(3x,2(1x,a,1x,es16.8))') 'Total Error =', sfSumDiff, 'Average Error =', sfSumDiff/dble(sfTestsN)
  deallocate(testBesselN, testBesselX, testBessel)
  write(*,*)

  write(*,'(2x,a)') 'TESTING MODIFIED BESSEL FUNCTION OF THE SECOND KIND:'
  open(unit=1,file='tests/test_kn.dat',status='old')
  read(1,*) sfTestsN
  allocate(testBesselN(sfTestsN), testBesselX(sfTestsN), testBessel(sfTestsN))
  do i=1,sfTestsN
    read(1,*) testBesselN(i), testBesselX(i), testBessel(i)
  end do
  close(1)
  sfSumDiff = 0.d0
  do i=1, sfTestsN
    sfSumDiff = sfSumDiff+abs(besselK(testBesselN(i),testBesselX(i))-testBessel(i))
  end do
  write(*,'(3x,2(1x,a,1x,es16.8))') 'Total Error =', sfSumDiff, 'Average Error =', sfSumDiff/dble(sfTestsN)
  deallocate(testBesselN, testBesselX, testBessel)
  write(*,*)

  write(*,'(2x,a)') 'TESTING SPHERICAL BESSEL FUNCTION OF THE FIRST KIND:'
  open(unit=1,file='tests/test_sphjn.dat',status='old')
  read(1,*) sfTestsN
  allocate(testBesselN(sfTestsN), testBesselX(sfTestsN), testBessel(sfTestsN))
  do i=1,sfTestsN
    read(1,*) testBesselN(i), testBesselX(i), testBessel(i)
  end do
  close(1)
  sfSumDiff = 0.d0
  do i=1, sfTestsN
    sfSumDiff = sfSumDiff+abs(sphBesselJ(int(testBesselN(i)),testBesselX(i))-testBessel(i))
  end do
  write(*,'(3x,2(1x,a,1x,es16.8))') 'Total Error =', sfSumDiff, 'Average Error =', sfSumDiff/dble(sfTestsN)
  deallocate(testBesselN, testBesselX, testBessel)
  write(*,*)

  ! write(*,'(2x,a)') 'TESTING SPHERICAL BESSEL FUNCTION OF THE SECOND KIND:'
  ! open(unit=1,file='tests/test_sphyn.dat',status='old')
  ! read(1,*) sfTestsN
  ! allocate(testBesselN(sfTestsN), testBesselX(sfTestsN), testBessel(sfTestsN))
  ! do i=1,sfTestsN
  !   read(1,*) testBesselN(i), testBesselX(i), testBessel(i)
  ! end do
  ! close(1)
  ! sfSumDiff = 0.d0
  ! do i=1, sfTestsN
  !   sfSumDiff = sfSumDiff+abs(sphBesselY(int(testBesselN(i)),testBesselX(i))-testBessel(i))
  ! end do
  ! write(*,'(3x,2(1x,a,1x,es16.8))') 'Total Error =', sfSumDiff, 'Average Error =', sfSumDiff/dble(sfTestsN)
  ! deallocate(testBesselN, testBesselX, testBessel)
  ! write(*,*)

  write(*,'(2x,a)') 'TESTING AIRY FUNCTION OF THE FIRST KIND:'
  open(unit=1,file='tests/test_airyA.dat',status='old')
  read(1,*) sfTestsN
  allocate(testAiryX(sfTestsN), testAiry(sfTestsN))
  do i=1,sfTestsN
    read(1,*) testAiryX(i), testAiry(i)
  end do
  close(1)
  sfSumDiff = 0.d0
  do i=1, sfTestsN
    sfSumDiff = sfSumDiff+abs(airyA(testAiryX(i))-testAiry(i))
  end do
  write(*,'(3x,2(1x,a,1x,es16.8))') 'Total Error =', sfSumDiff, 'Average Error =', sfSumDiff/dble(sfTestsN)
  deallocate(testAiryX, testAiry)
  write(*,*)

  write(*,'(2x,a)') 'TESTING AIRY FUNCTION OF THE SECOND KIND:'
  open(unit=1,file='tests/test_airyB.dat',status='old')
  read(1,*) sfTestsN
  allocate(testAiryX(sfTestsN), testAiry(sfTestsN))
  do i=1,sfTestsN
    read(1,*) testAiryX(i), testAiry(i)
  end do
  close(1)
  sfSumDiff = 0.d0
  do i=1, sfTestsN
    sfSumDiff = sfSumDiff+abs(airyB(testAiryX(i))-testAiry(i))
  end do
  write(*,'(3x,2(1x,a,1x,es16.8))') 'Total Error =', sfSumDiff, 'Average Error =', sfSumDiff/dble(sfTestsN)
  deallocate(testAiryX, testAiry)
  write(*,*)

  write(*,'(2x,a)') 'TESTING GAMMA FUNCTION:'
  print *, gamma(10.21d0)
  print *, gammaFunc(10.21d0)
  write(*,*)

contains

  real(dp) function integralFunction1(n,c,x)
    integer :: n
    real(dp) :: c(n)
    real(dp) :: x
    c = 0.d0
    integralFunction1 = tan(x)*cos(x)
  end function integralFunction1

  real(dp) function integralFunction2(n,c,x)
    integer :: n
    real(dp) :: c(n)
    real(dp) :: x
    integralFunction2 = c(1)*x*x*x-c(2)*x*x*cos(c(3)*x)
  end function integralFunction2

end program libraryTest
