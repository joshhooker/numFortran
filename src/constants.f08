!*********************************************************
!>
! A number of constants used in math and physics
!

module constants
  implicit none

  private :: dp
  integer, parameter :: dp = kind(0.d0)

  ! Math constants

  ! pi
  real(dp), parameter :: pi_ = 2.d0*acos(0.d0)

  ! eulers number
  real(dp), parameter :: e_ = exp(1.d0)

  ! euler-mascheroni constant
  real(dp), parameter :: em_ = 0.57721566490153286d0

  ! complex i
  complex(dp), parameter :: i_ = (0,1)

  ! Physics constants

  ! speed of light in m/s
  real(dp), parameter :: c_ = 2.99792458d8
  ! speed of light in cm/s (cgs units)
  real(dp), parameter :: c_cgs_ = 2.99792458d10

  ! planck constant in J*s
  real(dp), parameter :: h_ = 6.626070040d-34
  ! reduced planck constant, h/2*pi, in J*s
  real(dp), parameter :: hbar_ = 1.054571800d-34
  ! planck constant in eV*s
  real(dp), parameter :: h_eV_ = 4.135667662d-15
  ! reduced planck constant, \(h/2\pi\), in \(eV\) \(s\)
  real(dp), parameter :: hbar_eV_ = 6.582119514d-16
  ! planck constant in \(erg\) \(s\) (cgs units)
  real(dp), parameter :: h_cgs_ = 6.626070040d-27
  ! reduced planck constant, \(h/2\pi\), in \(erg\) \(s\) (cgs units)
  real(dp), parameter :: hbar_cgs_ = 1.054571800d-27
  ! reduced planck constant, \(hc/2\pi\), in \(MeV\) \(fm\)
  real(dp), parameter :: hbarc_ = 197.3269788d0

  ! gravitational constant in m^3 kg^-1 s^-2
  real(dp), parameter :: G_ = 6.6743d-11
  ! gravitational constant in cm^3 g^-1 s^-2 (dyne cm^2 g^-2) (cgs units)
  real(dp), parameter :: G_cgs_ = 6.6743d-8

  ! mass of electron in g (cgs units)
  real(dp), parameter :: me_ = 9.10938356d-28
  ! mass of proton in g (cgs units)
  real(dp), parameter :: mp_ = 1.672621898d-24
  ! mass of neutron in g (cgs units)
  real(dp), parameter :: mn_ = 1.674927471d-24
  ! mass of alpha particle in g (cgs units)
  real(dp), parameter :: malpha_ = 6.644657230d-24
  ! atomic mass unit in g (cgs units)
  real(dp), parameter :: amu_ = 1.660539d-24

  ! Bohr radius in cm (cgs units)
  real(dp), parameter :: bohrr_ = 5.2917721067d-9
  ! fine-structure constant
  real(dp), parameter :: alpha_ = 7.2973525664d-3
  ! hartree energy in J
  real(dp), parameter :: eHartree_ = 4.359744650d-18
  ! hartree energy in eV
  real(dp), parameter :: eHartree_eV_ = 27.21138602d0

  ! Avagadro's number
  real(dp), parameter :: na_ = 6.0221367d23
  ! Boltzmann constant in J K^-1
  real(dp), parameter :: k_ = 1.38064852d-23
  ! Faraday constant in C mol^-1
  real(dp), parameter :: F_ = 96485.33289d0
  ! Rydberg constant in m^-1
  real(dp), parameter :: ryd_ = 10973731.568508d0
  ! Stefan-Boltzmann constant in W m^-2 K^-4
  real(dp), parameter :: stefboltz_ = 5.670367d-8

  ! electric charge in C
  real(dp), parameter :: ec_ = 1.6021766208d-19
  ! electric charge in cm^{3/2} g^{1/2} s^-1 (cgs units)
  real(dp), parameter :: ec_cgs_ = 4.8032d-10
  ! electron volt in J
  real(dp), parameter :: ev_ = 1.6021766208d-19
  ! electron volt in erg (cgs units)
  real(dp), parameter :: ev_cgs_ = 1.60217733d-12

  ! electric constant (vacuum permittivity) in F m^-1
  real(dp), parameter :: epsilon0_ = 8.854187817d-12
  ! magnetic constant (vacuum permeability) in N A^-2
  real(dp), parameter :: mu0_ = 1.2566370614d-6

end module constants
