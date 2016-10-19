!*********************************************************
!>
! A number of constants used in math and physics
!

module nfConstants
  implicit none

  integer, parameter :: sp = selected_real_kind(6, 37)
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: qp = selected_real_kind(33, 4931)

  ! Math constants

  ! Pi
  real(dp), parameter :: pi_ = 2.d0*acos(0.d0)

  ! Eulers number
  real(dp), parameter :: e_ = exp(1.d0)

  ! Euler-Mascheroni constant
  real(dp), parameter :: em_ = 0.57721566490153286d0

  ! Complex i
  complex(dp), parameter :: i_ = (0,1)

  ! Physics constants

  ! Speed of light in m/s
  real(dp), parameter :: c_ = 2.99792458d8
  ! Speed of light in cm/s (cgs units)
  real(dp), parameter :: c_cgs_ = 2.99792458d10

  ! Planck constant in J*s
  real(dp), parameter :: h_ = 6.626070040d-34
  ! Reduced Planck constant, h/2*pi, in J*s
  real(dp), parameter :: hbar_ = 1.054571800d-34
  ! Planck constant in eV*s
  real(dp), parameter :: h_eV_ = 4.135667662d-15
  ! Reduced Planck constant, \(h/2\pi\), in \(eV\) \(s\)
  real(dp), parameter :: hbar_eV_ = 6.582119514d-16
  ! Planck constant in \(erg\) \(s\) (cgs units)
  real(dp), parameter :: h_cgs_ = 6.626070040d-27
  ! Reduced Planck constant, \(h/2\pi\), in \(erg\) \(s\) (cgs units)
  real(dp), parameter :: hbar_cgs_ = 1.054571800d-27
  ! Reduced Planck constant, \(hc/2\pi\), in \(MeV\) \(fm\)
  real(dp), parameter :: hbarc_ = 197.3269788d0

  ! Gravitational constant in m^3 kg^-1 s^-2
  real(dp), parameter :: G_ = 6.6743d-11
  ! Gravitational constant in cm^3 g^-1 s^-2 (dyne cm^2 g^-2) (cgs units)
  real(dp), parameter :: G_cgs_ = 6.6743d-8

  ! Mass of electron in g (cgs units)
  real(dp), parameter :: me_ = 9.10938356d-28
  ! Mass of muon in g (cgs units)
  real(dp), parameter :: mmu_ = 1.883531594d-25
  ! Mass of tau in g (cgs units)
  real(dp), parameter :: mtau_ = 3.16747d-24
  ! Mass of proton in g (cgs units)
  real(dp), parameter :: mp_ = 1.672621898d-24
  ! Mass of neutron in g (cgs units)
  real(dp), parameter :: mn_ = 1.674927471d-24
  ! Mass of alpha particle in g (cgs units)
  real(dp), parameter :: malpha_ = 6.644657230d-24
  ! Atomic mass unit in g (cgs units)
  real(dp), parameter :: amu_ = 1.660539d-24

  ! Bohr radius in cm (cgs units)
  real(dp), parameter :: bohrr_ = 5.2917721067d-9
  ! Fine-structure constant
  real(dp), parameter :: alpha_ = 7.2973525664d-3
  ! Hartree energy in J
  real(dp), parameter :: eHartree_ = 4.359744650d-18
  ! Hartree energy in eV
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
  ! Wien displacement law constant in m K
  real(dp), parameter :: wiend_ = 2.8977729d-3

  ! Electric charge in C
  real(dp), parameter :: ec_ = 1.6021766208d-19
  ! Electric charge in cm^{3/2} g^{1/2} s^-1 (cgs units)
  real(dp), parameter :: ec_cgs_ = 4.8032d-10
  ! Electron volt in J
  real(dp), parameter :: ev_ = 1.6021766208d-19
  ! Electron volt in erg (cgs units)
  real(dp), parameter :: ev_cgs_ = 1.60217733d-12

  ! Electric constant (vacuum permittivity) in F m^-1
  real(dp), parameter :: epsilon0_ = 8.854187817d-12
  ! Magnetic constant (vacuum permeability) in N A^-2
  real(dp), parameter :: mu0_ = 1.2566370614d-6

  ! Solar mass in g
  real(dp), parameter :: msolar_ = 1.988435d33
  ! Solar radius in cm
  real(dp), parameter :: rsolar_ = 6.955d10
  ! Solar luminosity in erg s^-1
  real(dp), parameter :: lumsolar_ = 3.839d33
  ! Astronomical Unit in cm
  real(dp), parameter :: au_ = 1.49597870700d13

end module nfConstants