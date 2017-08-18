module variables
  real(kind=8) :: mbh=0.0     ! central black hole mass
  real(kind=8) :: radius  ! radius
  real(kind=8) :: m200=0.0  ! virial mass
  real(kind=8) :: mvir=0.0  ! virial mass    
  real(kind=8) :: mhalo=0.0 ! halo mass
  real(kind=8) :: v200=0.0  ! virial velocity
  integer(kind=4) :: halotype ! Type of halo
  real(kind=8) :: rs      ! NFW scale radius
  real(kind=8) :: r200    ! virial radius
  real(kind=8) :: c200    ! NFW concentration
  real(kind=8) :: sigma0  ! central surface density
  real(kind=8) :: rdisc=0.0   ! disc scale length
  real(kind=8) :: bdisc=0.0   ! disc shape parameter (Miyamoto-Nagai disc)
  real(kind=8) :: mbulge=0.0  ! bulge mass
  real(kind=8) :: fbulge=0.0  ! bulge extent
  real(kind=8) :: rbulge=0.0  ! bulge extent
  real(kind=8) :: mstar=0.0   ! stellar mass
  real(kind=8) :: mdisc=0.0   ! disc mass
  real(kind=8) :: fdisc=0.0   ! stellar mass
  real(kind=8) :: beta0
  real(kind=8) :: r_a

  integer(kind=4) :: nrbin
  real(kind=8), allocatable :: rtab(:)
  real(kind=8), allocatable :: sigma_r(:)
  real(kind=8), allocatable :: sigma_p(:)
  real(kind=8), allocatable :: sigma_m(:)
end module variables


