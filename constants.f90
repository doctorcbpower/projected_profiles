module constants
  real(kind=8), parameter :: pi=3.14159
  real(kind=8), parameter :: Mpc_in_m=3.086d22
  real(kind=8), parameter :: Msol_in_kg=1.989d30
  real(kind=8), parameter :: kms_in_ms=1000.
  real(kind=8), parameter :: lunit=1d-3,vunit=1,munit=1.d10
  real(kind=8), parameter :: gmks=6.67d-11
  !  real(kind=8), parameter :: 
  real(kind=8), parameter :: ggdt=gmks*(munit*Msol_in_kg)/(lunit*Mpc_in_m)/(vunit*kms_in_ms)/(vunit*kms_in_ms)
  real(kind=8), parameter :: kb=1.38d-23  ! Boltzmann constant in MKS
  real(kind=8), parameter :: mp=1.67d-27  ! Proton mass in kg
  real(kind=8), parameter :: gamma=5./3.
  real(kind=8), parameter :: u_to_temp=((gamma-1)*mp/kb)*(kms_in_ms)**2
  real(kind=8), parameter :: joules_to_ergs=1d7
  real(kind=8), parameter :: rhocrit0=27.755d-9 ! Critical density at z=0
  real(kind=8) :: molecular_weight
  real(kind=8), parameter :: deltavir=200.!101.14  
end module constants

