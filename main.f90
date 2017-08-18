program main
  use constants
  use structure, only: getmvir,getrvir
  implicit none
  
  real(kind=8) :: virial_mass,virial_velocity,virial_radius
  real(kind=8) :: nfw_concentration,virial_overdensity
  integer(kind=4) :: halo_type
  real(kind=8) :: disc_mass_fraction,disc_scale_length_fraction
  real(kind=8) :: disc_mass,disc_radius
  real(kind=8) :: bulge_mass_fraction,bulge_scale_length_fraction
  real(kind=8) :: bulge_mass,bulge_radius
  real(kind=8) :: black_hole_mass_fraction
  real(kind=8) :: black_hole_mass
  real(kind=8) :: velocity_anisotropy_amplitude,velocity_anisotropy_radius
  real(kind=8) :: fmax=1000.,fcut=0.5
  
  logical :: isbh                         ! Is a black hole present?

  integer(kind=4) :: nr                   ! Dummy index for looping over r,z arrays
  real(kind=8) :: rmax,rcut
  real(kind=8) :: lr,lrmin,lrmax,dlr   ! Log radius
  
  integer(kind=4) :: i,j,k,l,m,n,nn
  real(kind=8) :: sum

  character(kind=1,len=132) :: outfile,paramfile
  character(kind=1,len=10) :: istring
  logical:: fexist

  integer(kind=4) :: nr_tab,nr_max
  real(kind=8), allocatable :: r_tab(:)
  real(kind=8), allocatable :: sigr_tab(:)
  real(kind=8), allocatable :: sigp_tab(:)  
  real(kind=8), allocatable :: sigm_tab(:)
  
  isbh=.false.
  
  if(command_argument_count().eq.0) stop 'Usage: compute_sigma.exe <parameter_file>'
  
  call get_command_argument(1,paramfile)

  inquire(file=paramfile,exist=fexist)

  if(fexist.eqv..false.) stop 'Error: parameter file does not exist'
  call set_parameters(paramfile,virial_mass,virial_velocity,nfw_concentration,&
     & virial_overdensity,halo_type,disc_mass_fraction,disc_scale_length_fraction,&
     & bulge_mass_fraction,bulge_scale_length_fraction,black_hole_mass_fraction,&
     & velocity_anisotropy_amplitude,velocity_anisotropy_radius)
  if(black_hole_mass_fraction.gt.0.0) isbh=.true.

  ! Define dimension and bounds of the mesh

  outfile='projected_props.'
  if(halo_type.eq.0) then
     outfile=trim(outfile)//'hern_halo.'
  else
     outfile=trim(outfile)//'nfw_halo.'
  end if

  if(bulge_mass.gt.0.0) outfile=trim(outfile)//'hern_bulge.'

  if(velocity_anisotropy_radius.le.0.0) then
     if(velocity_anisotropy_amplitude.lt.0) then
        write(istring,'(f5.2)') velocity_anisotropy_amplitude
     else
        write(istring,'(f4.2)') velocity_anisotropy_amplitude
     endif
  else
     istring='var'
  end if

  outfile=trim(outfile)//trim(istring)//'.txt'

  if(virial_overdensity.eq.0.0) virial_overdensity=deltavir
  
  if(virial_velocity.gt.0.0) then
     if(virial_mass.ne.getmvir(virial_velocity,virial_overdensity,rhocrit0))&
          & virial_mass=getmvir(virial_velocity,virial_overdensity,rhocrit0)
  end if

  bulge_mass = bulge_mass_fraction
  bulge_radius = bulge_scale_length_fraction 
  disc_mass = disc_mass_fraction
  disc_radius = disc_scale_length_fraction   
  black_hole_mass = black_hole_mass_fraction

  virial_radius=getrvir(virial_mass,virial_overdensity,rhocrit0)  
  nr_tab=500
  lrmin=dlog10(virial_radius)-2.7
  lrmax=dlog10(virial_radius)+1.1
  dlr=(lrmax-lrmin)/real(nr_tab-1)
  
  allocate(r_tab(nr_tab))       ! Tabulated radius
  allocate(sigr_tab(nr_tab))    ! Tabulated mass
  allocate(sigp_tab(nr_tab))    ! Tabulated mass  
  allocate(sigm_tab(nr_tab))    ! Tabulated mass

  lr=lrmin
  
  do nr=1,nr_tab
     r_tab(nr)=10**lr
     lr=lr+dlr
  end do
  
  write(*,*) 'Radial limits [kpc] :',r_tab(1),r_tab(nr_tab)
  write(*,*)
  
  call compute_sigma(nr_tab,r_tab,&
       & virial_mass,virial_radius,nfw_concentration,&
       & virial_overdensity,halo_type,&
       & bulge_mass,bulge_radius,&
       & disc_mass,disc_radius,&
       & black_hole_mass,fmax,fcut,&
       & velocity_anisotropy_amplitude,velocity_anisotropy_radius,&
       & sigr_tab,sigp_tab,sigm_tab)
  
  open(32,file=outfile,status='unknown')
  write(32,*) virial_mass,virial_radius,virial_velocity
  write(32,*) bulge_mass,bulge_radius
  write(32,*) black_hole_mass
  write(32,*) velocity_anisotropy_amplitude,velocity_anisotropy_radius
  do i=1,nr_tab
     write(32,*) r_tab(i),sigr_tab(i),sigp_tab(i),sigm_tab(i)
  end do
  close(32)
  
end program main



