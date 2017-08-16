subroutine set_parameters(infile,virial_mass,virial_velocity,nfw_concentration,&
     & virial_overdensity,disc_mass_fraction,disc_scale_length_fraction,&
     & bulge_mass_fraction,bulge_scale_length_fraction,black_hole_mass_fraction,&
     & velocity_anisotropy_amplitude,velocity_anisotropy_radius)

  implicit none
  
  character(kind=1,len=*) :: infile              ! input parameter file

  real(kind=8) :: virial_mass,virial_velocity,nfw_concentration,virial_overdensity
  real(kind=8) :: disc_mass_fraction,disc_scale_length_fraction
  real(kind=8) :: bulge_mass_fraction,bulge_scale_length_fraction
  real(kind=8) :: black_hole_mass_fraction
  real(kind=8) :: velocity_anisotropy_amplitude,velocity_anisotropy_radius

  ! Namelists

  ! Halo properties
  namelist /halo/ virial_mass,&   ! Virial mass, Msol
       & virial_velocity,&        ! Virial velocity, in km/s
       & nfw_concentration,&      ! NFW Concentration parameter
       & virial_overdensity       ! Overdensity parameter
  ! Disc properties
  namelist /disc/ disc_mass_fraction,&     ! Disc mass as fraction of virial mass
       & disc_scale_length_fraction        ! Disc scale length as fraction of halo
                                           ! scale radius
  ! Bulge properties
  namelist /bulge/ bulge_mass_fraction,&   ! Bulge mass as fraction of virial mass
       & bulge_scale_length_fraction       ! Bulge scale length as fraction of halo
                                           ! scale radius
  ! Black hole properties
  namelist /black_hole/ black_hole_mass_fraction  ! Black hole mass as fraction of
                                           ! halo virial mass
  ! Kinematics 
  namelist /kinematics/ velocity_anisotropy_amplitude,& ! Velocity anisotropy
       & velocity_anisotropy_radius        ! Anisotropy radius
  
  open(1,file=infile,status='old')
  rewind(1)
  read(1,nml=halo)
  rewind(1)
#ifdef STELLAR_DISC_EXPONENTIAL
  read(1,nml=disc)
  rewind(1)
#endif
#ifdef HERNQUIST_BULGE
  read(1,nml=bulge)
  rewind(1)
#endif
  read(1,nml=black_hole)
  rewind(1)
  read(1,nml=kinematics)
  rewind(1)    
  return
  
end subroutine set_parameters

