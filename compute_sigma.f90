subroutine compute_sigma(nr_tab,r_tab,virial_mass,virial_radius,&
     & nfw_concentration,virial_overdensity,halo_type,&
     & bulge_mass,bulge_radius,disc_mass,disc_radius,black_hole_mass,&
     & fmax,fcut,velocity_anisotropy_amplitude,velocity_anisotropy_radius,&
     & sigr_tab,sigp_tab,sigm_tab)

  use kinematics
  
  implicit none
  
  integer(kind=4) :: nr_tab
  real(kind=8), dimension(nr_tab), intent(inout) :: r_tab
  real(kind=8), dimension(nr_tab), intent(inout) :: sigr_tab
  real(kind=8), dimension(nr_tab), intent(inout) :: sigp_tab
  real(kind=8), dimension(nr_tab), intent(inout) :: sigm_tab    
  
  real(kind=8) :: virial_mass,virial_radius,virial_velocity
  real(kind=8) :: nfw_concentration,virial_overdensity
  integer(kind=4) :: halo_type
  real(kind=8) :: disc_mass,disc_radius
  real(kind=8) :: bulge_mass,bulge_radius
  real(kind=8) :: black_hole_mass
  real(kind=8) :: velocity_anisotropy_amplitude,velocity_anisotropy_radius
  real(kind=8) :: fmax,rmax,fcut,rcut
  
  logical :: isbh                         ! Is a black hole present?
  
  integer(kind=4) :: i,j,k,l,m,n,nn
  real(kind=8) :: sum
  
  nrbin=nr_tab
  
  allocate(rtab(nrbin))       ! Tabulated radius
  allocate(sigma_r(nrbin))      ! Tabulated mass
  allocate(sigma_p(nrbin))      ! Tabulated mass  
  allocate(sigma_m(nrbin))      ! Tabulated mass
  
  do n=1,nrbin
     rtab(n)=r_tab(n)
     sigma_r(n)=0.0
     sigma_p(n)=0.0
     sigma_m(n)=0.0
  end do

  halotype=halo_type
  m200=virial_mass
  mbulge=bulge_mass*m200
  mdisc=disc_mass*m200    
  mbh=black_hole_mass*m200
  
  c200=nfw_concentration
  virial_radius=getrvir(m200,virial_overdensity,rhocrit0)
  r200=virial_radius
  rs=r200/c200                                  ! NFW Scale radius    
  fbulge=bulge_radius
  fdisc=disc_radius    

  if(mbulge.gt.0.0 .and. fbulge.le.0.0) stop 'Error: fbulge.le.0.0'
  if(mdisc.gt.0.0 .and. fdisc.le.0.0) stop 'Error: fdisc.le.0.0'
  
  ! Compute component masses, in units of 1e10 solar masses
  
  if(mbh.gt.0.0) then
     write(*,*) 'Assuming a central black hole...'
     write(*,*) 'MBH: ',mbh
     write(*,*)
  end if
  
  if(halo_type.eq.0) then
     write(*,*) 'Assuming a Hernquist halo...'
     write(*,*) 'Halo properties:'
     write(*,*) 'M200: ',m200
     write(*,*) 'R200: ',r200      
     write(*,*) 'Mtot: ',m200*(1+(1/c200)*sqrt(2*(dlog(1+c200)-c200/(1+c200))))
     write(*,*) 'Scale radius:',rs*sqrt(2*(dlog(1+c200)-c200/(1+c200)))
     write(*,*)
  end if
  
  if(halo_type.eq.1) then
     write(*,*) 'Assuming a NFW halo...'
     write(*,*) 'Halo properties:'
     write(*,*) 'M200: ',m200
     write(*,*) 'R200: ',r200      
     write(*,*) 'Scale radius:',rs
     write(*,*)
  end if
  
  if(mbulge.gt.0.0) then
     rbulge=fbulge*rs
     if(halo_type.eq.0) &
          & rbulge=rbulge*sqrt(2*(dlog(1+c200)-c200/(1+c200)))     ! in units of kiloparsecs      
     write(*,*) 'Assuming a Hernquist bulge...'
     write(*,*) 'Bulge properties:'
     write(*,*) 'Mtot: ',mbulge
     write(*,*) 'Scale radius:',rbulge
     write(*,*)
  end if
  
  mhalo=m200-mbulge-mdisc-mbh
  rmax=fmax*r200
  rcut=fcut*r200
  
  do i=1,nrbin
     call qromo(sgfunc,dlog(rtab(i)),dlog(rmax),sum,midpnt)
     sigma_r(i)=-rtab(i)**(-2*get_vel_anisotropy(rtab(i)))*sum/rho(rtab(i))
  end do
  
  rmax=0.98*rtab(nrbin)
  
  i=1
  
  do
     if(rtab(i).gt.rcut) exit
     radius=rtab(i)
     call qromo(funk,dlog(rtab(i)),dlog(rmax),sum,midpnt)
     sigma_p(i)=sum
     i=i+1
  end do
  
  j=i-1
  
  do i=1,j
     radius=rtab(i)
     call qromo(sigma_mass,dlog(radius),dlog(rmax),sum,midpnt)
     sigma_m(i)=sum
     sigma_p(i)=2.*sigma_p(i)/sigma_m(i)
  end do
  
  nr_tab=j
  
  do i=1,nr_tab
     sigr_tab(i)=sqrt(sigma_r(i))
     sigp_tab(i)=sqrt(sigma_p(i))
     sigm_tab(i)=sigma_m(i)
  end do
  
  deallocate(rtab,sigma_r,sigma_p,sigma_m)
  
  return
  
end subroutine compute_sigma

  
