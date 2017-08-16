module kinematics
  use nrutils_modules 
  use constants
  use structure
  use variables
contains
  subroutine compute_sigma(nr_tab,r_tab,&
       & virial_mass,nfw_concentration,virial_overdensity,&
       & bulge_mass,bulge_radius,disc_mass,disc_radius,black_hole_mass,&
       & fmax,fcut,velocity_anisotropy_amplitude,velocity_anisotropy_radius,&
       & sigr_tab,sigp_tab,sigm_tab)
    
    implicit none

    integer(kind=4) :: nr_tab
    real(kind=8), dimension(nr_tab), intent(inout) :: r_tab
    real(kind=8), dimension(nr_tab), intent(inout) :: sigr_tab
    real(kind=8), dimension(nr_tab), intent(inout) :: sigp_tab
    real(kind=8), dimension(nr_tab), intent(inout) :: sigm_tab    
    
    real(kind=8) :: virial_mass,virial_velocity,nfw_concentration,virial_overdensity
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

    m200=virial_mass
    mbulge=bulge_mass*m200
    mdisc=disc_mass*m200    
    mbh=black_hole_mass*m200

    c200=nfw_concentration
    r200=getrvir(m200,virial_overdensity,rhocrit0)
    rs=r200/c200                                  ! NFW Scale radius    
    fbulge=bulge_radius
    fdisc=disc_radius    

    isbh=.false.    
    if(mbh.gt.0) isbh=.true.
    
    ! Compute component masses, in units of 1e10 solar masses
    
    if(isbh.eqv..true.) then
       write(*,*) 'Assuming a central black hole...'
       write(*,*) 'MBH: ',mbh
       write(*,*)
    end if
    
#ifdef HERNQUIST_HALO
    write(*,*) 'Assuming a Hernquist halo...'
    write(*,*) 'Halo properties:'
    write(*,*) 'M200: ',m200
    write(*,*) 'R200: ',r200      
    write(*,*) 'Mtot: ',m200*(1+(1/c200)*sqrt(2*(dlog(1+c200)-c200/(1+c200))))
    write(*,*) 'Scale radius:',rs*sqrt(2*(dlog(1+c200)-c200/(1+c200)))
    write(*,*)
#endif
    
#ifdef NFW_HALO
    write(*,*) 'Assuming a NFW halo...'
    write(*,*) 'Halo properties:'
    write(*,*) 'M200: ',m200
    write(*,*) 'R200: ',r200      
    write(*,*) 'Scale radius:',rs
    write(*,*)
#endif            
    
#ifdef HERNQUIST_BULGE
    rbulge=fbulge*rs
#ifdef HERNQUIST_HALO  
    rbulge=rbulge*sqrt(2*(dlog(1+c200)-c200/(1+c200)))     ! in units of kiloparsecs      
#endif
    write(*,*) 'Assuming a Hernquist bulge...'
    write(*,*) 'Bulge properties:'
    write(*,*) 'Mtot: ',mbulge
    write(*,*) 'Scale radius:',rbulge
    write(*,*)
#endif  
    
    mhalo=m200-mbulge-mdisc-mbh
    rmax=fmax*r200
    rcut=fcut*r200

    do i=1,nrbin
       call qromo(sgfunc,dlog(rtab(i)),dlog(rmax),sum,midpnt)
       sigma_r(i)=-rtab(i)**(-2*get_vel_anisotropy(rtab(i)))*sum/rho(rtab(i))
    end do

#ifdef VERBOSE    
    write(*,*) 'Computed radial velocity dispersion...'
#endif    
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

#ifdef VERBOSE    
    write(*,*) 'Computed projected velocity dispersion...'  
#endif

    nr_tab=j

    do i=1,nr_tab
       sigr_tab(i)=sigma_r(i)
       sigp_tab(i)=sigma_p(i)
       sigm_tab(i)=sigma_m(i)
    end do

    deallocate(rtab,sigma_r,sigma_p,sigma_m)

    return
    
  end subroutine compute_sigma

  real(kind=8) function get_vel_anisotropy(r)
    
    implicit none
    
    real(kind=8) :: r
    
    if(r_a.le.0.0) then
       get_vel_anisotropy = beta0;  
    else
       get_vel_anisotropy = beta0*r_a*r_a/(r_a+r)**2
    end if
    
    return
    
  end function get_vel_anisotropy

  real(kind=8) function sgfunc(logr)
    implicit none
    real(kind=8) :: logr,r
    
    r=dexp(logr)
    
    sgfunc=-r*r**(2.*get_vel_anisotropy(r))*rho(r)*dphidr(r)
    return
  end function sgfunc
  
  real(kind=8) function funk(logr)
    implicit none
    
    real(kind=8) :: logr,r,rmin,rmax,sigr,fmin,fmax,fmid
    integer(kind=4) :: i,imin,imax
    
    r=exp(logr)
    
    imin=1
    imax=nrbin
    
    if((rtab(imin)-r)*(rtab(imax)-r).gt.0) stop
    
    do
       i=(imin+imax)/2
       if(r.lt.rtab(i)) imax=i
       if(r.gt.rtab(i)) imin=i
       if((imax-imin).le.1) exit
    end do
    
    sigr=(r-rtab(i-1))/(rtab(i)-rtab(i-1))*sigma_r(i)+&
         & (rtab(i)-r)/(rtab(i)-rtab(i-1))*sigma_r(i-1)
    
    funk = (1.0-get_vel_anisotropy(r)*radius*radius/r/r)*rho(r)*sigr*r*r/sqrt(r*r-radius*radius)
    
    return
    
  end function funk
   
end module kinematics



