module kinematics
  use nrutils_modules 
  use constants
  use structure
  use variables
  use integrals
contains
  subroutine compute_sigma(nrbin,rtab,m200,c200,mbulge,fbulge,mdisc,fdisc,mbh,fmax,fcut,&
       & sigma_r,sigma_p,sigma_m)
    
    implicit none
    
    logical :: isbh                         ! Is a black hole present?
    
    real(kind=8) :: fmax,rmax,fcut,rcut
    integer(kind=4) :: i,j,k,l,m,n,nn
    real(kind=8) :: sum
    
    isbh=.false.
    
#if defined(HERNQUIST_HALO) && defined(NFW_HALO)
    write(*,*) 'Error: multiple DM halos!'
#endif      
    
    r200=getrvir(m200,deltavir,rhocrit0)
    rs=r200/c200                                  ! NFW Scale radius
    
    ! Compute component masses, in units of 1e10 solar masses
    
    if(isbh.eqv..true.) write(*,*) 'Assuming a central black hole...'
    
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
    mbulge=mbulge*m200          ! Bulge mass, as fraction of M200
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
    
    write(*,*) 'Computed radial velocity dispersion...'
    
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
    
    write(*,*) 'Computed projected velocity dispersion...'  
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
  
end module kinematics



