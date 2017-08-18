module structure
  use constants
  use variables
contains
  real(kind=8) function menc(r)
    
    implicit none
    
    real(kind=8) :: r,z  ! Project radius, height above midplane
    real(kind=8) :: fofc,phi
    
    real(kind=8) :: x,a,mtotal
    
    menc=0.0

    if(halotype.eq.0) then    
       a=rs*sqrt(2*(dlog(1+c200)-c200/(1+c200)))
       mtotal=mhalo*(1+a/r200)*(1+a/r200)
       
       x=r/a      
       menc=menc+m200*(1+a/r200)*(1+a/r200)*x*x/(1+x)/(1+x)
    end if
    
    if(halotype.eq.1) then
       x=r/rs
       
       fofc=dlog(1+c200)-c200/(1+c200)
       
       menc = menc+(mhalo/fofc)*(log(1+x)-x/(1+x))
    end if
      
    if(mbulge.gt.0.0) then
       a=fbulge*rs
       mtotal=mbulge
       
       x=r/a
       
       menc=menc+mbulge*x*x/(1+x)/(1+x)      
    end if

    if(mbh.gt.0.0) menc=menc+mbh
      
    if(mdisc.gt.0.0) then
       a=rdisc

       x=r/a

       menc=menc+mdisc*(1-exp(-x)*(1+x))
    end if

    return
    
  end function menc
  
  real(kind=8) function mhenc(r)
    implicit none
    
    real(kind=8) :: r,z  ! Project radius, height above midplane
    real(kind=8) :: fofc,phi
    
    real(kind=8) :: x,a,mtotal
    
    if(halotype.eq.0) then
       a=rs*sqrt(2*(dlog(1+c200)-c200/(1+c200)))
       mtotal=mhalo*(1+a/r200)*(1+a/r200)
       
       x=r/a      
       mhenc=mtotal*x*x/(1+x)/(1+x)
    end if
    
    if(halotype.eq.1) then
       x=r/rs
       
       fofc=dlog(1+c200)-c200/(1+c200)
       
       mhenc = (mhalo/fofc)*(log(1+x)-x/(1+x))
    end if

    return
    
  end function mhenc
    
  real(kind=8) function mbenc(r)
    implicit none
    
    real(kind=8) :: r,z  ! Project radius, height above midplane
    real(kind=8) :: fofc,phi
    
    real(kind=8) :: x,a,mtotal
    
    a=fbulge*rs
    mtotal=mbulge
    
    x=r/a
    
    mbenc=mbulge*x*x/(1+x)/(1+x)      

    return
    
  end function mbenc
  
  real(kind=8) function getrvir(mvir,dvir,rhoc)
    use constants
    implicit none
    real(kind=8) :: mvir,dvir,rhoc
    getrvir=(3*mvir/(4*pi) * (1/dvir) * (1/rhoc))**(1./3.)
    return
  end function getrvir
  
  real(kind=8) function getmvir(vvir,dvir,rhoc)
    use constants
    implicit none
    real(kind=8) :: vvir,dvir,rhoc
    getmvir=sqrt((3./4./pi) * (1/dvir) * (1/rhoc))
    getmvir=getmvir*(vvir/sqrt(ggdt))**3.0
    return
  end function getmvir
  
    
  real(kind=8) function rhodm(r)
    
    implicit none
    
    real(kind=8) :: a,x,r,rho_scale
    
    if(halotype.eq.0) then
       a=rs*sqrt(2*(dlog(1+c200)-c200/(1+c200)))
       
       x=r/a
       
       rho_scale = (mhalo/2/pi/a**3)*(r200/a)**2/(1+r200/a)**2
       
       rhodm = rho_scale/(x*(1+x)*(1+x)*(1+x))
    end if
    
    if(halotype.eq.1) then
       x=r/rs
       
       rho_scale=mhalo/(dlog(1+c200)-c200/(1+c200))
       
       rho_scale = rho_scale/(4.*pi*rs**3)
       
       rhodm = rho_scale/(x*(1+x)*(1+x))    
    end if
    
    return
    
  end function rhodm
  
  real(kind=8) function rhobulge(r)
    
    implicit none
    
    real(kind=8) :: a,x,r,rho_scale,mtotal

    rhobulge=0.0
    
    if(mbulge.le.0.0) return
       
    a=fbulge*rs
    mtotal=mbulge
    
    x=r/a      
    
    rho_scale = (mtotal/2/pi/a**3)
    
    rhobulge = rho_scale/(x*(1+x)*(1+x)*(1+x))
    
    return
    
  end function rhobulge
  
  real(kind=8) function rho(r)
    
    implicit none
    
    real(kind=8) :: r
    
    rho = rhodm(r)

    rho = rho + rhobulge(r)
      
    return
    
  end function rho
  
  real(kind=8) function dphidr(r)
    implicit none
    
    real(kind=8) :: r,a,mtotal
    
    if(halotype.eq.0) then
       a=rs*sqrt(2*(dlog(1+c200)-c200/(1+c200)))
       mtotal=m200*(1+a/r200)*(1+a/r200)-mbulge-mbh
       
       dphidr = ggdt*mtotal/(r+a)/(r+a)
    end if
    
    if(halotype.eq.1) then
       mtotal=m200-mbulge-mbh
       dphidr = (ggdt*mtotal/rs/rs/(dlog(1+c200)-c200/(1+c200))) * (rs/r)*(dlog(1+r/rs)/(r/rs)-1./(1+r/rs))
    end if
    
    if(mbulge.gt.0.0) then
       a=fbulge*rs
       
       dphidr = dphidr+ggdt*mbulge/(r+a)/(r+a)
    end if

    if(mbh.gt.0.0) dphidr = dphidr+ggdt*mbh/r
    
    return
    
  end function dphidr
  
  real(kind=8) function sigma_mass(logr)
    implicit none
    real(kind=8) :: logr,r,lognum
    
    r=dexp(logr)
    
    lognum=dlog(2.d0)+2.*dlog(r)+dlog(rho(r))-0.5*dlog(r*r-radius*radius)
    sigma_mass=dexp(lognum)
    return
  end function sigma_mass
  
end module structure
