module kinematics
  use nrutils_modules 
  use constants
  use structure
  use variables
contains
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



