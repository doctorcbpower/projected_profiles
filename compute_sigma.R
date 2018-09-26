compute_sigma <- function(ntab,rtab,m200=100,c200=10,
                          deltavir=200.,halotype=0,
                          mbulge=0.05,fbulge=0.01,mdisc=0,fdisc=0,
                          mbh=0.e-4,fmax=1.e2,fcut=1.5,beta=0.0,r_a=0.0) {
      if(is.null(ntab))
        ntab=500
      if(is.null(rtab)) { 
        lr=seq(-3,3,6/ntab)
        rtab=10**lr
      }
      
      dyn.load("~/Projected_Profiles/compute_sigma.so")
      out<-.Fortran("compute_sigma",
                    nr_tab=as.integer(ntab),
                    rtab=as.double(rtab),
                    m200=as.double(m200),
                    r200=as.double(1),
                    c200=as.double(c200),
                    deltavir=as.double(deltavir),
                    halotype=as.integer(halotype),
                    mbulge=as.double(mbulge),
                    fbulge=as.double(fbulge),
                    mdisc=as.double(mdisc),
                    fdisc=as.double(fdisc),
                    mbh=as.double(mbh),
                    fmax=as.double(fmax),
                    fcut=as.double(fcut),
                    beta=as.double(beta),
                    r_a=as.double(r_a),
                    sigr=as.double(rep(0,ntab)),
                    sigp=as.double(rep(0,ntab)),
                    sigm=as.double(rep(0,ntab))
                  )
     return(list(nrtab=out$nr_tab,
                 r200=out$r200,
                 m200=out$m200,
                 rtab=out$rtab[1:out$nr_tab],
                 sigma_r=out$sigr[1:out$nr_tab],
                 sigma_p=out$sigp[1:out$nr_tab],
                 sigma_m=out$sigm[1:out$nr_tab]))
}

