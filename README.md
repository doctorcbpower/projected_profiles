To build, simply type

% make

The default is to have a Hernquist halo, bulge, and a central black hole.

To run, type

% ./compute_sigma.exe ./parameters.txt

The format of the parameter file is

&HALO
v200=2.e2    - circular velocity at R200
c200=10      - c200
/

&DISC
mdisc=0.035    - disc mass, as a fraction of the halo mass
fdisc=0.0885   - disc scale radius
mstar=0.9      - stellar fraction (currently unused)
/

&BULGE
mbulge=0.1     - bulge mass, as a fraction of the halo mass
fbulge=0.25    - bulge scale radius
/

&BLACK_HOLE
mbh=5.0e-2     - black hole mass, in units of 1e10 Msol
/

&KINEMATICS
beta0=0.0      - velocity anisotropy
r_a=0.0        - anisotropy radius, for the Osipkov-Merritt model

The ouutput is in ASCII format, in a file with the

    projected_props.XXX_halo.(hern_bulge.)(beta0/var).txt

where XXX can be hern or nfw and hern_bulge is added in if the code is compiled with the
-DHERNQUIST_BULGE flag.

The format of the ASCII file is
   m200,r200,v200
   mbulge,rbulge
   mbh
   beta0,r_a
   i=1,n
     rtab(i),sqrt(sigma_r(i)),sqrt(sigma_p(i)),sigma_m(i)

Radii are in kpc, velocities in km/s, surface mass density in 1e10 Msol/kpc^2.
