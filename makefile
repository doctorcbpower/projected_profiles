#OPT += -DNFW_HALO
OPT += -DHERNQUIST_HALO
#OPT += -DHERNQUIST_BULGE
#OPT += -DSTELLAR_DISC_EXPONENTIAL
#OPT += -DVERBOSE

FC = gfortran
RC = R CMD SHLIB

SRCS = main.f90 read_params.f90 compute_sigma.f90
SOBJ = $(SRCS:.f90=.o)

FILE = constants.f90 variables.f90 structure.f90 nrutils_modules.f90 kinematics.f90
FOBJ = $(FILE:.f90=.o)

compute_sigma.exe: $(FOBJ) read_params.o compute_sigma.o main.o
	$(FC) $^ -o compute_sigma.exe 

compute_sigma.so : $(FOBJ) compute_sigma.f90
	$(RC) $^ -o compute_sigma.so

%.o: %.f90
	$(FC) -cpp $(OPT) -c $<
	touch $*.o $*.mod

clean: 
	rm *.mod *.o compute_sigma.exe
	touch makefile
