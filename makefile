#OPT += -DNFW_HALO
OPT += -DHERNQUIST_HALO
OPT += -DHERNQUIST_BULGE
#OPT += -DSTELLAR_DISC_EXPONENTIAL
#OPT += -DVERBOSE

FC = gfortran

SRCS = main.f90 read_params.f90
SOBJ = $(SRCS:.f90=.o)

FILE = constants.f90 variables.f90 structure.f90 nrutils_modules.f90 kinematics.f90
FOBJ = $(FILE:.f90=.o)

compute_sigma.exe: $(FOBJ) read_params.o main.o
	$(FC) $^ -o compute_sigma.exe 

%.o: %.f90
	$(FC) -cpp $(OPT) -c $<
	touch $*.o $*.mod

clean: 
	rm *.mod *.o compute_sigma.exe
	touch makefile
