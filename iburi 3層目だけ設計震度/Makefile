SPH = sph3d.exe
OBJS =  mod_utils.o mod_const.o mod_simpara.o mod_type.o mod_vtk.o mod_kernel.o mod_box.o \
	mod_extforce.o mod_time_integral.o mod_density.o mod_xsph.o mod_particle.o sph3d.o

ifeq ($(FC),pgfortran)
	FFLAGS = -byteswapio -mp=allcores -Mbounds -g
else ifeq ($(FC),ifort)
	FFLAGS =  -convert big_endian -g
#FFLAGS =  -O2 -convert big_endian -openmp -p #for gprof
else
	FC = gfortran
#	FFLAGS =  -O2 -fconvert=big-endian -finline-functions -fopenmp
#	FFLAGS =  -O2 -fconvert=big-endian -finline-functions -fopenmp -DPAIR_LIST -DNO_UPDATE_POS -DNO_GRAVITY
	FFLAGS =  -O2 -fconvert=big-endian -finline-functions -fopenmp 
#	FFLAGS =  -O2 -fconvert=big-endian -finline-functions -fopenmp -DPAIR_LIST -DNO_UPDATE_POS
#	FFLAGS =  -O2 -fconvert=big-endian -finline-functions -fopenmp -DNO_UPDATE_POS -DNO_GRAVITY
#FFLAGS =  -O2 -fconvert=big-endian -finline-functions -fopenmp -DXSPH
endif

.SUFFIXES : .F90

ALL: $(SPH)

$(SPH): $(OBJS)
	$(FC) -o $@ $(OBJS) $(FFLAGS) $(TESTS) $(OPTN)

%.o: %.F90
	$(FC) -c $< $(FFLAGS) $(TESTS) $(OPTN)

%.mod: %.F90 %.o
	@:

.PHONY: clean
clean:
	rm -f *.o sph3d.exe *.mod *~
