.SUFFIXES: .f90

include ./make.sys

default: all

modm = nr_files.o control_mod.o
modn = mod_mpi.o 
mpiu = global_parameters.o utility_mod.o
dmft = SCFT_constants.o SCMFT_print.o fftw3_mod.o wall_density.o wfieldfill.o
core = initialize_mod.o MDE_mod.o SCFT_iterate.o W_update.o WALL_density.o SCFT_Density_FE.o SCFT_dump.o SCFT_minimal.o
main = SCFT_main.o

fftwlib=/home/tangjz/Soft/fftw3.3.3/lib
fftwinc=/home/tangjz/Soft/fftw3.3.3/include
fftwCFLAG= -lfftw3_mpi -lfftw3 
MKLFLAG= -mkl

objects = $(modm) $(modn) $(mpiu) $(dmft) $(core) $(main)

all: SCFT

SCFT: $(objects)
	$(LINKER) $(objects) -o SCFT $(LFLAGS) -I$(fftwinc) -L$(fftwlib) $(fftwCFLAG)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f SCFT

clean-dat:
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
