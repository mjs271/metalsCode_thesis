# fortran compiler
# FC90 = gfortran
FC90 = ifort

# OpenMP-compiled library (intel)
LIB = -L/Users/Michael/phreeqcRM-3.4-OpenMP-intel/lib -lphreeqcrm

# optimization level
FLAGS = -O3

# for debugging with ifort
# FLAGS = -O0 -g -traceback -check all
# generate warnings and direct to compile.log
WARN = -warn all > compile.log 2>&1

# add LIB to FLAGS
FLAGS += $(LIB)
# add WARN to FLAGS
FLAGS += $(WARN)

FD_hMetal_files = FD_hMetal.f90 FD_hMetal_mod.f90 RM_interface.F90

all: exe_FD_hMetal

exe_FD_hMetal: $(FD_hMetal_files)
	$(FC90) RM_interface.F90 -c $(FLAGS)
	$(FC90) FD_hMetal_mod.f90 -c $(FLAGS)
	$(FC90) $(FD_hMetal_files) -o $@ $(FLAGS)

clean:
	rm exe_* *.mod *.o *.out *.err myrun* error.inp
