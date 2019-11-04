# fortran compiler
# FC90 = gfortran
FC90 = ifort

# OpenMP-compiled library (intel)
LIB = -L/Users/Michael/phreeqcRM-3.4-OpenMP-intel/lib -lphreeqcrm

# optimization level
# FLAGS = -O3

# for debugging with ifort
FLAGS = -O0 -g -traceback -check all
# generate warnings and direct to compile.log
WARN = -warn all > compile.log 2>&1

# add LIB to FLAGS
FLAGS += $(LIB)
# add WARN to FLAGS
FLAGS += $(WARN)

PT_hMetal_files = PT_hMetal.f90 mod_PT_hMetal.f90 RM_interface.F90 kdtree2.f90

executable = exe_hMetal_PT
binDir = ./
dirExecutable = $(addprefix $(binDir),$(executable))
all: $(dirExecutable)

$(dirExecutable): $(PT_hMetal_files)
	$(FC90) RM_interface.F90 -c $(FLAGS)
	$(FC90) kdtree2.f90 -c $(FLAGS)
	$(FC90) mod_PT_hMetal.f90 -c $(FLAGS)
	$(FC90) $(PT_hMetal_files) -o $@ $(FLAGS)

clean:
	rm -r exe_* *.mod *.o *.out *.err myrun* *.dSYM *.inp compile.log
