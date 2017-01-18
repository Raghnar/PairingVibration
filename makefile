#$preamble

# The compiler
FC = gfortran
# flags for debugging or for maximum performance, comment as necessary
#FCFLAGS = -g -fbounds-check -Wall -Wno-tabs #debug
FCFLAGS = -O3 -funroll-loops #production

# flags forall (e.g. look for system .mod files, required in gfortran)
#FCFLAGS += -I/usr/include

# libraries needed for linking, unused in the examples
#LDFLAGS = -li_need_this_lib

# List of executables to be built within the package
PROGRAMS = PairVibrations #Scattering-k #Scattering-r 

PairVibrations: PairVibrations.o

# ==== Not modify after this ==== #

# "make" builds all
all: $(PROGRAMS)

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@.x $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)

