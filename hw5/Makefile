FC = gfortran
FLAGS = -ffree-line-length-none -ffixed-line-length-none -llapack -lblas

UTILS = fortran_utils.f dsolve.f
OUTS = steady_state transient

# Default target
all: $(OUTS)

# Rule for compiling utility files into object files
utils: $(UTILS)
	$(FC) -c $(UTILS) $(FLAGS) 

# Rule for each program
steady_state: dsolve.o fortran_utils.o steady_state.f
	$(FC) -o steady_state steady_state.f dsolve.o fortran_utils.o $(FLAGS)

transient: dsolve.o fortran_utils.o transient.f
	$(FC) -o transient transient.f dsolve.o fortran_utils.o $(FLAGS)

# Clean up
clean:
	rm -f *.o $(OUTS)

