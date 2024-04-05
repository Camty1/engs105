FC = gfortran
FLAGS = -ffree-line-length-none -ffixed-line-length-none -llapack -lblas

UTILS = fortran_utils.f
OUTS = get_truth get_meas invert_data

# Default target
all: $(OUTS)

# Rule for compiling utility files into object files
utils.o: $(UTILS)
	$(FC) -c $(UTILS) -o utils.o $(FLAGS) 

# Rule for each program
get_truth: utils.o get_truth.f
	$(FC) -o get_truth get_truth.f utils.o $(FLAGS)

get_meas: utils.o get_meas.f
	$(FC) -o get_meas get_meas.f utils.o $(FLAGS)

invert_data: utils.o invert_data.f
	$(FC) -o invert_data invert_data.f utils.o $(FLAGS)

# Clean up
clean:
	rm -f *.o $(OUTS)

