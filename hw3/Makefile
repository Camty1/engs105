FC = gfortran
OBJs = hw3.o fortran_utils.o dsolve.o solve.o
FLAGs = -ffree-line-length-none -ffixed-line-length-none

out: $(OBJs)
	$(FC) -o $@ $^

%.o: %.f
	$(FC) -c $< $(FLAGs)

clean:
	rm -rf *.o *.mod out
