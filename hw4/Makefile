FC = gfortran
OBJs = hw4.o dsolve.o fortran_utils.o
FLAGs = -ffree-line-length-none -ffixed-line-length-none

out: $(OBJs)
	$(FC) -o $@ $^

%.o: %.f
	$(FC) -c $< $(FLAGs)

clean:
	rm -rf *.o *.mod out
