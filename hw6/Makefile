FC = gfortran
OBJs = hw6.f dsolve.o fortran_utils.o
FLAGS = -ffree-line-length-none -ffixed-line-length-none -llapack -lblas

out: $(OBJs)
	$(FC) -o $@ $^ $(FLAGS)

%.o: %.f
	$(FC) -c $< $(FLAGS)

clean:
	rm -rf *.o *.mod out
