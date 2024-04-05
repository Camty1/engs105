FC = gfortran
OBJs = dsolve.o hw2.o
FLAGs = -ffree-line-length-none -ffixed-line-length-none

hw2: $(OBJs)
	$(FC) -o $@ $^

%.o: %.f
	$(FC) -c $< $(FLAGs)

clean:
	rm -rf *.o *.mod hw2
