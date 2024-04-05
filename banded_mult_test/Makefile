FC = gfortran
OBJs = banded_mult.o test.o

test: $(OBJs)
	$(FC) -o $@ $^

%.o: %.f
	$(FC) -c $<

clean:
	rm -rf *.o *.mod test
