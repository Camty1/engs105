FC = gfortran
OBJs = test.o

test: $(OBJs)
	$(FC) -o $@ $^

%.o: %.f
	$(FC) -c $<

clean:
	rm -rf *.o *.mod test
