FC = gfortran
OBJs = solve.o hw1.o

hw1: $(OBJs)
	$(FC) -o $@ $^

%.o: %.f
	$(FC) -c $<

clean:
	rm -rf *.o *.mod hw1
