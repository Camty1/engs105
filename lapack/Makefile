FC=gfortran
FLAGS=-Wall -pedantic
OBJS=test
LIBS=-llapack -lblas

all: clean test

test:
	$(FC) $(FLAGS) test.f $(LIBS) -o $@
	
clean:
	rm -rf *.o *.mod test
