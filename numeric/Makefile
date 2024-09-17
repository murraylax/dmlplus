CPP = g++
INCLUDES = -I/usr/include/eigen3 -I/usr/include/gsl
LINKS = -llapacke -llapack -lcblas -lblas -lgsl -lgslcblas -lm
OPTIONS = -std=c++23

dml_multiroot.o: dml_multiroot.h dml_multiroot.cpp
	$(CPP) $(OPTIONS) $(INCLUDES) -c dml_multiroot.cpp -o dml_multiroot.o

test_multiroot.out: dml_multiroot.o test_multiroot.cpp
	$(CPP) $(OPTIONS) $(INCLUDES) test_multiroot.cpp dml_multiroot.o -o test_multiroot.out $(LINKS)
	chmod a+x ./test_multiroot.out

clean:
	rm -f *.o *.out *.a



