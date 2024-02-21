CPP = g++
INCLUDES = -I/usr/include/eigen3
LINKS = -llapacke -llapack -lcblas -lblas
OPTIONS = -std=c++23

utils.o: utils.h utils.cpp
	$(CPP) $(OPTIONS) $(INCLUDES) -c utils.cpp -o utils.o 

qz.o: qz.h qz.cpp
	$(CPP) $(OPTIONS) $(INCLUDES) -c qz.cpp -o qz.o 

gensys.o: gensys.h gensys.cpp qz.o
	$(CPP) $(OPTIONS) $(INCLUDES) -c gensys.cpp -o gensys.o

run_eigen.out: run_eigen.cpp utils.o qz.o gensys.o
	$(CPP) $(OPTIONS) $(INCLUDES) run_eigen.cpp utils.o qz.o gensys.o -o run_eigen.out $(LINKS)
	chmod a+x ./run_eigen.out

setup_gensysR: gensys.h gensys.cpp qz.h qz.cpp gensysR.cpp gensysR.R
	rm ./gensysR/src/gensys.h ./gensysR/src/gensys.cpp ./gensysR/src/qz.h ./gensysR/src/qz.cpp ./gensysR/src/gensysR.cpp ./gensysR/R/gensysR.R
	cp ./gensys.h ./gensys.cpp ./qz.h ./qz.cpp ./gensysR.cpp ./gensysR/src/
	cp ./gensysR.R ./gensysR/R/

clean:
	rm -f *.o *.out *.a

all:
	rm -f *.o *.out *.a
	$(CPP) $(OPTIONS) $(INCLUDES) -c utils.cpp -o utils.o 
	$(CPP) $(OPTIONS) $(INCLUDES) -c qz.cpp -o qz.o 
	$(CPP) $(OPTIONS) $(INCLUDES) -c gensys.cpp -o gensys.o
	$(CPP) $(OPTIONS) $(INCLUDES) run_eigen.cpp utils.o qz.o -o run_eigen.out $(LINKS)
	chmod a+x ./run_eigen.out

