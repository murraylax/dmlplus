CPP = g++
DMLPLUS_LOCATION = /home/murraylax/gitrepo/code/cpp/dmlplus
DMLPLUS_INCLUDES = -I$(DMLPLUS_LOCATION)/gensys -I$(DMLPLUS_LOCATION)/numeric
INCLUDES = -I/usr/include/eigen3 -I/usr/include/gsl $(DMLPLUS_INCLUDES)
LINKS = -llapacke -llapack -lcblas -lblas -lgsl -lgslcblas -lm
OPTIONS = -std=c++23

utils.o: utils/utils.h utils/utils.cpp
	$(CPP) $(OPTIONS) $(INCLUDES) -c utils/utils.cpp -o utils/utils.o 

qz.o: gensys/qz.h gensys/qz.cpp
	$(CPP) $(OPTIONS) $(INCLUDES) -c gensys/qz.cpp -o gensys/qz.o 

gensys.o: gensys/gensys.h gensys/gensys.cpp gensys/qz.o
	$(CPP) $(OPTIONS) $(INCLUDES) -c gensys/gensys.cpp -o gensys/gensys.o

run_eigen.out: gensys/run_eigen.cpp utils/utils.o gensys/qz.o gensys/gensys.o
	$(CPP) $(OPTIONS) $(INCLUDES) gensys/run_eigen.cpp utils/utils.o gensys/qz.o gensys/gensys.o -o gensys/run_eigen.out $(LINKS)
	chmod a+x gensys/run_eigen.out

setup_gensysR: gensys/gensys.h gensys/gensys.cpp gensys/qz.h gensys/qz.cpp gensys/gensysR.cpp gensys/gensysR.R 
	rm gensysgensysR/src/gensys.h gensys/gensysR/src/gensys.cpp gensys/gensysR/src/qz.h gensys/gensysR/src/qz.cpp gensys/gensysR/src/gensysR.cpp gensys/gensysR/R/gensysR.R
	cp gensys/gensys.h gensys/gensys.cpp gensys/qz.h gensys/qz.cpp gensys/gensysR.cpp gensys/gensysR/src/
	cp gensys/gensysR.R gensys/gensysR/R/

dml_multiroot.o: numeric/dml_multiroot.h numeric/dml_multiroot.cpp
	$(CPP) $(OPTIONS) $(INCLUDES) -c numeric/dml_multiroot.cpp -o numeric/dml_multiroot.o

test_multiroot.out: numeric/dml_multiroot.o numeric/test_multiroot.cpp
	$(CPP) $(OPTIONS) $(INCLUDES) numeric/test_multiroot.cpp numeric/dml_multiroot.o -o numeric/test_multiroot.out $(LINKS)
	chmod a+x numeric/test_multiroot.out

clean:
	rm -f utils/*.o utils/*.out utils/*.a
	rm -f numeric/*.o numeric/*.out numeric/*.a
	rm -f gensys/*.o gensys/*.out gensys/*.a

all:
	rm -f *.o *.out *.a
	$(CPP) $(OPTIONS) $(INCLUDES) -c utils/utils.cpp -o utils/utils.o 
	$(CPP) $(OPTIONS) $(INCLUDES) -c gensys/qz.cpp -o gensys/qz.o 
	$(CPP) $(OPTIONS) $(INCLUDES) -c gensys/gensys.cpp -o gensys/gensys.o
	$(CPP) $(OPTIONS) $(INCLUDES) gensys/run_eigen.cpp utils/utils.o gensys/qz.o -o gensys/run_eigen.out $(LINKS)
	chmod a+x gensys/run_eigen.out

