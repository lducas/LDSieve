PREFLAG = -O3 -funroll-loops -std=c++11 -I${HOME}/sw/include 
POSTFLAG = -L${HOME}/sw/lib -lntl -lm -L/usr/local/lib/ 
OBJ = LinAlg.o Lattice.o LDGaussSieve.o ListDecode.o main.o


all: ldsieve

ldsieve: LinAlg.o Lattice.o LDGaussSieve.o ListDecode.o main.o
	g++  ${PREFLAG} ${OBJ} -o $@ ${POSTFLAG}

%.o: %.cpp LinAlg.h Lattice.h LDGaussSieve.h ListDecode.h Buckets.h
	g++  ${PREFLAG} -c $< -o $@ ${POSTFLAG}

clean:
	rm ${OBJ} ldsieve