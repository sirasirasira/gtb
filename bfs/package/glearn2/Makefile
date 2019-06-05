CXX = g++
VERSION = 0.01
NAME = glearn
CXXFLAGS = -std=c++11 -O3 -Wall -Wno-deprecated
TARGETS = glearn
OBJ = utils.o glearn.o main.o 

all: glearn

glearn: ${OBJ}
	${CXX} ${LDFLAGS} -o ${TARGETS} ${OBJ}

clean:
	rm -f ${OBJ} ${TARGETS} *~

