CXX = g++
VERSION = 0.01
NAME = gspan
CXXFLAGS = -std=c++11 -O3 -Wall -Wno-deprecated
TARGETS = gspan
OBJ = gspan.o  main.o

all: gspan

gspan: ${OBJ}
	${CXX} ${LDFLAGS} -o ${TARGETS} ${OBJ}  ${LDFLAGS}

clean:
	rm -f ${OBJ} ${TARGETS} *~

