CXX = g++
VERSION = 0.01
CXXFLAGS = -O3 -Wall -Wno-deprecated -std=c++11
TARGETS = finder
OBJ = finder.o

all: ${TARGETS}

${TARGETS}: ${OBJ}
	${CXX} ${LDFLAGS} -o ${TARGETS} ${OBJ}

clean:
	rm -f ${OBJ} ${TARGETS} *~ *.out


