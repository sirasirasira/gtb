CXX = g++
VERSION = 0.01
NAME = glearn
CXXFLAGS = -std=c++11 -O3 -Wall -Wno-deprecated
LDFLAGS = -lboost_system -lboost_timer
TARGETS = glearn
OBJ = utils.o glearn.o main.o 

all: glearn

glearn: ${OBJ}
	${CXX} ${LDFLAGS} -o ${TARGETS} ${OBJ}

clean:
	rm -f ${OBJ} ${TARGETS} *~

#g++ -pg -I../boost_local/include -L../boost_local/lib -lboost_system -lboost_timer *.cpp
