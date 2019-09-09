CXX = g++
D_DIR = .depend
CXXFLAGS_BASE = -W -Wall -std=c++11 -MMD -MP -MF $(D_DIR)/$*.d
# -Wno-deprecated 
#CXXFLAGS = -O3 -DNDEBUG $(CXXFLAGS_BASE)
CXXFLAGS = -pg -g $(CXXFLAGS_BASE)
LDFLAGS = -pg
TARGETS = a.out
SRCS = $(wildcard *.cpp)
OBJDIR = .objs
OBJS = $(addprefix $(OBJDIR)/, $(SRCS:.cpp=.o))

all : $(TARGETS)

-include $(D_DIR)/*

$(TARGETS) : $(OBJDIR)  $(OBJS)
	${CXX} ${LDFLAGS} -o ${TARGETS} ${OBJS}

$(OBJDIR)/%.o : %.cpp
	${CXX} ${CXXFLAGS} -c -o $@ $<

.PHONY: cln
cln : clean
.PHONY: clean
clean:
	rm -f ${TARGETS} ${OBJS} $(D_DIR)/* *~

.PHONY: re
re : clean all

$(OBJDIR) :
	mkdir -p $(OBJDIR)
	mkdir -p $(D_DIR)
