# GNUmakefile
# (c) R.Shirakawa 2019.1

COMPILER = g++
CFLAGS_BASE = -W -Wall -std=c++11 -MMD -MP -MF
CFLAGS = -O3 -DNDEBUG $(CFLAGS_BASE)
LDFLAGS = -pg
LIBS =

SRC_DIR = src
INC_DIR = include
OBJ_DIR = .obj
DPN_DIR = .depend
EXE_DIR = bin

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(addprefix $(OBJ_DIR)/, $(notdir $(SRCS:.cpp=.o)))
DPNS = $(addprefix $(DPN_DIR)/, $(notdir $(SRCS:.cpp=.d)))
TARGETS = $(EXE_DIR)/a.out

all: $(TARGETS)

$(TARGETS): $(OBJ_DIR) $(DPN_DIR) $(EXE_DIR) $(OBJS)
	$(COMPILER) $(LDFLAGS) -o $(TARGETS) $(OBJS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(COMPILER) $(CFLAGS) $(DPN_DIR)/$*.d -I$(INC_DIR) -o $@ -c $<

.PHONY: cln
cln: clean

.PHONY: clean
clean:
	rm -f $(OBJS) $(DPNS) $(TARGETS) *~

.PHONY: re
re: clean all

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(DPN_DIR):
	mkdir -p $(DPN_DIR)

$(EXE_DIR):
	mkdir -p $(EXE_DIR)

-include $(DPNS)
