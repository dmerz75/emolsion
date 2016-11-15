# Makefile

# Declaration of variables
CC       := gcc
CXX      := g++
CL       := clang --analyze
# LINKER := gcc -fPIC
# OBJDIR := obj

# Compiler Flags: Use $(CF) for generic/old architectures
CF       := -g
CC_FLAGS := -g -O3
CFLAGS   := -O2 -g -Wall
CFLAGS_1 := -ansi -std=gnu99
CFLAGS_2 := -ansi -pedantic -std=gnu99 -Wall -W
CFLAGS_2 := -ansi -pedantic
CFLAGS_3 := -Wconversion -Wshadow -Wcast-qual -Wwrite-strings -g -O2

# Valgrind
VAL      := valgrind --track-origins=yes -v
VALFULL  := valgrind --track-origins=yes --leak-check=full -v
VALMEM   := valgrind --track-origins=yes --tool=memcheck --leak-check=full --show-leak-kinds=all --show-reachable=yes --num-callers=20 --track-fds=yes -v
VALMASS  := valgrind --tool=massif prog

# c files
# SOURCES  := $(CFILES) # $(CFILES2) all CFILES
CPPFILES   := $(wildcard src/*.cpp)


# o files
OBJECTS  := $(SOURCES:.cpp=.o)
OBJECTS_A:= $(SOURCES:.cpp=_a.o)

# lib
LIB      := -pthread
# LIB    := -pthread -larmadillo
# LIB    := -pthread -lmongoclient -L lib -lboost_thread-mt -lboost_filesystem-mt -lboost_system-mt

# include
INC      := -Iinclude
INC_EIGEN:= -I/usr/include/eigen3


# DIRS
TOPDIR   := $(shell pwd)
SRCDIR   := src
BUILDDIR := build
BIN_DIR  := bin
BIN      := run_segment
SRCEXT   := cpp
ODIR     := build
TESTD    := test

# Executable:
EXEC     := run_readpdb
EXEF     := $(wildcard /usr/local/bin/$(EXEC)*)



# ---------------------------------------------------------------------
# Macros
MACRO = -D


# ---------------------------------------------------------------------
# Main target
$(EXEC): $(OBJECTS)
	$(CXX) $(OBJECTS) -lm -o $(EXEC)

# To remove generated files
clean:
	rm -f $(OBJECTS) $(OBJECTS_A)
# rm -f $(EXEC) $(OBJECTS)

go:
	./$(EXEC)

foo:
	@echo echoing foo
	@echo $(CPPFILES)
	@echo $(SOURCES)
	@echo $(OBJECTS)
	@echo $(OBJECTS_A)


# ---------------------------------------------------------------------
# Examples & Testing
# ---------------------------------------------------------------------
def: $(OBJECTS)
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DEFAULT) -o test/$(EXEC)_def
	cd test && $(EXEC)_def proto1.pdb proto2501.pdb 14 2
main:
# $(CXX) $(CPPFILES) $(CFLAGS_3) $(INC) $(LIB) -o test/$(EXEC)_def
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) -o test/$(EXEC)_def
	cd test && ./$(EXEC)_def 2KHO.pdb 2KHO.dcd
main0:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) -o test/$(EXEC)_def
	cd test && ./$(EXEC)_def 2KHO.pdb 2KHO.dcd
main1:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) -o test/$(EXEC)_def
	cd test && ./$(EXEC)_def 4EZW.pdb nil.dcd
main2:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) -o test/$(EXEC)_def
	cd test && ./$(EXEC)_def mt.ref.pdb mt_partial.dcd



# -----------------------------------------------------------------------------
# Make all.
all: \
	main
