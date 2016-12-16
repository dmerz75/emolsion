# Makefile

# Declaration of variables
CC       := gcc
CXX      := g++
CL       := clang --analyze
# LINKER := gcc -fPIC
# OBJDIR := obj

# Compiler Flags: Use $(CF) for generic/old architectures
CF       := -g -std=c++11
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
EXEC     := emol
EXEF     := $(wildcard /usr/local/bin/$(EXEC)*)



# ---------------------------------------------------------------------
# Macros
MACRO   = -D
DCD     = -DDCDREAD
DCDW    = -DDCDREAD -DDCD_WRITE_B -DDCD_WRITE -DDCD_WRITE_E
CONS    = -DGET_CONTACTS
CON_BDA = -DCONTACTS_BEFORE -DCONTACTS_DURING -DCONTACTS_AFTER
MT      = -DMTMAP
MT2     = -DMTMAP_PRE -DMTMAP2
# Macros: Analysis Before. During. After.



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
dcd0:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DCD) -o test/$(EXEC)_dcdreader
	cd test && ./$(EXEC)_def mt.ref.pdb mt_partial.dcd 0 100 1
dcdr:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DCD) -o test/$(EXEC)_dcdreader
	cd test && ./$(EXEC)_dcdreader mt.ref.pdb mt_partial.dcd 6 27 3 # 6-9 .. 21-24-27.
dcdw:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DCDW) -o test/$(EXEC)_dcdwriter
	cd test && ./$(EXEC)_dcdwriter mt.ref.pdb mt_partial.dcd 6 27 3 # 6-9 .. 21-24-27.
contacts0:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(CONS) -o test/$(EXEC)_contacts
	cd test && ./$(EXEC)_contacts mt.ref.pdb mt_partial.dcd 6 27 3 # 6-9 .. 21-24-27.
contacts1:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DCD) $(CON_BDA) -o test/$(EXEC)_contactsbda
	cd test && ./$(EXEC)_contactsbda mt.ref.pdb mt_partial.dcd 6 27 3 # 6-9 .. 21-24-27.
contacts2:
# $(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(MT) $(DCD) $(CON_BDA) -DNDEBUG -o test/$(EXEC)_mtcontacts
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(MT) $(DCD) $(CON_BDA) -o test/$(EXEC)_mtcontacts
	cd test && ./$(EXEC)_mtcontacts mt.ref.pdb mt_partial.dcd 6 30 6 # 6-9 .. 21-24-27.
contacts3:
# $(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(MT2) $(DCD) $(CON_BDA) -DNDEBUG -o test/$(EXEC)_mtcontacts2
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(MT2) $(DCD) $(CON_BDA) -o test/$(EXEC)_mtcontacts2
	cd test && ./$(EXEC)_mtcontacts2 mt.ref.pdb mt_partial.dcd 6 30 6 # 6-9 .. 21-24-27.


# -----------------------------------------------------------------------------
# Make all.
all: \
	main
