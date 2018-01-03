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
CF0      := -g -std=c++11 -O3
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


#  ---------------------------------------------------------  #
#  Macros                                                     #
#  Analysis Before. During. After.
#  ---------------------------------------------------------  #
MACRO   = -D
DCDTEST = -DDCDTEST
DCD     = -DDCDREAD
DCDW    = -DDCDREAD -DDCD_WRITE_B -DDCD_WRITE -DDCD_WRITE_E
MT2     = -DMTBUILDMAP -DMTMAP2_BEFORE -DMTMAP2_DURING -DMTMAP2_AFTER
PHIPSI  = -DPHIPSI_B -DPHIPSI_M -DPHIPSI_E
TOPOw   = -DTOPO -DTOPO_write -DTOPO_write_mt -DMTBUILDMAP -DMTMAP2_BEFORE
TOPObond= -DTOPO -DTOPO_write -DTOPO_write_bonds
TOPOwh7 = -DTOPO -DTOPO_write -DTOPO_write_hsp70
TOPOwh7n= -DTOPO -DTOPO_write -DTOPO_write_hsp70 -DHSP70_NUC
TOPOr   = -DTOPO -DTOPO_read
TOPOmt  = -DDCDREAD -DTOPO -DTOPO_read -DTOPO_mt_BEFORE \
	-DTOPO_mt_DURING -DTOPO_mt_AFTER -DMTBUILDMAP -DTOPO_mt_SORT
TOPOmt2 = -DDCDREAD -DTOPO -DTOPO_read -DTOPO_mt_BEFORE \
	-DMTBUILDMAP -DTOPO_mt_SORT -DMTMAP2_DURING -DMTMAP2_AFTER
# TOPOmt3: external only.
# TOPOmt3 = -DDCDREAD -DTOPO -DTOPO_read -DTOPO_mt_BEFORE \
# -DMTBUILDMAP -DTOPO_mt_SORT -DMTMAP2_DURING -DMTMAP2_AFTER -DTOPO_ext_only
TOPOmt3 = -DDCDREAD -DTOPO -DTOPO_read \
	-DMTBUILDMAP -DTOPO_mt_SORT -DMTMAP2_DURING -DMTMAP2_AFTER -DTOPO_ext_only
# PFBEND  = -DDCDREAD -DPFBEND_BEFORE -DPFBEND_DURING -DPFBEND_AFTER
PFBEND  = -DPFBEND_BEFORE -DPFBEND_DURING -DPFBEND_AFTER


#  ---------------------------------------------------------  #
#  Macros' Descriptions:                                      #
#  ---------------------------------------------------------  #
# MT2:
# Print Analysis of Contacts by Subdomain.
# output_global_contacts_by_subdomain(global_contacts);
# fp_contacts = fopen("emol_mtcontacts_by_subdomain.dat", "w+");
# fp_contacts3 = fopen("emol_mtcontacts_by_subdomain3.dat", "w+");
# fp_contacts3n = fopen("emol_mtcontacts_by_subdomain3n.dat", "w+");



#  ---------------------------------------------------------  #
#  Targets:                                                   #
#  ---------------------------------------------------------  #
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
test-all-atom-dcd-0:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) -o test/$(EXEC)_dcdread
	cd test && ./$(EXEC)_dcdread 2kho_implicit.pdb 2kho_implicit.dcd 10 100 2
test-all-atom-dcd-1: # passes: DCDTEST
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DCDTEST) -o test/$(EXEC)_dcdread
	cd test && ./$(EXEC)_dcdread 2kho_implicit.pdb 2kho_implicit.dcd 10 100 2
test-all-atom-dcd-2: # causes the 1 chain into vvpAtoms failure.
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DCD) -o test/$(EXEC)_dcdread
	cd test && ./$(EXEC)_dcdread 2kho_implicit.pdb 2kho_implicit.dcd 10 100 2
test-all-atom-dcd-3: # succeeds because > 1 chain into vvpAtoms.
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DCD) -o test/$(EXEC)_dcdread
	cd test && ./$(EXEC)_dcdread 4EZW.pdb nil.dcd 10 100 2
main1:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) -o test/$(EXEC)_def
	cd test && ./$(EXEC)_def 4EZW.pdb nil.dcd
main2:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) -o test/$(EXEC)_def
	cd test && ./$(EXEC)_def mt.ref.pdb mt_partial.dcd
dcd0:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DCD) -o test/$(EXEC)_dcdreader
	cd test && ./$(EXEC)_dcdreader mt.ref.pdb mt_partial.dcd 0 100 1
dcdr:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DCD) -o test/$(EXEC)_dcdreader
	cd test && ./$(EXEC)_dcdreader mt.ref.pdb mt_partial.dcd 6 27 3 # 6-9 .. 21-24-27.
dcdw:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DCDW) -o test/$(EXEC)_dcdwriter
	cd test && ./$(EXEC)_dcdwriter mt.ref.pdb mt_partial.dcd 6 27 3 # 6-9 .. 21-24-27.
mt3:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DCD) $(MT2) -o test/$(EXEC)_mtcontacts3
	cd test && ./$(EXEC)_mtcontacts3 mt.ref.pdb mt_partial.dcd 6 30 6 # 6-9 .. 21-24-27.
mt6:
# 6 dimer test system.
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DCD) $(MT2) -o test/$(EXEC)_mtcontacts3
	cd test && ./$(EXEC)_mtcontacts3 mtdimer6.pdb mtdimer6.dcd 0 16 1 # 6-9 .. 21-24-27.
topo-w:
# 6 dimer test system for topology writing.
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(TOPOw) -o test/$(EXEC)_topo_write
	cd test && ./$(EXEC)_topo_write mtdimer6.pdb nil.dcd 0 16 2 emol_topology.top # 6-9 .. 21-24-27.
topo-write-bonds:
# 6 dimer test system for topology writing.
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(TOPObond) -o test/$(EXEC)_topo_write_bonds
	cd test && ./$(EXEC)_topo_write_bonds atp4b393.pdb nil.dcd 0 16 2 emol_topology.top
topo-w-hsp70:
# 6 dimer test system for topology writing.
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(TOPOwh7) -o test/$(EXEC)_topo_writehsp70
# cd test && ./$(EXEC)_topo_writehsp70 final_state44_ca_renum.pdb nil.dcd 0 16 2 emol_topology.top
# cd test && ./$(EXEC)_topo_writehsp70 hsp704b3934eCA.pdb nil.dcd 0 16 2 emol_topology.top
# cd test && ./$(EXEC)_topo_writehsp70 fullhsp70.pdb nil.dcd 0 16 2 emol_topology.top
topo-w-hsp70-nuc-atp:
# 6 dimer test system for topology writing.
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(TOPOwh7n) -o test/$(EXEC)_topo_writehsp70nuc
	cd test && ./$(EXEC)_topo_writehsp70nuc atp4b393.pdb nil.dcd 0 16 2 emol_topology.top
topo-w-hsp70-nuc-adp:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(TOPOwh7n) -o test/$(EXEC)_topo_writehsp70nuc
	cd test && ./$(EXEC)_topo_writehsp70nuc adp4b393.pdb nil.dcd 0 16 2 emol_topology.top
# for the peptide-linker-58
# cd test && ./$(EXEC)_topo_writehsp70 fullhsp70_ca_renum.pdb nil.dcd 0 16 2 emol_topology.top
topo-r:
# 6 dimer test system for topology writing.
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(TOPOr) -o test/$(EXEC)_topo_read
	cd test && ./$(EXEC)_topo_read mtdimer6.pdb nil.dcd 0 16 2 emol_topology.top
topo-mt:
# 6 dimer test system for topology writing.
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(TOPOmt2) -o test/$(EXEC)_mtcontacts_topo
	cd test && ./$(EXEC)_mtcontacts_topo mt.ref.pdb mt_partial.dcd 7 29 2 MT_regular_example.top
topo-mt-ext:
# external only evaluated.
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(TOPOmt3) -o test/$(EXEC)_mtcontacts_topo_extonly
	cd test && ./$(EXEC)_mtcontacts_topo_extonly mt.ref.pdb mt_partial.dcd 7 29 2 MT_regular_example.top
topo-mtp:
# 6 dimer test system for topology writing.
	$(CXX) $(CPPFILES) $(CF) -fopenmp $(INC) $(LIB) $(TOPOmt2) -o test/$(EXEC)_topo_mtp
	cd test && ./$(EXEC)_topo_mtp mt.ref.pdb mt_partial.dcd 7 29 2 MT_regular_example.top
# cd test && ./$(EXEC)_topo_mt mt.ref.pdb mt_partial.dcd 7 29 2 MT_regular.top
# cd test && ./$(EXEC)_topo_mt mt.ref.pdb mt_partial.dcd 7 29 2 emol_topology_example.top
# cd test && ./$(EXEC)_topo_mt mt.ref.pdb mt_partial.dcd 7 29 2 emol_topology_mtref.top
mtcontactstest:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DCD) $(MT2) -o test/$(EXEC)_mtcontacts2
	cd test && ./$(EXEC)_mtcontacts2 mt_test1.pdb mt_test1.dcd 4 220 5
angles-all-atoms:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DCD) $(PHIPSI) -o test/$(EXEC)_phipsi_angles
	cd test && ./$(EXEC)_phipsi_angles 2kho_implicit.pdb 2kho_implicit.dcd 10 100 2
pfbend:
# protofilament bending
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(PFBEND) -o test/$(EXEC)_mtpfbend
	cd test && ./$(EXEC)_mtpfbend mt.ref.pdb mt_partial.dcd 6 32 2


# ----------------------------------------------------------- #
# Deployment:                                                 #
# ----------------------------------------------------------- #
mtcontacts:
	$(CXX) $(CPPFILES) $(CF0) $(INC) $(LIB) $(DCD) $(MT2) -DNDEBUG -o bin/$(EXEC)_mtcontacts3n


#  ---------------------------------------------------------  #
#  Release:                                                   #
#  ---------------------------------------------------------  #
test-topo-mt:
# 6 dimer test system for topology writing.
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(TOPOmt2) -o test/$(EXEC)_mtcontacts_topo

test-topo-mt-ext:
# external only evaluated.
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(TOPOmt3) -o test/$(EXEC)_mtcontacts_topo_extonly




# -----------------------------------------------------------------------------
# Make all.
# -----------------------------------------------------------------------------
all: \
# To add here: insert (-DNDEBUG), change to (bin/)
# mt6:
# $(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DCD) $(MT2) -DNDEBUG -o bin/$(EXEC)_mtcontacts3
# topo-mt
# $(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(TOPOmt2) -DNDEBUG -o bin/$(EXEC)_mtcontacts_topo
# topo-mt-ext: # external only evaluated.
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(TOPOmt3) -DNDEBUG -o bin/$(EXEC)_mtcontacts_topo_extonly
