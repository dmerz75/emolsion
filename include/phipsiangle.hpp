# // 1: c++, c: cpp
# // 2: name: phipsiangle
# // phipsiangle.cpp
#ifndef _phipsiangle_cpp_
#define _phipsiangle_cpp_

/* ---------------------------------------------------------
   libraries:
   --------------------------------------------------------- */
// #include <stdio.h>
// #include <stdlib.h> // strtod?, stod
// #include <assert.h>
// #include <string.h> // string
// #include <cctype.h>
/* #include <algorithm> // remove_if, count */
/* #include <iostream> */
/* #include <fstream> */
/* #include <ctime> */
// #include <list>        // std::list
/* #include <vector> */
// #include <iterator> // istream_iterator



/* ---------------------------------------------------------
   headers:
   --------------------------------------------------------- */
#include "debug.h"
#include "system.hpp"
/* #include "system.hpp" */
// #include "atom.hpp"



/* ---------------------------------------------------------
   Definitions:
   --------------------------------------------------------- */
/* #define BUFFERSIZE 900 */


typedef std::vector<std::pair<double,double>> PhiPsi;
typedef std::vector<PhiPsi> Global_PhiPsi;


/* ---------------------------------------------------------
   Classes:
   --------------------------------------------------------- */
// header_class



/* ---------------------------------------------------------
   function declarations:
   --------------------------------------------------------- */
// void ReadPDBfile(PDBfile *pdbfile,char filename[40]);
PhiPsi compute_phipsi(Dihedral dh);


#endif
