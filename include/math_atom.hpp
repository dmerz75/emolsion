# // 1: c++, c: atom_math
# // 2: name:
# // .atom_math
#ifndef _atom_math_
#define _atom_math_

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
/* #include "system.hpp" */
// #include "atom.hpp"
#include "system.hpp"


/* ---------------------------------------------------------
   Definitions:
   --------------------------------------------------------- */
/* #define BUFFERSIZE 900 */



/* ---------------------------------------------------------
   Classes:
   --------------------------------------------------------- */
// header_class



/* ---------------------------------------------------------
   function declarations:
   --------------------------------------------------------- */
// void ReadPDBfile(PDBfile *pdbfile,char filename[40]);
double matom_distance(Atom a, Atom b);

#endif
