# // 1: c++, c: protofilaments
# // 2: name:
# // .protofilaments
#ifndef _protofilaments_
#define _protofilaments_

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
// #include "atom.hpp"
// #include "microtubule.hpp"



/* ---------------------------------------------------------
   Definitions:
   --------------------------------------------------------- */
/* #define BUFFERSIZE 900 */

typedef std::pair<int,int> Dimer;
typedef std::vector<Dimer> Protofilament;
typedef std::vector<Protofilament> Protofilaments;


/* ---------------------------------------------------------
   Classes:
   --------------------------------------------------------- */
// header_class



/* ---------------------------------------------------------
   function declarations:
   --------------------------------------------------------- */
// void ReadPDBfile(PDBfile *pdbfile,char filename[40]);
Protofilaments determine_num_protofilaments(vAtoms aa);

#endif
