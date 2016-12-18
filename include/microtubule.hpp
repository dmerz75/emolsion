# // 1: c++, c: microtubule
# // 2: name: microtubule.hpp
# // .microtubule
#ifndef _microtubule_hpp_
#define _microtubule_hpp_

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
// #include "debug.h"
/* #include "system.hpp" */
// #include "atom.hpp"



/* ---------------------------------------------------------
   Definitions:
   --------------------------------------------------------- */
/* #define BUFFERSIZE 900 */
typedef std::vector<std::pair<int,int>> DimerList;
typedef std::vector<std::vector<int>> MtNeighbors;

// mt_matrix(dimers.size(), std::vector<int>(8,-1));
// std::vector<std::vector<int>>::iterator itmap;
// std::vector<int>::iterator itmap_n;


/* ---------------------------------------------------------
   Classes:
   --------------------------------------------------------- */
// header_class



/* ---------------------------------------------------------
   function declarations:
   --------------------------------------------------------- */
// void ReadPDBfile(PDBfile *pdbfile,char filename[40]);


#endif
