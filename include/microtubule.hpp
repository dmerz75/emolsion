// 1: c++, c: microtubule
// 2: name: microtubule.hpp
// .microtubule
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
// #include <vector>
// #include <iterator> // istream_iterator
#include <map>

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
// indices:0  alpha/beta/not (0,1,-1)
//         1  0-chain,
//         2  (1,2)N,(3,4)M,(5,6)C

typedef std::map<std::string,int> MtIndexMapEntry;
typedef std::vector<MtIndexMapEntry> MtIndexMap;

/* ---------------------------------------------------------
   Classes:
   --------------------------------------------------------- */
// header_class


/* ---------------------------------------------------------
   function declarations:
   --------------------------------------------------------- */
// void ReadPDBfile(PDBfile *pdbfile,char filename[40]);


#endif
