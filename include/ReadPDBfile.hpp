# // 1: c++, c: ReadPDBfile
# // 2: name:
# // .ReadPDBfile
#ifndef _ReadPDBfile_hpp_
#define _ReadPDBfile_hpp_

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


/* ---------------------------------------------------------
   Definitions:
   --------------------------------------------------------- */


/* ---------------------------------------------------------
   Classes:
   --------------------------------------------------------- */


/* ---------------------------------------------------------
   function declarations:
   --------------------------------------------------------- */
int ReadPDBfile(char filename[40],int total_atoms,Atom *aa);
vAtoms ReadPDBfile(char filename[40]);
// typedef std::vector<Atom> Atoms;


#endif
