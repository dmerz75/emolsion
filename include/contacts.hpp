# // 1: c++, c: contacts
# // 2: name:
# // .contacts
#ifndef __contacts_hpp_
#define __contacts_hpp_

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
#include "system.hpp"
#include "dcd.h"
// #include "dcdio.h" // inside of dcd.h

/* ---------------------------------------------------------
   headers:
   --------------------------------------------------------- */
#include "debug.h"
/* #include "system.hpp" */
// #include "atom.hpp"



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
// void get_contacts(Atom *a1,char *argv,int num_atoms);
void get_contacts(Atom *a1,Atom *a2,char dcdfilename[40],int num_atoms);
// void get_map_of_mtneighbors(std::vector<Atom*> chain_ref,std::vector<std::vector<int>> matrix);
void get_map_of_mtneighbors(std::vector<std::vector <Atom>> chain_ref,std::vector<std::vector<int>> matrix,
                            std::vector<std::pair<int,int>> dimers);

#endif
