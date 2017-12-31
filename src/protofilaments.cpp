// protofilaments.cpp

/* ---------------------------------------------------------
   standard libraries
   --------------------------------------------------------- */
// #include <stdio.h>
// #include <stdlib.h> // strtod?, stod
// #include <assert.h>
// #include <string.h> // string.
// #include <cctype>
// #include <algorithm> // remove_if, count
// #include <iostream>
// #include <fstream>
// #include <ctime>
// #include <list>        // std::list
// #include <vector>
// #include <iterator> // istream_iterator


/* ---------------------------------------------------------
   headers
   --------------------------------------------------------- */
#include "debug.h"
#include "protofilaments.hpp"



/* ---------------------------------------------------------
   functions
   --------------------------------------------------------- */
// void ReadPDBfile (PDBfile *pdbfile,char filename[40],System sys)
// protofilament determine_protofilaments(allatoms_ref);
Protofilaments determine_num_protofilaments(vAtoms aa){

    Dimer dimer;
    Protofilament pf;
    Protofilaments lst_allpf;
    int a,b;
    // typedef std::pair<int,int> Dimer;
    // typedef std::vector<Dimer> Protofilament;
    // typedef std::vector<Protofilament> Protofilaments;


    std::cout << "Inside function: determine_protofilaments." << std::endl;


    // 13 protofilamants
    for (int i=0; i < 13; i++)
    {
        a = 2 * i;
        b = 2 * i + 1;
        // std::cout << "protofilament: " << i << std::endl;
        // std::cout << a << b << std::endl;
        dimer = std::make_pair(a,b);
        pf.push_back(dimer);
        lst_allpf.push_back(pf);
        pf.clear();
    }
    return lst_allpf;
}
