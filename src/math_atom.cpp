// math_atom.cpp

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
#include "math_atom.hpp"


/* ---------------------------------------------------------
   functions
   --------------------------------------------------------- */
// void ReadPDBfile (PDBfile *pdbfile,char filename[40],System sys)
double matom_distance(Atom a, Atom b)
{
    // a.print_Coords();
    // b.print_Coords();

    /* Distance between 2 points! */
    double distance;
    double xdist, ydist, zdist;

    xdist = b.x - a.x;
    ydist = b.y - a.y;
    zdist = b.z - a.z;

    // xdist = v2.x - v1.x;
    // ydist = v2.y - v1.y;
    // zdist = v2.z - v1.z;

    distance = sqrt( xdist * xdist + ydist * ydist + zdist * zdist);
    return distance;

}
