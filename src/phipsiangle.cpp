// phipsiangle.cpp

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
// #include <iterator> // istream_iterator
// #include <cmath>
// #include <vector>
// #include "boost/tuple/tuple.hpp"


/* ---------------------------------------------------------
   headers
   --------------------------------------------------------- */
#include "debug.h"
#include "md.h"
#include "system.hpp"
// #include "microtubule.hpp"
// #include "dcd.h"
// #include "dcdio.h"

/* ---------------------------------------------------------
   functions
   --------------------------------------------------------- */
void compute_phipsi(Atoms bb)
{
    std::cout << "Computing phi / psi." << std::endl;

    // for(auto a: bb)
    // {
    //     std::cout << a.x << std::endl;
    // }


    for(int i=0; i<bb.size(); i+=4)
    {
        std::cout << i << std::endl;
        std::cout << "0: " << bb[i].x   << " " << bb[i].y   << " " << bb[i].z;
        std::cout << std::endl;
        std::cout << "1: " << bb[i+1].x << " " << bb[i+1].y << " " << bb[i+1].z;
        std::cout << std::endl;
        std::cout << "2: " << bb[i+2].x << " " << bb[i+2].y << " " << bb[i+2].z;
        std::cout << std::endl;
        std::cout << "3: " << bb[i+3].x << " " << bb[i+3].y << " " << bb[i+3].z;
        std::cout << std::endl;

    }


}
