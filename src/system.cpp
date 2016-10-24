// system.cpp

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
#include "system.hpp"


/* ---------------------------------------------------------
   functions
   --------------------------------------------------------- */
// void ReadPDBfile (PDBfile *pdbfile,char filename[40],System sys)
int select(Atom aa,char param[],char criterion[],int num)
{


}
// int select(Atom aa,Atom asel,string param,string criterion,int num)
// {


// }


void get_minmax(System sys)
{
    debug("you made it!\n");

    // float minx, miny, minz;
    // float maxx, maxy, maxz;

    // minx = miny = minz = 9999.0;
    // maxx = maxy = maxz = -9999.0;


    // for(int i=0; i<sys.num_segment; i++)
    // {
    //     // printf("%d\n",i);

    //     for(int j=0; j<sys.segment[i].num_residue; j++)
    //     {
    //         // printf("%d ",j);

    //         for(int k=0; k<sys.segment[i].residue[j].num_atom; k++)
    //         {
    //             // sys.segment[i].residue[j].atom[k].print_info();
    //             // sys.segment[i].residue[j].atom[k].print_coords();

    //             if (sys.segment[i].residue[j].atom[k].x < minx)
    //             {
    //                 minx = sys.segment[i].residue[j].atom[k].x;
    //             }
    //             if (sys.segment[i].residue[j].atom[k].y < miny)
    //             {
    //                 miny = sys.segment[i].residue[j].atom[k].y;
    //             }
    //             if (sys.segment[i].residue[j].atom[k].z < minz)
    //             {
    //                 minz = sys.segment[i].residue[j].atom[k].z;
    //             }

    //             if (sys.segment[i].residue[j].atom[k].x > maxx)
    //             {
    //                 maxx = sys.segment[i].residue[j].atom[k].x;
    //             }
    //             if (sys.segment[i].residue[j].atom[k].y > maxy)
    //             {
    //                 maxy = sys.segment[i].residue[j].atom[k].y;
    //             }
    //             if (sys.segment[i].residue[j].atom[k].z > maxz)
    //             {
    //                 maxz = sys.segment[i].residue[j].atom[k].z;
    //             }



    //         }
    //         // printf("\n");
    //     }
    //     // printf("\n");
    // }

    // sys.minx = minx;
    // sys.miny = miny;
    // sys.minz = minz;

    // sys.maxx = maxx;
    // sys.maxy = maxy;
    // sys.maxz = maxz;
}
