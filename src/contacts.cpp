// contacts.cpp

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
#include "contacts.hpp"
// #include "dcd.h"
// #include "dcdio.h"

/* ---------------------------------------------------------
   functions
   --------------------------------------------------------- */
// void ReadPDBfile (PDBfile *pdbfile,char filename[40],System sys)
void get_contacts(Atom *a1,Atom *a2,char dcdfilename[40],int natoms)
{
    debug("welcome to get_contacts. 2 (inter)\n");

    std::cout << "DCD: " << dcdfilename << endl;

    // 2. to read a dcd.
    // natoms = ;
    molfile_timestep_t timestep;
    void *v;
    dcdhandle *dcd;
    // int natoms; // from the opening the dcd.
    float sizeMB =0.0, totalMB = 0.0;
    double starttime, endtime, totaltime = 0.0;

    v = open_dcd_read(dcdfilename,"dcd",&natoms);
    if (!v)
    {
        fprintf(stderr, "main) open_dcd_read failed for file %s\n", dcdfilename);
        exit(1);
    }

    dcd = (dcdhandle *)v;
    sizeMB = ((natoms * 3.0) * dcd->nsets * 4.0) / (1024.0 * 1024.0);
    totalMB += sizeMB;
    printf("main) file: %s\n", dcdfilename);
    printf("  %d atoms, %d frames, size: %6.1fMB\n", natoms, dcd->nsets, sizeMB);
    timestep.coords = (float *)malloc(3*sizeof(float)*natoms);

    for (int i=0; i<dcd->nsets; i++)
    {
        int rc = read_next_timestep(v, natoms, &timestep);
        if (rc)
        {
            fprintf(stderr, "error in read_next_timestep on frame %d\n", i);
            exit(1);
        }
       // std::cout << i << std::endl;
    }

    close_file_read(v);
    printf("Overall Size: %6.1f MB\n", totalMB);
}

void get_contacts(Atom *a1,char *argv)
{
    debug("welcome to get_contacts. 1 (intra)\n");


}
