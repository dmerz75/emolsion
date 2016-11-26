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


    /* ---------------------------------------------------------
       Begin DCD Open.
       --------------------------------------------------------- */
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
    timestep.coords = (float *)malloc(3*sizeof(float)*natoms);
    printf("main) file: %s\n", dcdfilename);
    printf("  %d atoms, %d frames, size: %6.1fMB\n", natoms, dcd->nsets, sizeMB);




    /* ---------------------------------------------------------
       The Default list of contacts.
       --------------------------------------------------------- */
    // By index?
    // a1 --> a2
    printf("the_size_of_a1: %d\n",a1[0].num_atoms);





    /* ---------------------------------------------------------
       Begin DCD Loop.
       --------------------------------------------------------- */
    for (int i=0; i<dcd->nsets; i++)
    {
        int rc = read_next_timestep(v, natoms, &timestep);
        if (rc)
        {
            fprintf(stderr, "error in read_next_timestep on frame %d\n", i);
            exit(1);
        }

        // DCD-COORDS:
        // void load_dcd_to_atoms(dcdhandle *dcd,Atom *aa);
        load_dcd_to_atoms(dcd,a1);
        load_dcd_to_atoms(dcd,a2);


        for(int j=0; j<a1[0].num_atoms; j++)
        {
            if (j <= 5)
            {
                printf("%f %f %f\n",a1[j].x,a1[j].y,a1[j].z);
            }
        }
        printf("---\n");

        for(int k=0; k<a2[0].num_atoms; k++)
        {
            if (k <= 5)
            {
                printf("%f %f %f\n",a1[k].x,a1[k].y,a1[k].z);
            }
        }


    }

    close_file_read(v);
    printf("Overall Size: %6.1f MB\n", totalMB);
}

void get_contacts(Atom *a1,char *argv)
{
    debug("welcome to get_contacts. 1 (intra)\n");


}


void get_map_of_mtneighbors(std::vector<Atom*> chain_ref,std::vector<std::vector<int>> matrix)
{
    // std::vector<Atom*> chain_ref;

    printf("Welcome to get_map_of_mtneighbors!\n");
    std::cout << matrix.size() << std::endl;

    // for(auto vec: matrix)
    // {
    //     // std::cout << vec[0] << vec[1] << vec[2] << vec[3] << std::endl;
    //     for(auto x: vec)
    //     {
    //         std::cout << x << std::endl;
    //     }
    //     std::cout << "\n" << std::endl;
    // }


    for(auto ch: chain_ref)
    {
        // if <Atom>
        // std::cout << ch.num_atoms << std::endl;
        // if <Atom*>
        std::cout << ch->num_atoms << std::endl;
        // std::cout << ch->size() << std::endl;

        for(int i=0; i<ch->num_atoms; i++)
        {
            // std::cout << ch[i].num_atoms << std::endl;
            std::cout << ch[i].x << std::endl;
        }

        // std::cout << ch.num_atoms << std::endl;

        // for(auto a: ch)
        // {
        //     std::cout << aa->num_atoms << std::endl;
        // }

        // std::cout << "\n" << std::endl;
    }


}
