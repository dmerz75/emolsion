// main.cpp

// headers C
extern "C" {
// your functions here for the header

}

// my headers
#include "debug.h"
// #include "readfile.h"
// #include "chain.h"
#include "md.h"
#include "ReadPDBfile.hpp"
#include "system.hpp"
// #include "config.hpp"
// #include "ConfigFile.h"
// #include "Chameleon.h"

// // #include "contacts.h"
// #include "topology.h"
// #include "topology_charmm.h"
// #include "curvature.h"
// #include "chi.h"
// #include "dcd.h"

// headers C++
// #include <stdlib.h>
// #include <stdio.h> // printf
// #include <string.h> // strcpy, memcpy
// #include <new> // delete
// #include <ctype.h> // getopt - stuff
// #include <unistd.h> // getopt - stuff
// #include <vector>
// #include <algorithm> // bool & sort
#include <iostream>
// #include <iomanip> // setw
// #include <fstream> //
// #include <map> // map


// boost
// #include "boost/multi_array.hpp"


// VMD
// #ifndef _DCDIO_H_
// #define _DCDIO_H_
// #include "dcdio.h"
// #endif

// namespaces

int main(int argc, char *argv[]) {

    debug(">> Welcome to SOPCC! (debug)\n");
    printf(">> Welcome to SOPCC!\n");

    if (argc != 3) {
        std::cout << "Usage: " << argv[0] \
                  << " <Filename-reference-PDB>" \
                  << " <Filename-timelater-DCD>"
                  << std::endl;
        exit(1);
    }


    // Procedure:
    // 0. Create System.
    //    ()
    // 1. Read the config file. <conf.sop>
    //    () read_config_file
    //    1.1 Alter the default parameters. <def_param.h>
    //    () alter_default_parameters
    // 2. Read the PDB. <2KHO.pdb>
    //    () read_pdb
    // 3. Open the DCD.
    //    3.1 If DCD exists, read coordinates, velocities.
    //    3.2 If DCD does not exist, run minimization.
    // 4. Run MD.
    //    4.1 Write DCD every so often.
    //    4.2 Include PDBs
    // 5. Write final PDB.
    // 6. Close, Free up memory.

    int num_atoms;
    num_atoms = -1;


    // 0. Create the System.
    // System sys;
    // sys.print_prop();


    // 1. read the config file!


    // 2. Read the pdb.
    Atom a1[0];
    std::cout << "Currently, there are " << num_atoms << " atoms." << endl;
    num_atoms = ReadPDBfile(argv[1],num_atoms,a1);
    // delete a1;

    std::cout << "Currently, there are " << num_atoms << " atoms." << endl;
    // try Allocation Failure catch here?
    Atom allatoms[num_atoms];
    num_atoms = ReadPDBfile(argv[1],num_atoms,allatoms);



    // for(int i=0; i<num_atoms; i++)
    // {
    //     allatoms[i].print_coords();
    //     printf("%s  ",allatoms[i].chain.c_str());
    //     printf("%d\n",allatoms[i].resid);
    // }



    // Selection: H, precheck!
    int total_H = 0;
    for(int i=0; i<num_atoms; i++)
    {
        if(allatoms[i].chain.compare("H") == 0)
        {
            total_H += 1;
        }
    }
    std::cout << "total_H: " << total_H << endl;

    // Selection: H
    // overloaded, with/without return selection
    // select(fromthese,parameter-chain-resid,idn-H,num_select);
    int num_select;
    num_select = -1;
    num_select = select(allatoms,"chain","H",num_select);

    // // pointer.
    // Atom *selectionH;
    // // memory allocation.
    // try {
    //     selectionH = new Atom;
    // } catch (std::bad_alloc xa) {
    //     std::cout << "Allocation Failure\n";
    //     return 1;
    // }
    // num_select = select(allatoms,selectionH,"chain","H",num_select);





    std::cout << "\nclosing stdin,stdout,stderr.\n";
    fclose(stdin);
    fclose(stdout);
    fclose(stderr);
    return 0;
}
