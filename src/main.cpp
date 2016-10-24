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
// #include "config.hpp"
#include "ConfigFile.h"
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


    // read the config file!


    ConfigFile cf("conf.sop");

    std::string name,foo,bah;
    double number,bar;

    name   = cf.Value("section_1","name");
    number   = cf.Value("section_1","number");

    std::cout << name << std::endl;
    std::cout << number << std::endl;


    // water = cf.Value("section_2","water");
    // four  = cf.Value("section_2","four");

    // std::cout << foo   << std::endl;
    // std::cout << water << std::endl;
    // std::cout << four  << std::endl;


    std::cout << "\nclosing stdin,stdout,stderr.\n";
    fclose(stdin);
    fclose(stdout);
    fclose(stderr);
    return 0;
}
