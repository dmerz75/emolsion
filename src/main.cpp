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
    //    3.0 allocate for aa_later.
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
    std::cout << "1. Currently, there are " << num_atoms << " atoms." << endl;
    num_atoms = ReadPDBfile(argv[1],num_atoms,a1);
    // delete a1;

    std::cout << "2. Currently, there are " << num_atoms << " atoms." << endl;
    // try Allocation Failure catch here?

    // memory allocation.
    // #include <iostream>
    // Chain *chain_ref;
    // try {
        // chain_ref = new Chain [num_chains];
    Atom *aa_ref;
    try
    {
        aa_ref = new Atom[num_atoms];
    }
    catch (std::bad_alloc xa)
    {
        std::cout << "Allocation Failure\n";
        exit(1);
    }
    num_atoms = ReadPDBfile(argv[1],num_atoms,aa_ref);


    // Assign Total.
    for(int i=0; i<num_atoms; i++)
    {
        // aa_ref[i].print_coords();
        // printf("%s  ",aa_ref[i].chain.c_str());
        // printf("%d\n",aa_ref[i].resid);
        aa_ref[i].num_atoms = num_atoms;
    }
    // aa_ref[0].num_atoms = num_atoms;


    // Selection: H, precheck!
    int total_H = 0;
    for(int i=0; i<num_atoms; i++)
    {
        if(aa_ref[i].chain.compare("H") == 0)
        {
            total_H += 1;
        }
    }
    std::cout << "total_atoms_on_chain_H(4EZW-55?): " << total_H << endl;




    /* ---------------------------------------------------------
       Begin Selection:
       --------------------------------------------------------- */
    // overloaded, with/without return selection
    // select(fromthese,parameter-chain-resid,idn-H,num_select);
    int num_select;

    // CHAIN:
    // num_select = select(aa_ref,"chain","A",num_select);
    // Example 1.
    num_select = -1;
    num_select = system_select(aa_ref,"chain H",num_select);
    std::cout << "We found " << num_select << " atoms on chain H." << endl;

    // Example 2.
    num_select = -1;
    num_select = system_select(aa_ref,"chain D",num_select);
    std::cout << "We found " << num_select << " atoms on chain D." << endl;

    // Example 3.
    num_select = -1;
    num_select = system_select(aa_ref,"chain A",num_select);
    std::cout << "We found " << num_select << " atoms on chain A." << endl;


    // RESIDUE:
    // Example 4.
    num_select = -1;
    num_select = system_select(aa_ref,"resid 539 to 541",num_select);
    std::cout << "We found " << num_select << " atoms for residues 539 to 541" << endl;

    // Example 5.
    num_select = -1;
    num_select = system_select(aa_ref,"resid 397",num_select);
    std::cout << "We found " << num_select << " atoms for residues 397" << endl;


    // INDEX:
    // Example 4.
    num_select = -1;
    num_select = system_select(aa_ref,"index 150 to 350",num_select);
    std::cout << "We found " << num_select << " atoms for indices 150 to 350" << endl;

    // Example 5.
    num_select = -1;
    num_select = system_select(aa_ref,"index 1051",num_select);
    std::cout << "We found " << num_select << " atoms for index 1050" << endl;



    // Working Example.
    num_select = -1;
    const char *pickme = "resid 411 to 418";
    num_select = system_select(aa_ref,pickme,num_select);
    std::cout << "We found " << num_select << " atoms in this selection: " << pickme << endl;


    // Allocate for selection.
    Atom *aa_sel;
    try
    {
        aa_sel = new Atom[num_select];
    }
    catch (std::bad_alloc xa)
    {
        std::cout << "Allocation Failure\n";
        exit(1);
    }


    // Get selection.
    system_select(aa_ref,pickme,num_select,aa_sel);





    // Verify Selection.
    for(int i=0; i<num_select; i++)
    {
        aa_sel[i].num_atoms = num_select;
        // aa_ref[i].print_coords();
        // printf("%s  ",aa_ref[i].chain.c_str());
        // printf("%d\n",aa_ref[i].resid);
        std::cout << "select_i: " << ' '
                  << aa_sel[i].num_atoms << ' '
                  << aa_sel[i].index << ' '
                  << aa_sel[i].resid << ' '
                  << aa_sel[i].chain << ' '
                  << aa_sel[i].restype << ' '
                  << endl;
    }
    /* ---------------------------------------------------------
       End Selection.
       --------------------------------------------------------- */



    /* ---------------------------------------------------------
       For DCD. create frame0 and time-later reference states.
       --------------------------------------------------------- */
    Atom *aa_zero;
    try
    {
        aa_zero = new Atom[num_atoms];
    }
    catch (std::bad_alloc xa)
    {
        std::cout << "Allocation Failure\n";
        exit(1);
    }


    Atom *aa_lat;
    try
    {
        aa_lat = new Atom[num_atoms];
    }
    catch (std::bad_alloc xa)
    {
        std::cout << "Allocation Failure\n";
        exit(1);
    }


    // LOAD DCD.







    std::cout << "\nclosing stdin,stdout,stderr.\n";
    fclose(stdin);
    fclose(stdout);
    fclose(stderr);
    return 0;
}
