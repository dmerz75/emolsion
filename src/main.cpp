// main.cpp

// headers C
extern "C" {
// your functions here for the header

}

// headers C++
// #include <stdlib.h>
// #include <stdio.h> // printf
// #include <string.h> // strcpy, memcpy
// #include <new> // delete
// #include <ctype.h> // getopt - stuff
// #include <unistd.h> // getopt - stuff
#include <vector>
// #include <utility>
// #include <algorithm> // bool & sort
#include <iostream>
// #include <iomanip> // setw
// #include <fstream> //
// #include <map> // map


// boost
// #include "boost/multi_array.hpp"


// my headers
#include "debug.h"
// #include "readfile.h"
// #include "chain.h"
#include "md.h"
#include "ReadPDBfile.hpp"
#include "system.hpp"
#include "dcd.h"
#include "contacts.hpp"
// #include "mt.hpp"
// #include "config.hpp"
// #include "ConfigFile.h"
// #include "Chameleon.h"
// // #include "contacts.h"
// #include "topology.h"
// #include "topology_charmm.h"
// #include "curvature.h"
// #include "chi.h"


// VMD
// #ifndef _DCDIO_H_
// #define _DCDIO_H_
// #include "dcdio.h"
// #endif

// namespaces

int main(int argc, char *argv[]) {

    debug(">> Welcome to SOPCC! (debug)\n");
    printf(">> Welcome to SOPCC!\n");

    if (argc < 3) {
        std::cout << "Usage: " << argv[0] \
                  << " <Filename-reference-PDB>" \
                  << " <Filename-timelater-DCD>"
                  << std::endl;
        exit(1);
    }


    /* ---------------------------------------------------------
       Procedure:
       --------------------------------------------------------- */
    // 1. Read the config file. <conf.sop>
    //    () read_config_file
    //    () alter_default_parameters. Alter the default parameters. <def_param.h>
    // 2. Read the PDB. <2KHO.pdb>
    //    () ReadPDBfile / Count the atoms.
    //    () Allocate for the Reference System. aa_ref
    //    () Populate some System parameters, like the total num_atoms.
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
    // () ReadPDBfile
    // () Count the total "num_atoms".
    Atom a1[0];
    std::cout << "1. Currently, there are " << num_atoms << " atoms." << std::endl;
    num_atoms = ReadPDBfile(argv[1],num_atoms,a1);
    // delete[] a1;
    // delete a1;
    // delete ?
    std::cout << "2. Currently, there are " << num_atoms << " atoms." << std::endl;


    // 2.1 Allocate for the Reference System.
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

    // SUCCESS! Copy constructor implemented!
    Atom *aa_ref2 = aa_ref;
    Atom aa_ref2_1 = aa_ref[0];
    aa_ref2_1.print_coords();
    aa_ref[0].print_coords();
    printf("Copy Constructor example here.\n");
    // exit(0);


    // 2.2 Populate System Parameters
    int chainid, atoms_resid, cur_resid;
    std::string chain;

    chainid = 0;
    chain = aa_ref[0].chain.c_str(); // Get "A"
    // std::cout << "chain: " << chain << std::endl;
    // exit(0);

    // integer list of chains.
    std::vector <int> lst_chainid;
    std::vector <int>::iterator it;

    for(int i=0; i<num_atoms; i++)
    {
        // printf("atom: %d  chain: %d  resid: %d  atoms-resid: %d\n",
        //        aa_ref[i].index,aa_ref[i].chainid,aa_ref[i].resid,
        //        aa_ref[i].num_atoms_res);
        // printf("chain: %s\n",aa_ref[i].chain.c_str());

        if(aa_ref[i].chain.c_str() != chain)
        {
            chainid += 1;
            lst_chainid.push_back(chainid);
        }

        aa_ref[i].chainid = chainid;
        chain = aa_ref[i].chain.c_str();


        // printf("atom: %d  chain: %d  resid: %d  atoms-resid: %d\n",
        //        aa_ref[i].index,aa_ref[i].chainid,aa_ref[i].resid,
        //        aa_ref[i].num_atoms_res);
        // printf("chain: %s\n",aa_ref[i].chain.c_str());
    }

    chainid += 1;
    for(int i=0; i<num_atoms; i++)
    {
        // aa_ref[i].print_coords();
        // printf("%s  ",aa_ref[i].chain.c_str());
        // printf("%d\n",aa_ref[i].resid);
        aa_ref[i].num_atoms = num_atoms;
        aa_ref[i].num_chains = chainid;
    }
    // printf("frames: %d\n",dcd->nsets);




// #ifdef DEBUG // Selection
    /* ---------------------------------------------------------
       Begin Selection:
       Selection: H, precheck!
       --------------------------------------------------------- */
    // system_select() [overloaded] with/without return int (number of matches)
    // system_select(all_atom_sys,"selection string",num_select);
    //          opt: selected_atom_sys);
    int total_H = 0;
    for(int i=0; i<num_atoms; i++)
    {
        if(aa_ref[i].chain.compare("H") == 0)
        {
            total_H += 1;
        }
    }
    std::cout << "total_atoms_on_chain_H(4EZW-55?): " << total_H << std::endl;

    int num_select = -1;

    // CHAIN:
    // num_select = select(aa_ref,"chain","A",num_select);
    // Example 1.
    num_select = -1;
    num_select = system_select(aa_ref,"chain H",num_select);
    std::cout << "1. We found " << num_select << " atoms on chain H." << std::endl;

    // Example 2.
    num_select = -1;
    num_select = system_select(aa_ref,"chain D",num_select);
    std::cout << "2. We found " << num_select << " atoms on chain D." << std::endl;

    // Example 3.
    num_select = -1;
    num_select = system_select(aa_ref,"chain A",num_select);
    std::cout << "3. We found " << num_select << " atoms on chain A." << std::endl;


    // RESIDUE:
    // Example 4.
    num_select = -1;
    num_select = system_select(aa_ref,"resid 539 to 541",num_select);
    std::cout << "4. We found " << num_select << " atoms for residues 539 to 541." << std::endl;

    // Example 5.
    num_select = -1;
    num_select = system_select(aa_ref,"resid 397",num_select);
    std::cout << "5. We found " << num_select << " atoms for residues 397." << std::endl;


    // INDEX:
    // Example 6.
    num_select = -1;
    num_select = system_select(aa_ref,"index 150 to 350",num_select);
    std::cout << "6. We found " << num_select << " atoms for indices 150 to 350." << std::endl;

    // Example 7.
    num_select = -1;
    num_select = system_select(aa_ref,"index 1051",num_select);
    std::cout << "7. We found " << num_select << " atoms for index 1050." << std::endl;

    // Example 8.
    num_select = -1;
    num_select = system_select(aa_ref,"chainid 0",num_select);
    std::cout << "8. We found " << num_select << " atoms for chainid 0." << std::endl;

    // Example 9.
    num_select = -1;
    num_select = system_select(aa_ref,"chainid 3",num_select);
    std::cout << "9. We found " << num_select << " atoms for chainid 3." << std::endl;

    // // Example 10.
    // num_select = -1;
    // num_select = system_select(aa_ref,"chainid 3",num_select);
    // std::cout << "10. We found " << num_select << " atoms for chainid 3." << std::endl;

    // Working Example.
    num_select = -1;
    const char *pickme = "resid 411 to 418";
    num_select = system_select(aa_ref,pickme,num_select);
    std::cout << "We found " << num_select << " atoms in this selection: " << pickme << std::endl;




    /* ---------------------------------------------------------
       Allocate for selection.
       --------------------------------------------------------- */
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



    /* ---------------------------------------------------------
       Verify: aa_sel
       --------------------------------------------------------- */
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
                  << aa_sel[i].x << ' '
                  << aa_sel[i].y << ' '
                  << aa_sel[i].z << ' '
                  << std::endl;

        if(i>4)
        {
            break;
        }
    }

// #endif // Selection



    /* ---------------------------------------------------------
       Vector of chain-vectors of Atoms.
       --------------------------------------------------------- */
// #ifdef NDEBUG
    // Atom *aa_sel;
    // int num_select = 0;
// #endif

    std::vector<std::vector<Atom>> chain_ref;
    // std::vector<std::pair<int,int>> dimers;

    // ITERATORS: itchain, ita for chain, atom.
    std::vector<std::vector<Atom>>::iterator itchain;
    std::vector<Atom>::iterator ita;

    for(int i=0; i<aa_ref[0].num_chains; i++)
    {
        std::vector<Atom> atombin; // Empty?

        // Atom *aa_sel;
        std::string selection = "chainid " + std::to_string(i);
        const char *pickme = selection.c_str();
        // std::cout << i << ' ' << pickme << std::endl;
        num_select = system_select(aa_ref,pickme,num_select);

        // std::cout << i << ' '
        //           << "We found " << num_select
        //           << " atoms for this selection: "
        //           << selection
        //           << " : "
        //           << num_select
        //           << " atoms."
        //           << std::endl;
        // std::cout << "Number of chains in chain_ref(size): " << chain_ref.size() << std::endl;

        try
        {
            aa_sel = new Atom[num_select];
        }
        catch (std::bad_alloc xa)
        {
            std::cout << "Allocation Failure\n";
            exit(1);
        }
        system_select(aa_ref,pickme,num_select,aa_sel);


        for(int j=0; j<num_select; j++)
        {
            // aa_sel[j].print_coords();
            Atom a1 = aa_sel[j]; // Copy Constructor
            // a1.print_coords();
            atombin.push_back(a1);
            // atombin.push_back(&a1);
            // atombin.push_back(aa_sel[j]);

            // if(j>3)
            // {
            //     break;
            // }
        }
        chain_ref.push_back(atombin);
        // printf("\n");
        // std::cout << "binsize: " << atombin.size() << std::endl; // 439, 427, 2..
    }
    std::cout << "# of chains: " << chain_ref.size() << std::endl;

    // exit(0);


    // // EXAMPLE Iteration:
    // for(itchain = chain_ref.begin(); itchain != chain_ref.end(); itchain++)
    // {
    //     std::cout << "Atoms in chain: " << (*itchain).size() << std::endl;

    //     for(ita = (*itchain).begin(); ita != (*itchain).end(); ita++)
    //     {
    //         (*ita).print_coords();
    //     }
    // }
    // // exit(0);

    /* ---------------------------------------------------------
       Vector of chain-vectors of Atoms. END.
       --------------------------------------------------------- */





#ifdef MTMAP
    // Pairs: A and B.
    std::vector<std::pair<int,int>> dimers;

    // ACCESS chain_ref
    int amon, bmon;
    amon = 0;
    bmon = 1;

    for(itchain = chain_ref.begin(); itchain != chain_ref.end(); itchain++)
    {
        if((*itchain).size() >= 433)
        {
            dimers.push_back(std::make_pair(amon,bmon));
            amon += 2;
            bmon += 2;
        }
    }

    // // DIMERS
    // std::vector<std::pair<int,int>>::iterator itdimers;
    // for(itdimers = dimers.begin(); itdimers != dimers.end(); itdimers++)
    // {
    //     std::cout << (*itdimers).first << ' ' << (*itdimers).second << std::endl;
    // }
    // exit(0);


    // External Neighbor (chainid).
    // std::vector<std::vector<int>> mt_matrix(aa_ref[0].num_chains, std::vector<int>(8,-1));


    // std::vector<std::vector<Atom>>::iterator itchain;
    // std::vector<Atom>::iterator ita;

    std::vector<std::vector<int>> mt_matrix(dimers.size(), std::vector<int>(8,-1));

    // Get Map of MT neighbors.
    mt_matrix = get_map_of_mtneighbors(chain_ref,dimers);
    std::vector<std::vector<int>>::iterator itmap;
    std::vector<int>::iterator itmap_n;

    for(itmap = mt_matrix.begin(); itmap != mt_matrix.end(); itmap++)
    {
        // std::cout << (*itmap).size() << std::endl;

        for(itmap_n = (*itmap).begin(); itmap_n != (*itmap).end(); itmap_n++)
        {
            std::cout << (*itmap_n) << " ";
        }
                std::cout << std::endl;
    }



    // std::vector<std::vector<int>>::iterator itm_cid; // for chainid
    // std::vector<std::vector<int>>::iterator itm_n; // for the neighbor
    // std::cout << "0,0 \t" <<mt_matrix[0][0] << '\n'
    //           << "150,150 \t" << mt_matrix[150][150] << '\n'
    //           << std::endl;

    // for(itm_cid = mt_matrix.begin(); itm_cid != mt_matrix.end(); itm_cid++)
    // {
    //     std::cout << (*itm_cid).size() << std::endl;
    //     std::cout << (*itm_cid)[0] << std::endl;
    //     std::cout << (*itm_cid)[2] << std::endl;
    //     std::cout << (*itm_cid)[3] << std::endl;
    //     std::cout << (*itm_cid)[6] << std::endl;
    //     std::cout << (*itm_cid)[7] << std::endl;

    //     // for(itm_n = (*itm_cid)->begin(); itm_n != (*itm_cid).end(); itm_n++)
    //     // {
    //     //     printf("hello\n.");
    //     // }
    //     // std::cout <<
    // }





    // for(auto vec: mt_matrix)
    // {
    //     // std::cout << vec[0] << vec[1] << vec[2] << vec[3] << std::endl;
    //     for(auto x: vec)
    //     {
    //         std::cout << x << std::endl;
    //     }
    //     std::cout << "\n" << std::endl;
    // }






#endif // MTMAP
    /* ---------------------------------------------------------
       End Selection.
       --------------------------------------------------------- */





    /* ---------------------------------------------------------
       Create aa_zero, aa_later reference states.
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

    Atom *aa_later;
    try
    {
        aa_later = new Atom[num_atoms];
    }
    catch (std::bad_alloc xa)
    {
        std::cout << "Allocation Failure\n";
        exit(1);
    }

    //
    system_select(aa_ref,"all",num_atoms,aa_zero);
    system_select(aa_ref,"all",num_atoms,aa_later);


    // Verify aa_zero and aa_later.
    for(int i=0; i<num_atoms; i++)
    {
        aa_zero[i].num_atoms = num_atoms;
        aa_later[i].num_atoms = num_atoms;

        // aa_ref[i].print_coords();
        // printf("%s  ",aa_ref[i].chain.c_str());
        // printf("%d\n",aa_ref[i].resid);
        // std::cout << "select_i: " << ' '
        //           << aa_sel[i].num_atoms << ' '
        //           << aa_sel[i].index << ' '
        //           << aa_sel[i].resid << ' '
        //           << aa_sel[i].chain << ' '
        //           << aa_sel[i].restype << ' '
        //           << std::endl;
    }



#if defined (DCDREAD) || defined (DCD_WRITE_B) || defined (DCD_WRITE) || defined (DCD_WRITE_E)

    /* ---------------------------------------------------------
       Analysis Before. Start.
       --------------------------------------------------------- */
    int someindex = 0;
    someindex = aa_ref[0].num_atoms / 2;
    debug("middle-index: %d\n",someindex);
    debug("coords: %f %f %f\n",aa_ref[someindex].x,aa_ref[someindex].y,aa_ref[someindex].z);
    debug("coords: %f %f %f\n",aa_zero[someindex].x,aa_zero[someindex].y,aa_zero[someindex].z);
    debug("coords: %f %f %f\n",aa_later[someindex].x,aa_later[someindex].y,aa_later[someindex].z);

    debug("\n");
    debug("coords(8)[0]: %f %f %f\n",aa_zero[8].x,aa_zero[8].y,aa_zero[8].z);
    debug("coords(8)[later]: %f %f %f\n",aa_later[8].x,aa_later[8].y,aa_later[8].z);
    debug("coords(37)[0]: %f %f %f\n",aa_zero[37].x,aa_zero[37].y,aa_zero[37].z);
    debug("coords(37)[later]: %f %f %f\n",aa_later[37].x,aa_later[37].y,aa_later[37].z);
#endif // multi-dcd




#ifdef GET_CONTACTS
    // void get_contacts(Atom *a1,Atom *a2,char dcdfilename[40],int num_atoms);
    // get_contacts(aa_sel,aa_sel,argv[2],num_atoms);


    // get_contacts(aa_sel,argv[2],num_atoms);

#endif // GET_CONTACTS

    /* ---------------------------------------------------------
       Analysis Before. Finish.
       --------------------------------------------------------- */





    /* ---------------------------------------------------------
       DCD Read. Preload.
       --------------------------------------------------------- */
#ifdef DCDREAD
    int frame_position;
    frame_position = 0;
    int start,stop,step;

    // start,stop,step
    start = atoi(argv[3]);
    stop = atoi(argv[4]);
    step = atoi(argv[5]);



    /* ---------------------------------------------------------
       DCD Write. Preload.
       --------------------------------------------------------- */
#ifdef DCD_WRITE_B
    std::string str_dcd_read(argv[2]);
    std::size_t found;
    found=str_dcd_read.find(".dcd",0);

    std::string str_dcd_name = str_dcd_read.substr(0,found);
    std::string str_dcd_write = str_dcd_name + "_subset.dcd";
    char *fn_dcd_write = (char*) str_dcd_write.c_str();

    printf("\nreading dcd >>>  <%s>\n",str_dcd_read.c_str());
    printf("using name for dcd >>>  <%s>\n",str_dcd_name.c_str());
    printf("writing dcd >>>  <%s>\n",fn_dcd_write);

    // Write DCD
    molfile_timestep_t timestep_w;
    void *vw;
    dcdhandle *dcdw;
    int natoms_w = 0;

    dcdw = (dcdhandle *)vw;

    // // get atom total for writing.
    // for(int i=0; i<num_chains; i++){
    //     natoms_w += chain_ref[i].num_atoms_ca;
    // }
    // printf("for dcd writing >>>  <%d> atoms expected.\n",natoms_w);
    printf("for dcd writing >>>  <%d> atoms expected.\n",num_atoms);


    // vw = open_dcd_write(fn_dcd_write,"dcd",natoms_w);
    vw = open_dcd_write(fn_dcd_write,"dcd",num_atoms);
    if (!vw) {
        fprintf(stderr, "main) open_dcd_write failed for file %s\n", *fn_dcd_write);
        return 1;
    } else {
        printf("opened <%s> successfully!!\n\n",fn_dcd_write);
    }

    // timestep_w.coords = (float *)malloc(3*sizeof(float)*natoms_w);
    timestep_w.coords = (float *)malloc(3*sizeof(float)*num_atoms);

    // dcd = (dcdhandle *)v;
    // sizeMB = ((natoms * 3.0) * dcd->nsets * 4.0) / (1024.0 * 1024.0);
    // totalMB += sizeMB;
    // printf("main) file: %s\n", *argv);
    // printf("  %d atoms, %d frames, size: %6.1fMB\n", natoms, dcd->nsets, sizeMB);
    // timestep.coords = (float *)malloc(3*sizeof(float)*natoms);
#endif // DCD_WRITE_B


    /* ---------------------------------------------------------
       DCD Read.
       --------------------------------------------------------- */
    molfile_timestep_t timestep;
    void *v;
    dcdhandle *dcd;
    int natoms; // from the opening the dcd.
    float sizeMB =0.0, totalMB = 0.0;
    double starttime, endtime, totaltime = 0.0;
    // typedef struct {
    //     fio_fd fd;
    //     int natoms; 382
    //     int nsets; 17501 (0-17500)
    //     int setsread;
    //     int istart;
    //     int nsavc;
    //     double delta;
    //     int nfixed;
    //     float *x, *y, *z; ->x[0-381];
    //     int *freeind;
    //     float *fixedcoords;
    //     int reverse;
    //     int charmm;
    //     int first;
    //     int with_unitcell;
    // } dcdhandle;



    // // 1
    // while (--argc) {
    //     ++argv;
    //     natoms = 0;
    //     v = open_dcd_read(*argv, "dcd", &natoms);
    //     if (!v) {
    //         fprintf(stderr, "main) open_dcd_read failed for file %s\n", *argv);
    //         return 1;
    //     }
    //     dcd = (dcdhandle *)v;
    //     sizeMB = ((natoms * 3.0) * dcd->nsets * 4.0) / (1024.0 * 1024.0);
    //     totalMB += sizeMB;
    //     printf("main) file: %s\n", *argv);
    //     printf("  %d atoms, %d frames, size: %6.1fMB\n", natoms, dcd->nsets, sizeMB);

    //     // starttime = time_of_day();
    //     timestep.coords = (float *)malloc(3*sizeof(float)*natoms);
    //     for (i=0; i<dcd->nsets; i++) {
    //         int rc = read_next_timestep(v, natoms, &timestep);
    //         if (rc) {
    //             fprintf(stderr, "error in read_next_timestep on frame %d\n", i);
    //             return 1;
    //         }
    //     }
    //     // endtime = time_of_day();
    //     close_file_read(v);
    //     // totaltime += endtime - starttime;
    //     // printf("  Time: %5.1f seconds\n", endtime - starttime);
    //     // printf("  Speed: %5.1f MB/sec, %5.1f timesteps/sec\n", sizeMB \
    //     //        / (endtime - starttime), (dcd->nsets / (endtime - starttime)));
    // }
    // printf("Overall Size: %6.1f MB\n", totalMB);
    // // printf("Overall Time: %6.1f seconds\n", totaltime);
    // // printf("Overall Speed: %5.1f MB/sec\n", totalMB / totaltime);


    printf("----->  READING DCD  <-----\n");
    // int atoms_in_chain;
    // 2. to read a dcd.
    natoms = 0;
    v = open_dcd_read(argv[2],"dcd",&natoms);
    if (!v)
    {
        fprintf(stderr, "main) open_dcd_read failed for file %s\n", *argv);
        return 1;
    }
    dcd = (dcdhandle *)v;
    sizeMB = ((natoms * 3.0) * dcd->nsets * 4.0) / (1024.0 * 1024.0);
    totalMB += sizeMB;
    timestep.coords = (float *)malloc(3*sizeof(float)*natoms);
    // timestep.velocities = (float *)malloc(3*sizeof(float)*natoms);

    // printf("main) file: %s\n", *argv);
    // printf("  %d atoms, %d frames, size: %6.1fMB\n", natoms, dcd->nsets, sizeMB);
    // close_file_read(v);







    /* ---------------------------------------------------------
       Analysis Before. But with DCD Open. Begin.
       --------------------------------------------------------- */

#ifdef CONTACTS_BEFORE




#endif // CONTACTS_BEFORE




    /* ---------------------------------------------------------
       Analysis Before. But with DCD Open. End.
       --------------------------------------------------------- */


    frame_position = 1;
    advance_dcd(dcd->nsets,0,dcd,natoms,&timestep); // 1st advance. 1-vmd

    // THIS ONE
    // load_dcd_to_chain(dcd,aa_zero,num_chains);
    load_dcd_to_atoms(dcd,aa_zero);

    // printf("%f %f %f\n",)

    // printf("frame_position: %d\n",frame_position);
    // printf("ref-findex(%d): %f\n",chain_ref[0].findex,chain_ref[0].pos[chain_ref[0].findex].y);
    // printf("0-findex(%d): %f\n",chain_0[0].findex,chain_0[0].pos[chain_0[0].findex].y);
    // printf("later-findex(%d): %f\n",chain_later[0].findex,chain_later[0].pos[chain_later[0].findex].y);

    frame_position = 2;
    advance_dcd(dcd->nsets,0,dcd,natoms,&timestep); // 2nd. 2-vmd

    // THIS ONE
    // load_dcd_to_chain(dcd,chain_later,num_chains);
    load_dcd_to_atoms(dcd,aa_later);

    printf("frame_position: %d\n",frame_position);
    // exit(0);

    // Advancing Rules.
    // ----------------
    // example. step size -> 5.
    // int advance_size = atoi(argv[3]) - 1;

    // step.
    int advance_size = step - 1; // 0 counts, so advance_size of 4, advances by 5.

    // stop.
    if(stop > dcd->nsets){
        stop = dcd->nsets;
        printf("use stop value: %d\n",stop);
    }

    // // start.
    // debug("starting frame: %d\n",start);
    // if((start > 2) && (step < start)) {
    //     // for (int nset1=2; nset1<start; nset1 += advance_size + 1 ) {
    //     for (int nset1=2; nset1<start-step; nset1 += 1 ) {
    //         // for (int nset1=2; nset1<dcd->nsets; nset1 += step + 1) {
    //         debug("forwarding --> current: %d\n",nset1);
    //         frame_position += advance_dcd(dcd->nsets,0,dcd,natoms,&timestep);
    //         debug("forwarding --> frame_position: %d\n",frame_position);

    //         // if ( nset1 + advance_size >= start) {
    //         //     // if ( nset1 + step >= dcd->nsets ) {
    //         //     break;
    //         // }
    //         // frame_position is returned value;
    //         // frame_position += advance_dcd(dcd->nsets,step,dcd,natoms,&timestep);
    //         // load_dcd_to_chain(dcd,chain_later,num_chains);
    //     }
    // } else if ((start > 2) && (step > start - 2)) {
    //     for (int nset1=2; nset1<start; nset1 += 1 ) {
    //         debug("forwarding --> current: %d\n",nset1);
    //         frame_position += advance_dcd(dcd->nsets,0,dcd,natoms,&timestep);
    //         debug("forwarding --> frame_position: %d\n",frame_position);
    //     }
    // }

    // sleep(0.5);
    // exit(0);


    for (int nset1=2; nset1<start; nset1 += 1 ) {
        // for (int nset1=2; nset1<dcd->nsets; nset1 += step + 1) {
        // debug("forwarding --> current: %d\n",nset1);
        frame_position += advance_dcd(dcd->nsets,0,dcd,natoms,&timestep);


        // THIS ONE
        // load_dcd_to_chain(dcd,chain_later,num_chains);
        load_dcd_to_atoms(dcd,aa_later);

        debug("forwarding --> frame_position: %d\n",frame_position);
    }

    // Get initial starting point.
    printf("--> fast forwarded. to frame: %d\n",frame_position);


    /* ---------------------------------------------------------
       major for loop begin || doloop.
       --------------------------------------------------------- */
    // for (int nset2=frame_position; nset2<stop; nset2 += advance_size + 1) {
    //     // for (int nset=frame_position; nset<dcd->nsets; nset += advance_size + 1 ) {

    //     if ( nset2 + advance_size >= stop) {
    //         // if ( nset + advance_size >= dcd->nsets ) {
    //         break;
    //     }

    //     frame_position += advance_dcd(dcd->nsets,advance_size,dcd,natoms,&timestep);
    //     load_dcd_to_chain(dcd,chain_later,num_chains);
    //     debug("current: %d\n",nset2);
    //     printf("--> frame_position: %d\n",frame_position);
    //     // check.
    //     // printf("ref-findex(%d): %f\n",chain_ref[0].findex,chain_ref[0].pos[chain_ref[0].findex].y);
    //     // printf("0-findex(%d): %f\n",chain_0[0].findex,chain_0[0].pos[chain_0[0].findex].y);
    //     // printf("later-findex(%d): %f\n",chain_later[0].findex,chain_later[0].pos[chain_later[0].findex].y);
    //     // continue;

    int nset2;
    nset2 = frame_position;
    do {

#endif //DCDREAD



#if defined (DCDREAD) || defined (DCD_WRITE_B) || defined (DCD_WRITE) || defined (DCD_WRITE_E)
        /* ---------------------------------------------------------
           Analysis During. Start.
           --------------------------------------------------------- */
        debug("middle-index: %d\n",someindex);
        debug("coords: %f %f %f\n",aa_ref[someindex].x,aa_ref[someindex].y,aa_ref[someindex].z);
        debug("coords: %f %f %f\n",aa_zero[someindex].x,aa_zero[someindex].y,aa_zero[someindex].z);
        debug("coords: %f %f %f\n",aa_later[someindex].x,aa_later[someindex].y,aa_later[someindex].z);

        debug("\n");
        debug("coords(8)[0]: %f %f %f\n",aa_zero[8].x,aa_zero[8].y,aa_zero[8].z);
        debug("coords(8)[%d]: %f %f %f\n",frame_position,aa_later[8].x,aa_later[8].y,aa_later[8].z);
        debug("coords(37)[0]: %f %f %f\n",aa_zero[37].x,aa_zero[37].y,aa_zero[37].z);
        debug("coords(37)[%d]: %f %f %f\n",frame_position,aa_later[37].x,aa_later[37].y,aa_later[37].z);







#ifdef CONTACTS_DURING

#endif // CONTACTS_DURING







        /* ---------------------------------------------------------
           Analysis During. Finish.
           --------------------------------------------------------- */
#endif // multi-dcd



#ifdef DCD_WRITE
        // READ
        // static void *open_dcd_read(const char *path, const char *filetype,
        //                            int *natoms) {
        // int rc = read_next_timestep(v, natoms, timestep);

        // WRITE
        // static void *open_dcd_write(const char *path, const char *filetype,
        //                             int natoms) {
        // static void close_file_write(void *v) {

        // open_dcd_write(fn_dcd_write,"dcd",&natoms);

        // write_dcdstep
        // static int write_dcdstep(fio_fd fd, int curframe, int curstep, int N,
        //                          const float *X, const float *Y, const float *Z,
        //                          const double *unitcell, int charmm) {

        // static int write_dcdheader(fio_fd fd, const char *remarks, int N,
        //                            int ISTART, int NSAVC, double DELTA, int with_unitcell,
        //                            int charmm) {

        // load_chain_coords_to_timestep(dcdw,chain_later,chains_to_use,vw,&timestep_w,natoms_w);
        // load_chain_coords_to_timestep(chain_later,num_chains,&timestep_w);



        // THIS ONE
        // load_chain_to_timestep(chain_later,num_chains,&timestep_w);
        load_atom_to_timestep(&timestep_w,aa_later);



#ifdef DCD_WRITE_UNMOD
    // Write the DCD read in.
        write_timestep(vw,&timestep);

    // #elif DCD_WRITE_ROT
#else
        // Write modified coordinates.
        write_timestep(vw,&timestep_w);
#endif // A
#endif



#ifdef DCDREAD
    // } // DCD PRIMARY LOOP

        debug("current: %d\n",nset2);
        printf("frame: --> %d <-- was evaluated.\n",frame_position);

        if (nset2 + advance_size + 1 <= stop)
        {
            frame_position += advance_dcd(dcd->nsets,advance_size,dcd,natoms,&timestep);
            printf("frame: --> %d <-- loaded.\n",frame_position);

            // THIS ONE
            // load_dcd_to_chain(dcd,chain_later,num_chains);
            load_dcd_to_atoms(dcd,aa_later);

            // nset2 += advance_size + 1;
        }
        nset2 += advance_size + 1;
    } while (nset2<=stop);

    debug("\n..closing dcd..\n");
    close_file_read(v);


    printf("\n----->  READING DCD completed!  <-----\n");
    printf("\t\tThe maximum possible frame_position was: %d\n",stop);
    printf("\t\tThe last frame evaluated was: %d\n",frame_position);
#endif //DCDREAD



#if defined (DCDREAD) || defined (DCD_WRITE_B) || defined (DCD_WRITE) || defined (DCD_WRITE_E)
    /* ---------------------------------------------------------
       Analysis After. Start.
       --------------------------------------------------------- */
    debug("middle-index: %d\n",someindex);
    debug("coords: %f %f %f\n",aa_ref[someindex].x,aa_ref[someindex].y,aa_ref[someindex].z);
    debug("coords: %f %f %f\n",aa_zero[someindex].x,aa_zero[someindex].y,aa_zero[someindex].z);
    debug("coords: %f %f %f\n",aa_later[someindex].x,aa_later[someindex].y,aa_later[someindex].z);


    debug("\n");
    debug("coords(8)[0]: %f %f %f\n",aa_zero[8].x,aa_zero[8].y,aa_zero[8].z);
    debug("coords(8)[%d]: %f %f %f\n",frame_position,aa_later[8].x,aa_later[8].y,aa_later[8].z);
    debug("coords(37)[0]: %f %f %f\n",aa_zero[37].x,aa_zero[37].y,aa_zero[37].z);
    debug("coords(37)[%d]: %f %f %f\n",frame_position,aa_later[37].x,aa_later[37].y,aa_later[37].z);


#ifdef CONTACTS_AFTER

#endif // CONTACTS_AFTER



    /* ---------------------------------------------------------
       Analysis After. Finish.
       --------------------------------------------------------- */
#endif // multi-dcd


#ifdef DCD_WRITE_E
    // open_dcd_write(fn_dcd_write,"dcd",natoms);
    // static void close_file_write(void *v) {
    close_file_write(vw);
#endif // DCD_WRITE_E




    /* ---------------------------------------------------------
       Delete Malloc Systems.
       --------------------------------------------------------- */
    delete [] aa_ref;
    delete [] aa_zero;
    delete [] aa_later;
    delete [] aa_sel;



    /* ---------------------------------------------------------
       The End.
       --------------------------------------------------------- */
    std::cout << "\nclosing stdin,stdout,stderr.\n";
    fclose(stdin);
    fclose(stdout);
    fclose(stderr);
    return 0;
}
