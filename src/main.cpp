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
#include <iomanip> // setw
// #include <fstream> //
// #include <map> // map
#include "boost/tuple/tuple.hpp"
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
#include "microtubule.hpp"
#include "phipsiangle.hpp"
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

    printf(">> Welcome to Emolsion!\n");
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] \
                  << " <Filename-reference-PDB>" \
                  << " <Filename-timelater-DCD>"
                  << std::endl;
        exit(1);
    }

    /* ---------------------------------------------------------
       Procedure: Steps.
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
    // 4. Evaluate, run MD.
    //    4.1 Write DCD every so often.
    //    4.2 Include PDBs
    //    4.3 Evaluate.
    // 5. Write.
    //    5.1 Write PDBs?
    //    5.2 Write output files.
    // 6. Close, Free up memory.


    /* ---------------------------------------------------------
       Step 1. read the config file!
       --------------------------------------------------------- */


    /* ---------------------------------------------------------
       Step 2. Read the pdb.
       --------------------------------------------------------- */
    // 2.1.1 Copy for the Reference state. (PDB)
    // 2.1.2 Copy for DCD-0.
    // 2.1.3 Pointers for aa_later.

    // 2.1.1 allatoms_ref (PDB)
    vAtoms allatoms_ref;
    allatoms_ref = ReadPDBfile(argv[1]);

    int num_atoms;
    num_atoms = allatoms_ref.size();
    // std::cout << "Number of atoms in allatoms_ref: " << allatoms_ref.size() << std::endl;

    // Reserve space for the standards.
    // allatoms_ref  -->  allatoms_0, allatoms.
    vAtoms allatoms_0; // from DCD-0
    vAtoms allatoms; // evolve in time.
    allatoms_0.reserve(num_atoms);
    allatoms.reserve(num_atoms);

    // Set the chainid based on switch from chain A to chain B .. etc.
    allatoms_ref = set_chainid(allatoms_ref);

    // Check index:
    // for(auto a: allatoms_ref)
    // {
    //     std::cout << a.index << std::endl;
    // }

    // Check chainid:
    // for(auto a: allatoms_ref)
    // {
    //     std::cout << "chainid: " << a.chainid
    //               << " chain: " << a.chain
    //               << std::endl;
    // }
    // exit(0);


    // Duplication: Initial set of copies:
    // 1. allatoms_ref (PDB)
    // 1. allatoms (time-evolving)
    // 2. allatoms_0 (from DCD-0).
    for(int j=0; j<allatoms_ref.size(); j++)
    {
        Atom a = allatoms_ref[j];
        // a.print_coords();
        allatoms.push_back(a);
        allatoms_0.push_back(a);
    }

    // // Check 3 standards:
    // for(int j=0; j<allatoms_ref.size(); j++)
    // {
    //     // allatoms_ref[j].print_Atom();
    //     allatoms_ref[j].print_coords();
    //     allatoms[j].print_coords();
    //     allatoms_0[j].print_coords();
    // }
    // // exit(0);

    // // Check chainid.
    // for(auto a: allatoms)
    // {
    //     std::cout << "chainid: " << a.chainid
    //               << " chain: " << a.chain
    //               << std::endl;
    // }
    // exit(0);


#ifndef NDEBUG
    /* ---------------------------------------------------------
       Select: vector of Atoms
       --------------------------------------------------------- */
    // select()
    IndexGroup chainH;
    chainH = select(allatoms,"chain H");
    std::cout << "chain H: " << chainH.size() << std::endl;

    IndexGroup chainD;
    chainD = select(allatoms,"chain D");
    std::cout << "chain D: " << chainD.size() << std::endl;

    IndexGroup chainA;
    chainA = select(allatoms,"chain A");
    std::cout << "chain A: " << chainA.size() << std::endl;

    IndexGroup resid5;
    resid5 = select(allatoms,"resid 539 to 541");
    std::cout << "resid 539 to 541: " << resid5.size() << std::endl;

    IndexGroup resid3;
    resid3 = select(allatoms,"resid 37");
    std::cout << "resid 37: " << resid3.size() << std::endl;

    IndexGroup sel150;
    sel150 = select(allatoms,"index 150 to 350");
    std::cout << "Selection: index 150 to 350: " << sel150.size() << std::endl;

    IndexGroup sel597;
    sel597 = select(allatoms,"index 597");
    std::cout << "Selection: index 597: " << sel597.size() << std::endl;

    IndexGroup sel0;
    sel0 = select(allatoms,"chainid 0");
    std::cout << "Selection: chainid 0: " << sel0.size() << std::endl;

    IndexGroup sel1;
    sel1 = select(allatoms,"chainid 1");
    std::cout << "Selection: chainid 1: " << sel1.size() << std::endl;

    IndexGroup sel_all;
    sel_all = select(allatoms,"all");
    std::cout << "Selection: (all)  " << sel_all.size() << std::endl;
    // exit(0);

    IndexGroup sel_calpha;
    sel_calpha = select(allatoms,"atomtype CA");
    std::cout << "Selection: (CALPHA) " << sel_calpha.size() << std::endl;
#endif // NDEBUG Selection.

    // Primary Selection:
    // Organization by chains: isel_chain.

    // Build allatoms_chain, as indices from allatoms, but sorted by chain.
    vIndexGroup isel_chain;
    isel_chain = sort_segment_chain(allatoms);
    // for(auto c: isel_chain)
    // {
    //     std::cout << c.size() << std::endl;
    //     // for(auto i: c)
    //     // {
    //     //     std::cout << i << " ";
    //     // }
    // }
    // // exit(0);

    /* ---------------------------------------------------------
       End of Selection.
       --------------------------------------------------------- */


#if defined (DCDREAD) || defined (DCD_WRITE_B) || defined (DCD_WRITE) || defined (DCD_WRITE_E)

    /* ---------------------------------------------------------
       Analysis Before. Start.
       Before DCD Read is open. Just checking coordinates at this point.
       --------------------------------------------------------- */

#ifndef NDEBUG // Double Negative Define: Not No-debugging.
    // Position Check with
    int atomcheck = 0;
    atomcheck = allatoms_ref.size() * 0.5;
    if(atomcheck == 0)
    {
        atomcheck = 1;
    }
    else if(atomcheck > 5)
    {
        atomcheck = 5;
    }

    std::cout << "Allatoms-Ref:" << std::endl;
    for(int i=0; i<atomcheck; i++)
    {
        allatoms_ref[i].print_Coords();
    }

    std::cout << "Allatoms-t:" << std::endl;
    for(int i=0; i<atomcheck; i++)
    {
        allatoms[i].print_Coords();
        // std::cout << "\t" << allatoms[i].x << " "
        //           << "\t" << allatoms[i].y << " "
        //           << "\t" << allatoms[i].z << " "
        //           << std::endl;
    }
#endif // NDEBUG

#ifdef MTMAP2_BEFORE
    // Pairs: A(~439) and B(427-8).
    DimerList dimers;
    MtIndexMap mtmap_subdomain;
    MtIndexMapEntry mtentry;

    // ACCESS chain_ref
    int imonomer = 0;
    int betabool = -1;
    int i_count = -1;
    int high_index, low_index, Nterm2, Mterm1, Mterm2, Cterm1;


    // std::cout << "exiting .." << std::endl;
    // exit(0);
    // std::cout << "mtmap-size: " << mtmap_subdomain[0].size() << std::endl;

    for(auto c: isel_chain)
    {
        // std::cout << "ChainSize: " << c.size() << std::endl;
        high_index = low_index = Nterm2 = Mterm1 = Mterm2 = Cterm1 = -1;
        i_count += 1;
        betabool = -1;
        imonomer += 1;

        if((c.size() >= 433) and (c.size() <= 442))
        {
            // std::cout << "Alpha." << std::endl;
            dimers.push_back(std::make_pair(imonomer-1,imonomer));
            betabool = 0;
        }
        else if ((c.size() > 420) and (c.size() < 433))
        {
            // std::cout << "Beta." << std::endl;
            betabool = 1;
        }
        else
        {
            // std::cout << "Neither." << std::endl;
            continue;
        }
        // exit(0);

        // high_index = -1;
        for(auto i:c)
        {
            if(i > high_index)
            {
                high_index = i;
            }
        }

        low_index = high_index - c.size() + 1;

        if((c.size() > 420) and (c.size() < 442))
        {
            Nterm2 = low_index + 214;
            Mterm1 = Nterm2 + 1;   // 215
            Mterm2 = Mterm1 + 168; // 383
            Cterm1 = Mterm2 + 1;   // 427/438
        }

        // std::cout << "chainid: " << imonomer - 1 << " \n"
        //           << "chainsize: " << c.size() << " \n"
        //           << "low|high: " << low_index << " "
        //           << high_index << " "
        //           << "  " << betabool << "  (0-alpha,1-beta)\n"
        //           << "N-M-C: \n"
        //           << "\t" << low_index << " " << Nterm2 << " \n"
        //           << "\t" << Mterm1 << " " << Mterm2 << " \n"
        //           << "\t" << Cterm1 << " " << high_index << " \n"
        //           // << "chainid: " << c[0].chainid << " \n " // -1
        //           // <<
        //           << " "
        //           << std::endl;

        mtentry["chaintype"] = betabool;
        mtentry["index"] = low_index;
        mtentry["findex"] = high_index;
        mtentry["Nterm2"] = Nterm2;
        mtentry["Mterm1"] = Mterm1;
        mtentry["Mterm2"] = Mterm2;
        mtentry["Cterm1"] = Cterm1;

        mtmap_subdomain.push_back(mtentry);
        mtentry.clear();
        // mtdemarcations.push_back()
    }
    // exit(0);

    // // PRINT mtmap
    // for(auto m:mtmap_subdomain)
    // {
    //     std::cout << "chaintype: " << m["chaintype"] << " \n"
    //               << "index: " << m["index"] << " \n"
    //               << "Nterm2: " << m["Nterm2"] << " \n"
    //               << "Mterm1: " << m["Mterm1"] << " \n"
    //               << "Mterm2: " << m["Mterm2"] << " \n"
    //               << "Cterm1: " << m["Cterm1"] << " \n"
    //               << "findex: " << m["findex"] << " \n"
    //               << std::endl;
    // }
    // // exit(0);


    // Print DimerList dimers.
    // imonomer = 0;
    // for(auto d: dimers)
    // {
    //     std::cout << dimers[imonomer].first << " " << dimers[imonomer].second << std::endl;
    //     imonomer += 1;
    // }
    // exit(0);


    // // DIMERS
    // std::vector<std::pair<int,int>>::iterator itdimers;
    // for(itdimers = dimers.begin(); itdimers != dimers.end(); itdimers++)
    // {
    //     std::cout << (*itdimers).first << ' ' << (*itdimers).second << std::endl;
    // }
    // exit(0);


    /* ---------------------------------------------------------
       MTMAP: set a matrix for the map, populate that matrix.
       --------------------------------------------------------- */
    // Vector of Vector <int>
    MtNeighbors mt_matrix(dimers.size(), std::vector<int>(8,-1));
    // Get Map of MtNeighbors.
    // mt_matrix = get_map_of_mtneighbors(chain_ref,dimers);
    // mt_matrix = get_map_of_mtneighbors(&allatoms_chain,dimers);
    // mt_matrix = get_map_of_mtneighbors(allatoms_chain,dimers);
    std::cout << "Getting map for microtubule." << std::endl;
    mt_matrix = get_map_of_mtneighbors(isel_chain,allatoms_ref,dimers);

    // Print Map of MT neighbors.
    // std::cout << "Printing map of microtubule neighbors." << std::endl;
    // print_mt_map(mt_matrix);
    // exit(0);

    std::cout << "MTMAP2: Beginning contacts by sector." << std::endl;

    SetContacts contact_set;
    SetNeighbors neighbor_set;
    SetChains chain_set;
    SetGlobalContacts global_contacts;

    for(auto c: mt_matrix)
    {
        neighbor_set.clear();

        // Alpha, Beta, Alpha-Beta
        // contact_set = get_contacts_for_chain(chain_ref[c[0]],chain_ref[c[0]],8.0);
        // neighbor_set.push_back(contact_set);
        // contact_set.clear();
        // contact_set = get_contacts_for_chain(chain_ref[c[1]],chain_ref[c[1]],8.0);
        // neighbor_set.push_back(contact_set);
        // contact_set.clear();

        // contact_set = get_contacts_for_chain(allatoms_ref[c[0]],8.0,mtmap_subdomain,c[0]);
        // contact_set = get_contacts_for_chain(allatoms_ref[c[1]],8.0,mtmap_subdomain,c[1]);
        // contact_set = get_contacts_for_chain(allatoms_ref[c[0]],allatoms_ref[c[1]],8.0,
        //                                      mtmap_subdomain,
        //                                      c[0],
        //                                      c[1]);

        contact_set = get_contacts_for_chain(allatoms_ref,8.0,
                                             mtmap_subdomain,
                                             isel_chain[c[0]],
                                             c[0]);
        neighbor_set.push_back(contact_set);
        contact_set.clear();
        contact_set = get_contacts_for_chain(allatoms_ref,8.0,
                                             mtmap_subdomain,
                                             isel_chain[c[1]],
                                             c[1]);
        neighbor_set.push_back(contact_set);
        contact_set.clear();

        contact_set = get_contacts_for_chain(allatoms_ref,8.0,
                                             mtmap_subdomain,
                                             isel_chain[c[0]],
                                             isel_chain[c[1]],
                                             c[0],c[1]);
        neighbor_set.push_back(contact_set);
        contact_set.clear();
    // }
    // exit(0);

        for(int m=2; m<=4; m++)
        {
            if(c[m] < 0)
            {
                contact_set.clear();
                neighbor_set.push_back(contact_set);
                continue;
            }
            // contact_set = get_contacts_for_chain(allatoms_ref[c[0]],
            //                                      allatoms_ref[c[m]],
            //                                      8.0,
            //                                      mtmap_subdomain,
            //                                      c[0],
            //                                      c[m]);
            contact_set = get_contacts_for_chain(allatoms_ref,
                                                 8.0,
                                                 mtmap_subdomain,
                                                 isel_chain[c[0]],
                                                 isel_chain[c[m]],
                                                 c[0],
                                                 c[m]);
            neighbor_set.push_back(contact_set);
            contact_set.clear();
        }

        for(int m=5; m<=7; m++)
        {
            if(c[m] < 0)
            {
                contact_set.clear();
                neighbor_set.push_back(contact_set);
                continue;
            }
            contact_set = get_contacts_for_chain(allatoms_ref,
                                                 8.0,
                                                 mtmap_subdomain,
                                                 isel_chain[c[0]],
                                                 isel_chain[c[m]],
                                                 c[1],
                                                 c[m]);
            neighbor_set.push_back(contact_set);
            contact_set.clear();
        }

        chain_set.push_back(neighbor_set);
    }

    global_contacts.push_back(chain_set);
    chain_set.clear();
    // exit(0);

    std::cout << "Original Contacts obtained!" << std::endl;
    std::cout << global_contacts.size() << std::endl;

    // Print some of the original contacts.
    // std::cout << "Printing global_contacts, 0th frame: " << std::endl;
    // print_global_contacts(global_contacts);
    // print_global_contacts_count(global_contacts);
    // exit(0);
#endif // MTMAP2_BEFORE


#ifdef PHIPSI_B // PHIPSI Beginning section.
    std::cout << "Getting PHI / PSI Angles!" << std::endl;
    // aa_ref: crystal structure.
    // aa_zero: dcd-0
    // aa_later: dcd-time-later.

    vvAtoms vec_dihedrals; // vector<vector<Atom>>


    Atom *aa_temp_backbone;
    try
    {
        aa_temp_backbone = new Atom[num_atoms];
    }
    catch (std::bad_alloc xa)
    {
        std::cout << "Allocation Failure\n";
        exit(1);
    }

    // Select backbone
    num_select = system_select_atomtype(aa_ref,"backbone",num_atoms,aa_temp_backbone);
    int jk;
    // for(int j=0; j<num_select; j+=4)
    for(int j=0; j<num_select; j+=3)
    {
        vAtoms backbone(3); // vector<Atom>
        for(int k=0; k<3; k++)
        {
            jk = j + k;
            Atom a1 = aa_temp_backbone[jk];

            // std::cout << "atomtype: " << a1.atomtype << std::endl;
            if(a1.atomtype == "N")
            {
                backbone[0] = a1;
            }
            else if(a1.atomtype == "CA")
            {
                backbone[1] = a1;
            }
            else if(a1.atomtype == "C")
            {
                backbone[2] = a1;
            }
            // else if(a1.atomtype == "O")
            // {
            //     backbone[3] = a1;
            // }
        }
        vec_dihedrals.push_back(backbone);
        backbone.clear();
    }

    // // Select CA
    // num_select = system_select_atomtype(aa_ref,"CA",num_atoms,aa_temp_backbone);
    // // aa_temp_backbone[0].print_coords();
    // for(int j=0; j<num_select; j++)
    // {
    //     Atom a1 = aa_temp_backbone[j]; // Copy Constructor
    //     aa_backbone.push_back(a1);
    // }
    // num_select = system_select_atomtype(aa_ref,"C",num_atoms,aa_temp_backbone);
    // for(int j=0; j<num_select; j++)
    // {
    //     Atom a1 = aa_temp_backbone[j]; // Copy Constructor
    //     aa_backbone.push_back(a1);
    // }
    // num_select = system_select_atomtype(aa_ref,"O",num_atoms,aa_temp_backbone);
    // for(int j=0; j<num_select; j++)
    // {
    //     Atom a1 = aa_temp_backbone[j]; // Copy Constructor
    //     aa_backbone.push_back(a1);
    // }
    // num_select = system_select_atomtype(aa_ref,"N",num_atoms,aa_temp_backbone);
    // for(int j=0; j<num_select; j++)
    // {
    //     Atom a1 = aa_temp_backbone[j]; // Copy Constructor
    //     aa_backbone.push_back(a1);
    // }

    // Not needed .. cleared.
    // std::cout << "backbone_size: " << backbone.size() << std::endl;

    Global_PhiPsi global_phipsi;
    PhiPsi local_phipsi;
    local_phipsi = compute_phipsi(vec_dihedrals);
    global_phipsi.push_back(local_phipsi);
    local_phipsi.clear();


    for(auto p: global_phipsi)
    {
        for(auto ph: p)
        {
            std::cout << ph.first << " \t " << ph.second << std::endl;
        }
    }

#endif // PHIPSI_B End.


    /* ---------------------------------------------------------
       Analysis Before. Finish.
       --------------------------------------------------------- */

#endif // multi-dcd



#ifdef DCDREAD
    /* ---------------------------------------------------------
       Step 3.0 DCD Read. Preload.
       --------------------------------------------------------- */
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
       Step 3.2 DCD Read.
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


    /* ---------------------------------------------------------
       Step 3.3 Open DCD.
       --------------------------------------------------------- */
    // 2. to read a dcd.
    printf("----->  READING DCD  <-----\n");

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
    // printf("  %d atoms, %d frames, size: %6.1fMB\n", natoms, dcd->nsets, sizeMB);


    /* ---------------------------------------------------------
       Step 3.3 Load DCD.
       --------------------------------------------------------- */
    // Load first frame into allatoms_0. (0-position)
    frame_position = 1;
    advance_dcd(dcd->nsets,0,dcd,natoms,&timestep); // 1st advance. 1-vmd
    allatoms_0 = load_dcd_to_atoms(dcd,allatoms_0);

    // Load 2nd frame into allatoms. (time-later)
    frame_position = 2;
    advance_dcd(dcd->nsets,0,dcd,natoms,&timestep); // 2nd. 2-vmd
    allatoms = load_dcd_to_atoms(dcd,allatoms);

    printf("frame_position: %d\n",frame_position);


    /* ---------------------------------------------------------
       Step 3.4 Advance DCD.
       --------------------------------------------------------- */
    // Advancing Rules.
    // ----------------
    // example. step size -> 5.
    // int advance_size = atoi(argv[3]) - 1;
    int advance_size = step - 1; // 0 counts, so advance_size of 4, advances by 5.

    // stop. Check that the stopping frame isn't beyond the actual frame count.
    if(stop > dcd->nsets)
    {
        stop = dcd->nsets;
        std::cout << "Use stop value: " << stop << std::endl;
    }


    /* ---------------------------------------------------------
       Step 3.5 Fast Forward DCD
       --------------------------------------------------------- */
    for (int nset1=2; nset1<start; nset1 += 1 )
    {
        frame_position += advance_dcd(dcd->nsets,0,dcd,natoms,&timestep);
        allatoms = load_dcd_to_atoms(dcd,allatoms);
        debug("forwarding --> frame_position: %d\n",frame_position);
    }

    // Get initial starting point.
    printf("FF: -->\nframe_position: %d\n",frame_position);


    /* ---------------------------------------------------------
       Step 3.6 Loop through DCD, evaluate. (do loop)
       --------------------------------------------------------- */
    int nset2;
    nset2 = frame_position;
    do {
#endif //DCDREAD


#if defined (DCDREAD) || defined (DCD_WRITE_B) || defined (DCD_WRITE) || defined (DCD_WRITE_E)
        /* ---------------------------------------------------------
           DURING DCD
           LOAD DCD
           Analysis Evaluation.
           --------------------------------------------------------- */


        /* ---------------------------------------------------------
           Check that the points are evolving in time.
           --------------------------------------------------------- */
#ifndef NDEBUG // Double Negative Define: Not No-debugging.
        if(frame_position <=10)
        {
            std::cout << "Allatoms-Ref:" << std::endl;
            for(int i=0; i<atomcheck; i++)
            {
                allatoms_ref[i].print_Coords();
            }

            std::cout << "Allatoms-t:" << std::endl;
            for(int i=0; i<atomcheck; i++)
            {
                allatoms[i].print_Coords();
            }
        }
#endif // NDEBUG


        /* ---------------------------------------------------------
           Analysis During:
           --------------------------------------------------------- */
#ifdef MTMAP2_DURING
        std::cout << "MTMAP2: Evaluating contacts by sector." << std::endl;
        std::cout << "Global Contacts: " << std::endl;
        std::cout << global_contacts[0].size() << std::endl;
        // SetContacts contact_set;
        // SetNeighbors neighbor_set;
        // SetChains chain_set;
        // SetGlobalContacts global_contacts;

        // Clear as precaution:
        contact_set.clear();
        neighbor_set.clear();
        chain_set.clear();

        // c: chain (but actually 156 dimers in 312 chains for the MT)
        // n: 9 neighboring interaction types
        int it_c, it_n;
        it_c = it_n = 0;

        for(auto c: global_contacts[0])
        {
            it_n = 0; // 0-8, the 9 interaction types: alpha, beta, alpha-beta,
            // alpha-s,w,e..
            // std::cout << "c: " << c.size() << std::endl; // ~ 9

            for(auto n: c)
            {
                // aa_later, allatoms_chain, allatoms, sel_all
                contact_set = get_contacts_for_chain_later(allatoms,
                                                           8.0,2.0,
                                                           global_contacts[0][it_c][it_n]);

                // std::cout << contact_set.size() << std::endl;
                neighbor_set.push_back(contact_set);
                contact_set.clear();
                it_n += 1;
            }
            chain_set.push_back(neighbor_set);
            neighbor_set.clear();
            it_c += 1;
        }
        global_contacts.push_back(chain_set);


        // for(auto c: mt_matrix)
        // {
        //     // neighbor_set.clear();
        //     std::cout << " "
        //               << c[0] << " "
        //               << c[1] << " "
        //               << c[2] << " "
        //               << c[3] << " "
        //               << c[4] << " "
        //               << c[5] << " "
        //               << c[6] << " "
        //               << c[7] << " "
        //               << std::endl;

        //     for(int it_n=0; it_n < c.size(); it_n++)
        //     {
        //         if(c[it_n] >= 0)
        //         {
        //             contact_set = get_contacts_for_chain_later(aa_later,
        //                                                        8.0,2.0,
        //                                                        global_contacts[0][c[it_n]][it_n]);
        //             neighbor_set.push_back(contact_set);
        //             contact_set.clear();
        //         }
        //         else
        //         {
        //             contact_set.clear();
        //             neighbor_set.push_back(contact_set);
        //         }

        //     }

        // Checkpoint.
        // std::cout << global_contacts[0] // frame(0) --> 156 --> 8
        // std::cout << global_contacts[0][0][4].size() << std::endl; // frame-156-8-sc
        // std::cout << global_contacts[0][1][4].size() << std::endl; // frame-156-8-sc
        // std::cout << global_contacts[0][2][1].size() << std::endl; // frame-156-8-sc
        // exit(0);
#endif // MTMAP2_DURING


#ifdef PHIPSI_M // Begin.
        std::cout << "Getting the phi / psi angles." << std::endl;

        // Reload.
        load_dcd_to_atoms(dcd,aa_temp_backbone);
        vec_dihedrals.clear();
        for(int j=0; j<num_select; j+=3)
        {
            vAtoms backbone(3); // vector<Atom>

            for(int k=0; k<3; k++)
            {
                jk = j + k;
                Atom a1 = aa_temp_backbone[jk];

                // std::cout << "atomtype: " << a1.atomtype << std::endl;
                if(a1.atomtype == "N")
                {
                    backbone[0] = a1;
                }
                else if(a1.atomtype == "CA")
                {
                    backbone[1] = a1;
                }
                else if(a1.atomtype == "C")
                {
                    backbone[2] = a1;
                }
                // else if(a1.atomtype == "O")
                // {
                //     backbone[3] = a1;
                // }
            }
            vec_dihedrals.push_back(backbone);
            backbone.clear();
        }

        local_phipsi = compute_phipsi(vec_dihedrals);
        global_phipsi.push_back(local_phipsi);
        local_phipsi.clear();

#endif // PHIPSI End.


        /* ---------------------------------------------------------
           Step 4. Analysis During. Finish.
           DCD LOAD complete.
           No more analysis. (unless it's just printing.)
           --------------------------------------------------------- */
#endif // multi-dcd


#ifdef DCD_WRITE
        /* ---------------------------------------------------------
           Step 3. DCD Write. Finish.
           --------------------------------------------------------- */
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
        // load_atom_to_timestep(&timestep_w,aa_later);
        load_atom_to_timestep(&timestep_w,allatoms);

#ifdef DCD_WRITE_UNMOD
        // Write the DCD read in.
        write_timestep(vw,&timestep);
#else
        // Write modified coordinates.
        write_timestep(vw,&timestep_w);
#endif // DCD_WRITE_UNMOD

#endif // DCD_WRITE


#ifdef DCDREAD
        /* ---------------------------------------------------------
           Step 3.6 DCD READ
           --------------------------------------------------------- */
        debug("current: %d\n",nset2);
        printf("frame: --> %d <-- was evaluated.\n",frame_position);

        if (nset2 + advance_size + 1 <= stop)
        {
            frame_position += advance_dcd(dcd->nsets,advance_size,dcd,natoms,&timestep);
            printf("frame: --> %d <-- loaded.\n",frame_position);

            // THIS ONE
            allatoms = load_dcd_to_atoms(dcd,allatoms);
        }
        nset2 += advance_size + 1;
    } while (nset2<=stop);

    debug("\n..closing dcd..\n");
    close_file_read(v);

    printf("\n----->  READING DCD completed!  <-----\n");
    printf("\t\tThe maximum possible frame_position was: %d\n",stop);
    printf("\t\tThe last frame evaluated was: %d\n",frame_position);
#endif //DCDREAD END



#if defined (DCDREAD) || defined (DCD_WRITE_B) || defined (DCD_WRITE) || defined (DCD_WRITE_E)
    /* ---------------------------------------------------------
       Analysis After. Start.
       --------------------------------------------------------- */

#ifdef MTMAP2_AFTER
    std::cout << "MTMAP2: Contacts by sector complete." << std::endl;

    // Print Analysis of Contacts File.
    output_global_contacts(global_contacts);

    // Print Analysis of Contacts by Subdomain.
    output_global_contacts_by_subdomain(global_contacts);

    // Notes:
    // fp_contacts = fopen("emol_mtcontacts_by_subdomain.dat", "w+");
    // fp_contacts3 = fopen("emol_mtcontacts_by_subdomain3.dat", "w+");
    // fp_contacts3n = fopen("emol_mtcontacts_by_subdomain3n.dat", "w+");
#endif // MTMAP2_AFTER


#ifdef PHIPSI_E // PHIPSI Final section.
    std::cout << "Getting phi/psi angles completed." << std::endl;
    delete [] aa_temp_backbone;

    int n = 0;
    for(auto p: global_phipsi)
    {
        n += 1;
        std::cout << "Frame: " << n << std::endl;
        for(auto ph: p)
        {
            std::cout << ph.first << " \t " << ph.second << std::endl;
        }
    }
#endif // PHIPSI_E End.

#endif // multi-dcd


#ifdef DCD_WRITE_E
    // open_dcd_write(fn_dcd_write,"dcd",natoms);
    // static void close_file_write(void *v) {
    close_file_write(vw);
#endif // DCD_WRITE_E



    /* ---------------------------------------------------------
       Step 6. Delete Malloc Systems.
       --------------------------------------------------------- */
    // delete [] a_initial;
    // delete [] aa_ref;
    // delete [] aa_zero;
    // delete [] aa_later;
    // delete [] aa_sel;


    /* ---------------------------------------------------------
       The End.
       --------------------------------------------------------- */
    std::cout << "\nClosing stdin,stdout,stderr.." << std::endl;
    std::cout << "Emolsion Out." << std::endl;
    fclose(stdin);
    fclose(stdout);
    fclose(stderr);
    return 0;
}
