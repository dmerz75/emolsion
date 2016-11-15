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
#include "dcd.h"

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

    if (argc < 3) {
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
#ifdef DCDREAD
    int frame_position;
    frame_position = 0;
    int start,stop,step;

    // start,stop,step
    start = atoi(argv[3]);
    stop = atoi(argv[4]);
    step = atoi(argv[5]);



    /* ---------------------------------------------------------
       write dcd
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

    // get atom total for writing.
    for(int i=0; i<num_chains; i++){
        natoms_w += chain_ref[i].num_atoms_ca;
    }
    printf("for dcd writing >>>  <%d> atoms expected.\n",natoms_w);


    vw = open_dcd_write(fn_dcd_write,"dcd",natoms_w);
    if (!vw) {
        fprintf(stderr, "main) open_dcd_write failed for file %s\n", *fn_dcd_write);
        return 1;
    } else {
        printf("opened <%s> successfully!!\n\n",fn_dcd_write);
    }

    timestep_w.coords = (float *)malloc(3*sizeof(float)*natoms_w);

    // dcd = (dcdhandle *)v;
    // sizeMB = ((natoms * 3.0) * dcd->nsets * 4.0) / (1024.0 * 1024.0);
    // totalMB += sizeMB;
    // printf("main) file: %s\n", *argv);
    // printf("  %d atoms, %d frames, size: %6.1fMB\n", natoms, dcd->nsets, sizeMB);
    // timestep.coords = (float *)malloc(3*sizeof(float)*natoms);
#endif // DCD_WRITE_B


    /* ---------------------------------------------------------
       DCD Preface.
       --------------------------------------------------------- */
    // Start DCD
    molfile_timestep_t timestep;
    void *v;
    dcdhandle *dcd;
    // int i, natoms;
    int natoms;
    float sizeMB =0.0, totalMB = 0.0;
    double starttime, endtime, totaltime = 0.0;

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



    // int atoms_in_chain;
    // 2. to read a dcd.
    natoms = 0;
    v = open_dcd_read(argv[2], "dcd", &natoms);
    if (!v) {
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
    // END DCD
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
    printf("--------START HERE-------------\n");

    // dcd
    // 0: (pdb) | ref | chain_0 (from dcd) | chain_later
    // int advance_dcd(int numframes,int frame,dcdhandle *v,int natoms,molfile_timestep_t *timestep);

    // printf("ref: %f\n",chain_ref[0].pos[105].y);
    // printf("0-105: %f\n",chain_0[0].pos[105].y);

    // frame_position = advance_dcd(dcd->nsets,0,dcd,natoms,&timestep); // 1
    // load_dcd_to_chain(dcd,chain_0,num_chains);
    // printf("frame_position(0): %d\n",frame_position);
    // printf("0-105: %f\n",chain_0[0].pos[105].y);

    // frame_position = advance_dcd(dcd->nsets,1,dcd,natoms,&timestep); // 1
    // load_dcd_to_chain(dcd,chain_0,num_chains);
    // printf("frame_position: %d\n",frame_position);
    // printf("0-105: %f\n",chain_0[0].pos[105].y);

    // frame_position = advance_dcd(dcd->nsets,10,dcd,natoms,&timestep); // 1
    // load_dcd_to_chain(dcd,chain_0,num_chains);
    // printf("frame_position: %d\n",frame_position);
    // printf("0-105: %f\n",chain_0[0].pos[105].y);


    // printf("beginning new loop.\n");
    // int countd;

    // countd = 0;
    // for(int df=start; df<=stop; df+=step){
    //     countd +=1;
    //     // dummy = my_dcd_read(df);
    //     frame_position = advance_dcd(dcd->nsets,df,dcd,natoms,&timestep);
    //     printf("frame_position:--->  %d  <--- %d %d\n",frame_position,countd,df);
    //     load_dcd_to_chain(dcd,chain_later,num_chains);
    //     printf("0-105: %f\n",chain_later[0].pos[105].y);

    // }
    // exit(0);

    // nil.

    frame_position = 1;
    advance_dcd(dcd->nsets,0,dcd,natoms,&timestep); // 1st advance. 1-vmd

    // THIS ONE
    // load_dcd_to_chain(dcd,chain_0,num_chains);

    // printf("frame_position: %d\n",frame_position);
    // printf("ref-findex(%d): %f\n",chain_ref[0].findex,chain_ref[0].pos[chain_ref[0].findex].y);
    // printf("0-findex(%d): %f\n",chain_0[0].findex,chain_0[0].pos[chain_0[0].findex].y);
    // printf("later-findex(%d): %f\n",chain_later[0].findex,chain_later[0].pos[chain_later[0].findex].y);


    frame_position = 2;
    advance_dcd(dcd->nsets,0,dcd,natoms,&timestep); // 2nd. 2-vmd


    // THIS ONE
    // load_dcd_to_chain(dcd,chain_later,num_chains);


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

        if (nset2 + advance_size + 1 <= stop) {
            frame_position += advance_dcd(dcd->nsets,advance_size,dcd,natoms,&timestep);
            printf("frame: --> %d <-- loaded.\n",frame_position);

            // THIS ONE
            // load_dcd_to_chain(dcd,chain_later,num_chains);


            // nset2 += advance_size + 1;
        }
        nset2 += advance_size + 1;
    } while (nset2<=stop);

    debug("..closing dcd..\n");
    close_file_read(v);
    // END DCD
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

    // delete [] chain_ref;
    // delete [] chain;
    // delete [] chain_0;
    // delete [] chain_later;

    printf("DCDREAD complete.\n\tThe maximum possible frame_position was: %d\n",stop);
    printf("\tThe last frame evaluated was: %d\n",frame_position);
#endif //DCDREAD



#ifdef DCD_WRITE_E
    // open_dcd_write(fn_dcd_write,"dcd",natoms);
    // static void close_file_write(void *v) {
    close_file_write(vw);
#endif // DCD_WRITE_E



    /* ---------------------------------------------------------
       The End.
       --------------------------------------------------------- */
    std::cout << "\nclosing stdin,stdout,stderr.\n";
    fclose(stdin);
    fclose(stdout);
    fclose(stderr);
    return 0;
}
