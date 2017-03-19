// dcd.cpp

/* ---------------------------------------------------------
   standard libraries
   --------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h> // strtod?, stod


// namespaces
// using namespace

// extern
extern int endoffile;

/* ---------------------------------------------------------
   headers
   --------------------------------------------------------- */
// #include "readfile.h"
// #include "chain.h"
#include "md.h"
#include "dcdio.h"
#include "system.hpp"

/* ---------------------------------------------------------
   functions
   --------------------------------------------------------- */
int advance_dcd(int numframes,int frame,dcdhandle *v,int natoms,molfile_timestep_t *timestep) {
    for (int i=0; i<numframes; i++) {
        int rc = read_next_timestep(v, natoms, timestep);
        if (rc) {
            fprintf(stderr, "error in read_next_timestep on frame %d\n", i);
            break;
            // exit(1);
        }

        if ( i == frame) {
            break;
        }
    }
    return frame + 1;
    // return frame;
}
void load_chain_to_timestep(System *chain,int num_chains,const molfile_timestep_t *ts)
{

    // debug("the number of chains to load: %d\n",num_chains);
    // debug("maximum_index: %d\n",chain[num_chains-1].findex);
    // debug("chain_count: %d\n",num_chains);

    // Works! - not necessary
    // dcdhandle *dcd1;
    // dcd1 = (dcdhandle *)malloc(sizeof(dcdhandle));
    // memset(dcd1, 0, sizeof(dcdhandle));
    // dcd1->x = (float *)malloc(natoms * sizeof(float));
    // dcd1->y = (float *)malloc(natoms * sizeof(float));
    // dcd1->z = (float *)malloc(natoms * sizeof(float));

    // int x,c,sindex,findex,total,natoms;
    int natoms,count,count1;
    natoms = count = count1 = 0;

    for(int h=0; h<num_chains; h++){

        // natoms += chain[h].num_atoms_ca;
        natoms += chain[h].num_atoms;

        // for(int i=0; i<chain[h].num_atoms_ca; i++){
        for(int i=0; i<chain[h].num_atoms; i++){
            // ts->coords[count] = (float)chain[h].pos[i].x;
            // ts->coords[count+1] = (float)chain[h].pos[i].y;
            // ts->coords[count+2] = (float)chain[h].pos[i].z;

            // FIX

            // ts->coords[count] = (float)chain[h].pos[i].x;
            // ts->coords[count+1] = (float)chain[h].pos[i].y;
            // ts->coords[count+2] = (float)chain[h].pos[i].z;

            count += 3;
            count1 ++;
        }
    }
    debug("number of atoms: %d\n",natoms);
    debug("total coordinates written: %d  %d\n",count1,count);

    // CHECK TS
    // for(int i=0; i<natoms; i++){
    //     printf("i: %d   %f %f %f\n",i*3,ts->coords[i*3],ts->coords[i*3+1],ts->coords[i*3+2]);
    // }

    return;
}

void load_atom_to_timestep(const molfile_timestep_t *ts,Atom *aa)
{
    /* ---------------------------------------------------------
       For DCD writing.
       --------------------------------------------------------- */
    // int natoms,count,count1;
    // natoms = count = count1 = 0;
    int count;
    count = 0;

    for (int i=0; i<aa[0].num_atoms; i++)
    {
        // ts->coords[count] = (float)chain[h].pos[i].x;
        // ts->coords[count+1] = (float)chain[h].pos[i].y;
        // ts->coords[count+2] = (float)chain[h].pos[i].z;

        // FIX

        ts->coords[count] = (float)aa[i].x;
        ts->coords[count+1] = (float)aa[i].y;
        ts->coords[count+2] = (float)aa[i].z;

        count += 3;
        // count1 ++;
    }
    debug("total coordinates written: %d\n",count);

    // CHECK TS
    // for(int i=0; i<natoms; i++){
    //     printf("i: %d   %f %f %f\n",i*3,ts->coords[i*3],ts->coords[i*3+1],ts->coords[i*3+2]);
    // }
    return;
}

void load_dcd_to_chain(dcdhandle *dcd,Chain *chain,int num_chains)
{

    for (int i=0; i<num_chains; i++) {

        // for ( int j=0; j<chain[i].num_atoms_ca; j++ ) {
        for ( int j=0; j<chain[i].num_atoms; j++ ) {

            // printf("%d ",chain[i].indices[j]);
            // chain[i].pos[j].x = dcd->x[chain[i].indices[j]];
            // chain[i].pos[j].y = dcd->y[chain[i].indices[j]];
            // chain[i].pos[j].z = dcd->z[chain[i].indices[j]];

            // FIX
            // chain[i].pos[j].x = dcd->x[chain[i].indices[j]];
            // chain[i].pos[j].y = dcd->y[chain[i].indices[j]];
            // chain[i].pos[j].z = dcd->z[chain[i].indices[j]];
        }
    }
}

void load_dcd_to_atoms(dcdhandle *dcd,Atom *aa)
{
    /* ---------------------------------------------------------
       for DCD reading.
       natoms: total atoms from DCD
       --------------------------------------------------------- */
    debug("coordinate loading..\n");
    debug("this many: %d\n",aa[0].num_atoms);
    // debug("this many: %d\n",aa[0]);

    // for the selection:
    for (int i=0; i<aa[0].num_atoms; i++)
    {
        // printf("index: %d\n",aa[i].index);
        aa[i].x = dcd->x[aa[i].index];
        aa[i].y = dcd->y[aa[i].index];
        aa[i].z = dcd->z[aa[i].index];
    }


    // 2nd rendition.
    // for (int i=0; i<aa[0].num_atoms; i++)
    // {
    //         // printf("%d ",chain[i].indices[j]);
    //         // chain[i].pos[j].x = dcd->x[chain[i].indices[j]];
    //         // chain[i].pos[j].y = dcd->y[chain[i].indices[j]];
    //         // chain[i].pos[j].z = dcd->z[chain[i].indices[j]];

    //     // debug("this far: %d\n",i);
    //     // debug("%f %f %f\n",dcd->x[i],dcd->y[i],dcd->z[i]);
    //     // FIX
    //     aa[i].x = dcd->x[i];
    //     aa[i].y = dcd->y[i];
    //     aa[i].z = dcd->z[i];
    //         // chain[i].pos[j].x = dcd->x[chain[i].indices[j]];
    //         // chain[i].pos[j].y = dcd->y[chain[i].indices[j]];
    //         // chain[i].pos[j].z = dcd->z[chain[i].indices[j]];
    // }

}

std::vector<Atom> load_dcd_to_atoms(dcdhandle *dcd,std::vector<Atom> atoms)
// void load_dcd_to_atoms(dcdhandle *dcd,std::vector<Atom> chain)
{
    // Vector of Atoms...
    std::vector<Atom> movedatoms;
    Atom ma;

    /* ---------------------------------------------------------
       for DCD reading.
       natoms: total atoms from DCD
       --------------------------------------------------------- */
    // debug("coordinate loading for Vector of Atoms..\n");

    for(auto aa: atoms)
    {
        // aa.print_Atom();
        aa.x = dcd->x[aa.index];
        aa.y = dcd->y[aa.index];
        aa.z = dcd->z[aa.index];
        ma = aa;
        // aa.print_Atom();
        // ma.print_Atom();
        // std::cout << std::endl;

        // ma.x = dcd->x[aa.index];
        // ma.y = dcd->y[aa.index];
        // ma.z = dcd->z[aa.index];
        movedatoms.push_back(ma);
    }


    // for the selection:
    // for (int i=0; i<chain.size(); i++)
    // {
    //     // printf("index: %d\n",chain[i].index);
    //     chain[i].x = dcd->x[chain[i].index];
    //     chain[i].y = dcd->y[chain[i].index];
    //     chain[i].z = dcd->z[chain[i].index];
    // }


    return movedatoms;
}

void load_dcd_to_atoms(dcdhandle *dcd,std::vector<Atom *> atoms)
{
    // Vector of Atoms...
    std::cout << "Moving the Atoms at Pointers." << std::endl;

    /* ---------------------------------------------------------
       for DCD reading.
       natoms: total atoms from DCD
       --------------------------------------------------------- */
    // debug("coordinate loading for Vector of Pointers to Atoms..\n");

    for(auto aa: atoms)
    {
        // aa.print_Atom();
        aa->x = dcd->x[aa->index];
        aa->y = dcd->y[aa->index];
        aa->z = dcd->z[aa->index];
    }
    return;
}


// void load_dcd_to_atoms(dcdhandle *dcd,Atom aa)
// {
//     /* ---------------------------------------------------------
//        for DCD reading.
//        natoms: total atoms from DCD
//        --------------------------------------------------------- */
//     debug("coordinate loading..\n");
//     debug("this many: %d\n",aa[0]->num_atoms);

//     // for the selection:
//     for (int i=0; i<aa[0]->num_atoms; i++)
//     {
//         // printf("index: %d\n",aa[i].index);
//         aa[i]->x = dcd->x[aa[i].index];
//         aa[i]->y = dcd->y[aa[i].index];
//         aa[i]->z = dcd->z[aa[i].index];
//     }


//     // 2nd rendition.
//     // for (int i=0; i<aa[0].num_atoms; i++)
//     // {
//     //         // printf("%d ",chain[i].indices[j]);
//     //         // chain[i].pos[j].x = dcd->x[chain[i].indices[j]];
//     //         // chain[i].pos[j].y = dcd->y[chain[i].indices[j]];
//     //         // chain[i].pos[j].z = dcd->z[chain[i].indices[j]];

//     //     // debug("this far: %d\n",i);
//     //     // debug("%f %f %f\n",dcd->x[i],dcd->y[i],dcd->z[i]);
//     //     // FIX
//     //     aa[i].x = dcd->x[i];
//     //     aa[i].y = dcd->y[i];
//     //     aa[i].z = dcd->z[i];
//     //         // chain[i].pos[j].x = dcd->x[chain[i].indices[j]];
//     //         // chain[i].pos[j].y = dcd->y[chain[i].indices[j]];
//     //         // chain[i].pos[j].z = dcd->z[chain[i].indices[j]];
//     // }

// }
