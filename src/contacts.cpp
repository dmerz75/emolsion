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
#include <cmath>


/* ---------------------------------------------------------
   headers
   --------------------------------------------------------- */
#include "debug.h"
#include "contacts.hpp"
#include "md.h"
#include "system.hpp"
// #include "dcd.h"
// #include "dcdio.h"

/* ---------------------------------------------------------
   functions
   --------------------------------------------------------- */
// void ReadPDBfile (PDBfile *pdbfile,char filename[40],System sys)
void get_contacts(Atom *a1,Atom *a2,char dcdfilename[40],int natoms)
{
    debug("welcome to get_contacts. 2 (inter)\n");
    std::cout << "DCD: " << dcdfilename << std::endl;


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


// void get_map_of_mtneighbors(std::vector<std::vectorAtom> chain_ref,std::vector<std::vector<int>> matrix,
//                             std::)
// std::vector<std::vector<int>> get_map_of_mtneighbors(std::vector<std::vector <Atom>> chain_ref,std::vector<std::vector<int>> matrix,
//                             std::vector<std::pair<int,int>> dimers)
std::vector<std::vector<int>> get_map_of_mtneighbors(std::vector<std::vector <Atom>> chain_ref,
                                                     std::vector<std::pair<int,int>> dimers)
{
    // printf("Welcome to get_map_of_mtneighbors!\n");
    // std::cout << matrix.size() << std::endl;
    std::vector<std::vector<int>> matrix(dimers.size(), std::vector<int>(8,-1));

    int ic;
    ic = 0;

    Vector centroid;
    // centroid.print_Vector();
    std::vector<Vector> centroids;


    for(auto ch: chain_ref)
    {
        // std::cout << ch.size() << std::endl;
        centroid = get_centroid(ch);
        // centroid.print_Vector();
        centroids.push_back(centroid);


        // std::cout << "Return_Centroid: "
        //           << centroid[0]
        //           << " "
        //           << centroid[1]
        //           << " "
        //           << centroid[2]
        //           << std::endl;


        ic += 1;
    }
    // std::cout << "# of centroids: " << centroids.size() << std::endl;


    int id;
    id = 0;

    Vector vdist;
    Vector vdist_n;
    float vdist_mag;
    vdist_mag = 0.0;


    std::vector<int> monomers;

    for(auto d: dimers)
    {
        // std::cout << d.first
        //           << " "
        //           << d.second
        //           << " "
        //           << std::endl;

        monomers.push_back(d.first);
        monomers.push_back(d.second);

        // centroids[d.first].print_Vector();
        // centroids[d.second].print_Vector();

        // vdist = get_vector(centroids[d.first],centroids[d.second]);
        // vdist.print_Vector();
        // vdist_mag = magnitude(vdist);
        // std::cout << "Magnitude: " << vdist_mag << std::endl;

        matrix[id][0] = d.first;
        matrix[id][1] = d.second;

        id += 1;
    }


    // 83.0 Angstroms;
    std::vector<std::vector<int>> chain_candidates;
    std::vector<int> candidates;

    for(auto m1: monomers)
    {
        // std::cout << "monomer: " << m1 << std::endl;

        for(auto m2: monomers)
        {
            if(m2 == m1)
            {
                continue;
            }
            else
            {

                // std::cout << "comparing: m2-m1 " << m2 << " <-> "<< m1 << std::endl;
                vdist = get_vector(centroids[m1],centroids[m2]);
                // vdist.print_Vector();
                vdist_mag = magnitude(vdist);
                // std::cout << "Magnitude: " << vdist_mag << std::endl;
                if (vdist_mag < 83.0)
                {
                    candidates.push_back(m2);
                }
            }
        }
        chain_candidates.push_back(candidates);
        candidates.clear();
    }
    // std::cout << "Candidates acquired." << std::endl;


    int ican,pdim;
    ican = pdim = 0;

    // Axis of dimer.
    Vector avec;
    Vector avec_n;
    float avec_mag;
    avec_mag = 0.0;

    // SINE and COSINE
    double dsin;
    double dcos;


    for(auto cc: chain_candidates)
    {
        // std::cout << "monomer: " << ican << std::endl;

        if(ican % 2 == 0)
        {
            pdim = ican / 2;
        }
        else
        {
            pdim = (ican - 1) / 2;
        }

        for(auto can: cc)
        {
            avec = get_vector(centroids[dimers[pdim].second],centroids[dimers[pdim].first]);
            avec_n = normalize(avec);
            avec_mag = magnitude(avec);

            // std::cout << dimers[pdim].first
            //           << "-"
            //           << dimers[pdim].second
            //           << "   "
            //           << avec_mag
            //           << " "
            //           << std::endl;




            vdist = get_vector(centroids[ican],centroids[can]);
            vdist_n = normalize(vdist);
            vdist_mag = magnitude(vdist);
            // std::cout << can << " " << vdist_mag << " " << std::endl;


            dsin = get_sintheta(avec_n,vdist_n);
            dcos = get_costheta(avec_n,vdist_n);

            // std::cout << "sin: " << dsin << std::endl;
            // std::cout << "cos: " << dcos << std::endl;

            // std::cout << std::endl;


            //     West
            //     4  7
            //  2  0  1  5
            //     3  6
            //     East

            // find intra-dimer.
            if((dsin < 0.1) && (std::abs(dcos) > 0.95))
            {
                if((can != dimers[pdim].first) && (can != dimers[pdim].second))
                {
                    if(ican % 2 == 0)
                    {
                        // std::cout << "same-pf-south: " << can << std::endl;
                        matrix[pdim][2] = can;
                        continue;
                    }
                    else
                    {
                        // std::cout << "same-pf-north: " << can << std::endl;
                        matrix[pdim][5] = can;
                        continue;
                    }
                }
                // else
                // {
                //     std::cout << "intra-dimer: "<< can << std::endl;
                // }
            }
            // else if ((std::abs(dcos) > 0.3) && (std::abs(dcos) < 0.6))
            // {
            //     if(dsin > 0.0)
            //     {
            //         std::cout << "45-east-west^ " << std::endl;
            //     }
            //     else // dsin < 0.0
            //     {
            //         std::cout << "45-east-west^ " << std::endl;
            //     }
            // }
            else if ((dsin > 0.95) && (std::abs(dcos) < 0.28))
            {
                if(dcos > 0)
                {
                    // std::cout << "WEST " << std::endl;
                    if(ican % 2 == 0)
                    {
                        matrix[pdim][4] = can;
                    }
                    else
                    {
                        matrix[pdim][7] = can;
                    }
                }
                else
                {
                    // std::cout << "EAST " << std::endl;
                    if(ican % 2 == 0)
                    {
                        matrix[pdim][3] = can;
                    }
                    else
                    {
                        matrix[pdim][6] = can;
                    }
                }
            }
        } // monomer-candidates

        // std::cout << std::endl;

        ican += 1;

    } // monomers.

    return matrix;
}
