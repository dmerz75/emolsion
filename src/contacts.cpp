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
// #include <iterator> // istream_iterator
#include <cmath>
#include <vector>
// #include "boost/tuple/tuple.hpp"
#include <iomanip> // setw
#include <fstream> // ifstream
#include <string.h> // string.
#include <boost/algorithm/string.hpp> // std::string, boost


/* ---------------------------------------------------------
   headers
   --------------------------------------------------------- */
#include "debug.h"
#include "contacts.hpp"
#include "md.h"
#include "system.hpp"
#include "microtubule.hpp"
#include "math_atom.hpp"
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

SetContacts get_contacts_for_chain(vAtoms aa,
                                   double cutoff,
                                   IndexGroup ig1,
                                   IndexGroup ig2)
{
    // SM, Single Molecule
    // INTERCHAIN, Not MT

    // std::cout << "Welcome to get_contacts_for_chain!" << std::endl;
    // std::cout << "Cutoff: " << cutoff << std::endl;
    // std::vector<boost::tuple<int,int,double>> contacts;
    SetContacts contacts;

    Vector a1v,a2v;
    double dist = 0.0;
    double eh = 0.0;

    int subd1 = -1; // subdomains 1 & 2
    int subd2 = -1;

    for(auto i: ig1)
    {
        // std::cout << i1 << " " << a.index << std::endl;
        // std::cout << a.index << std::endl;
        a1v.x = aa[i].x;
        a1v.y = aa[i].y;
        a1v.z = aa[i].z;

        for(auto j: ig2)
        {
            // Directive: Exclude an index being in contact with itself,
            // or its +/- 2 neighbor.
            // if((j - 2 >= i) and (j + 2 <= i))
            if((j >= i - 2) and (j <= i + 2))
            {
                continue;
            }
            if(j < i)
            {
                continue;
            }

            a2v.x = aa[j].x;
            a2v.y = aa[j].y;
            a2v.z = aa[j].z;
            dist = distance(a1v,a2v);

            if(dist <= cutoff)
            {
                contacts.push_back(boost::tuple<int,int,double,int,int,double>(
                                       i,
                                       j,
                                       dist,
                                       subd1,
                                       subd2,
                                       eh));
            }
        }
    }
    return contacts;
}

SetContacts get_contacts_for_chain(vAtoms aa,
                                   double cutoff,
                                   MtIndexMap map,
                                   IndexGroup ig1,
                                   IndexGroup ig2,
                                   int cid1,int cid2)
{
    // INTERCHAIN, w/MT
    // std::cout << "Welcome to get_contacts_for_chain! MT, inter-chain." << std::endl;
    // std::cout << "Cutoff: " << cutoff << std::endl;
    // std::vector<boost::tuple<int,int,double>> contacts;
    SetContacts contacts;
    // std::cout << "chain1: " << chain1.size() << std::endl;
    // std::cout << "chain2: " << chain2.size() << std::endl;

    Vector a1v,a2v;
    double dist = 0.0;
    double eh = 0.0;

    int subd1 = -1; // subdomains 1 & 2
    int subd2 = -1;

    for(auto i: ig1)
    {
        // a1v.x = a1->x;
        // a1v.y = a1->y;
        // a1v.z = a1->z;
        a1v.x = aa[i].x;
        a1v.y = aa[i].y;
        a1v.z = aa[i].z;

        for(auto j: ig2)
        {
            // Directive: Exclude an index being in contact with itself,
            // or its +/- 2 neighbor.
            // std::cout << i << "---" << j << std::endl;
            if((j >= i - 2) and (j <= i + 2))
            {
                continue;
            }
            if(j < i)
            {
                continue;
            }

            a2v.x = aa[j].x;
            a2v.y = aa[j].y;
            a2v.z = aa[j].z;
            dist = distance(a1v,a2v);
            // a.print_coords();
            // a2.print_coords();

            if(dist <= cutoff)
            {
                // contacts.push_back(boost::tuple<int,int,int,double>(a1.index,a2.index,1,dist));
                // total_contacts += 1;
                if((i >= map[cid1]["index"]) and (i <= map[cid1]["Nterm2"]))
                {
                    subd1 = 0;
                }
                else if((i >= map[cid1]["Mterm1"]) and (i <= map[cid1]["Mterm2"]))
                {
                    subd1 = 1;
                }
                else if((i >= map[cid1]["Cterm1"]) and (i <= map[cid1]["findex"]))
                {
                    subd1 = 2;
                }

                if((j >= map[cid2]["index"]) and (j <= map[cid2]["Nterm2"]))
                {
                    subd2 = 0;
                }
                else if((j >= map[cid2]["Mterm1"]) and (j <= map[cid2]["Mterm2"]))
                {
                    subd2 = 1;
                }
                else if((j >= map[cid2]["Cterm1"]) and (j <= map[cid2]["findex"]))
                {
                    subd2 = 2;
                }

                if((subd1 == -1) or (subd2 == -1))
                {
                    // std::cout << cid1 << " " << cid2 << " " << std::endl;
                    // std::cout << a1.index << " " << a2.index << " " << std::endl;
                    // std::cout << map[cid1] << " " << map[cid2] << " " << std::endl;

                    // std::cout << "chaintype: " << map[cid1]["chaintype"] << " \n"
                    //           << "index: " << map[cid1]["index"] << " \n"
                    //           << "Nterm2: " << map[cid1]["Nterm2"] << " \n"
                    //           << "Mterm1: " << map[cid1]["Mterm1"] << " \n"
                    //           << "Mterm2: " << map[cid1]["Mterm2"] << " \n"
                    //           << "Cterm1: " << map[cid1]["Cterm1"] << " \n"
                    //           << "findex: " << map[cid1]["findex"] << " \n"
                    //           << std::endl;

                    // std::cout << "chaintype: " << map[cid2]["chaintype"] << " \n"
                    //           << "index: " << map[cid2]["index"] << " \n"
                    //           << "Nterm2: " << map[cid2]["Nterm2"] << " \n"
                    //           << "Mterm1: " << map[cid2]["Mterm1"] << " \n"
                    //           << "Mterm2: " << map[cid2]["Mterm2"] << " \n"
                    //           << "Cterm1: " << map[cid2]["Cterm1"] << " \n"
                    //           << "findex: " << map[cid2]["findex"] << " \n"
                    //           << std::endl;



                    // std::cout << subd1 << " " << subd2 << " " << std::endl;

                    exit(0);
                }

                contacts.push_back(boost::tuple<int,int,double,int,int,double>
                                   (i,
                                    j,
                                    dist,
                                    subd1,
                                    subd2,
                                    eh));
            }
        }
    }
    // std::cout << "contacts_size: " << contacts.size() << std::endl;
    // std::cout << "contacts_counted: " << total_contacts << std::endl;
    return contacts;
}

SetContacts get_contacts_for_chain(vAtoms aa,
                                   double cutoff,
                                   MtIndexMap map,
                                   IndexGroup ig1,
                                   int cid1)
{
    // MT - intra
    // std::cout << "Welcome to get_contacts_for_chain! MT, inter-chain." << std::endl;
    // std::cout << "Cutoff: " << cutoff << std::endl;
    // std::vector<boost::tuple<int,int,double>> contacts;
    SetContacts contacts;
    // std::cout << "chain1: " << chain1.size() << std::endl;
    // std::cout << "chain2: " << chain2.size() << std::endl;

    Vector a1v,a2v;
    double dist = 0.0;
    double eh = 0.0;

    int subd1 = -1; // subdomains 1 & 2
    int subd2 = -1;

    for(auto i: ig1)
    {
        a1v.x = aa[i].x;
        a1v.y = aa[i].y;
        a1v.z = aa[i].z;

        for(auto j: ig1)
        {
            // Directive: Exclude an index being in contact with itself,
            // or its +/- 2 neighbor.
            // if((j - 2 >= i) and (j + 2 <= i))
            if((j >= i - 2) and (j <= i + 2))
            {
                continue;
            }
            if(j < i)
            {
                continue;
            }

            a2v.x = aa[j].x;
            a2v.y = aa[j].y;
            a2v.z = aa[j].z;
            dist = distance(a1v,a2v);
            // a.print_coords();
            // a2.print_coords();

            if(dist <= cutoff)
            {
                // contacts.push_back(boost::tuple<int,int,int,double>(a1.index,a2.index,1,dist));
                // total_contacts += 1;
                if((i >= map[cid1]["index"]) and (i <= map[cid1]["Nterm2"]))
                {
                    subd1 = 0;
                }
                else if((i >= map[cid1]["Mterm1"]) and (i <= map[cid1]["Mterm2"]))
                {
                    subd1 = 1;
                }
                else if((i >= map[cid1]["Cterm1"]) and (i <= map[cid1]["findex"]))
                {
                    subd1 = 2;
                }

                if((j >= map[cid1]["index"]) and (j <= map[cid1]["Nterm2"]))
                {
                    subd2 = 0;
                }
                else if((j >= map[cid1]["Mterm1"]) and (j <= map[cid1]["Mterm2"]))
                {
                    subd2 = 1;
                }
                else if((j >= map[cid1]["Cterm1"]) and (j <= map[cid1]["findex"]))
                {
                    subd2 = 2;
                }

                if((subd1 == -1) or (subd2 == -1))
                {
                    // std::cout << cid1 << " " << cid2 << " " << std::endl;
                    // std::cout << a1.index << " " << a2.index << " " << std::endl;


                    // std::cout << map[cid1] << " " << map[cid2] << " " << std::endl;



                    // std::cout << "chaintype: " << map[cid1]["chaintype"] << " \n"
                    //           << "index: " << map[cid1]["index"] << " \n"
                    //           << "Nterm2: " << map[cid1]["Nterm2"] << " \n"
                    //           << "Mterm1: " << map[cid1]["Mterm1"] << " \n"
                    //           << "Mterm2: " << map[cid1]["Mterm2"] << " \n"
                    //           << "Cterm1: " << map[cid1]["Cterm1"] << " \n"
                    //           << "findex: " << map[cid1]["findex"] << " \n"
                    //           << std::endl;


                    // std::cout << "chaintype: " << map[cid2]["chaintype"] << " \n"
                    //           << "index: " << map[cid2]["index"] << " \n"
                    //           << "Nterm2: " << map[cid2]["Nterm2"] << " \n"
                    //           << "Mterm1: " << map[cid2]["Mterm1"] << " \n"
                    //           << "Mterm2: " << map[cid2]["Mterm2"] << " \n"
                    //           << "Cterm1: " << map[cid2]["Cterm1"] << " \n"
                    //           << "findex: " << map[cid2]["findex"] << " \n"
                    //           << std::endl;



                    // std::cout << subd1 << " " << subd2 << " " << std::endl;

                    exit(0);
                }

                contacts.push_back(boost::tuple<int,int,double,int,int,double>
                                   (i,
                                    j,
                                    dist,
                                    subd1,
                                    subd2,
                                    eh));
            }
        }
    }
    // std::cout << "contacts_size: " << contacts.size() << std::endl;
    // std::cout << "contacts_counted: " << total_contacts << std::endl;
    return contacts;
}




// SetContacts get_contacts_for_chain(vAtoms aa,
//                                    double cutoff)
// {
//     // SINGLE, INTRACHAIN, Not MT

//     // std::cout << "Welcome to get_contacts_for_chain! MT, intra." << std::endl;
//     // std::cout << "Cutoff: " << cutoff << std::endl;
//     // std::vector<boost::tuple<int,int,double>> contacts;
//     SetContacts contacts;
//     // std::cout << "aa: " << aa.size() << std::endl;
//     // std::cout << "chain2: " << chain2.size() << std::endl;


//     // int i1, i2;
//     // i1 = i2 = 0;
//     // int i = 0;

//     Vector a1v,a2v;
//     double dist = 0.0;

//     int subd1 = -1; // subdomains 1 & 2
//     int subd2 = -1;

//     // int total_contacts = 0;

//     for(auto a1: aa)
//     {
//         // std::cout << i1 << " " << a1->index << std::endl;
//         // std::cout << a1->index << std::endl;
//         a1v.x = a1->x;
//         a1v.y = a1->y;
//         a1v.z = a1->z;

//         for(auto a2: aa)
//         {
//             // Directive: Exclude an index being in contact with itself,
//             // or its +/- 2 neighbor.
//             // if((a2->index - 2 >= a1->index) and (a2->index + 2 <= a1->index))
//             // {
//             //     continue;
//             // }
//             if (a2->index - 2 <= a1->index)
//             {
//                 continue;
//             }


//             // std::cout << "\t" << i2 << "\t" << a2->index << std::endl;
//             // std::cout << "\t" << a2->index << std::endl;

//             a2v.x = a2->x;
//             a2v.y = a2->y;
//             a2v.z = a2->z;
//             dist = distance(a1v,a2v);

//             // a1->print_coords();
//             // a2->print_coords();

//             if(dist <= cutoff)
//             {
//                 // contacts.push_back(boost::tuple<int,int,int,double>(a1->index,a2->index,1,dist));
//                 // total_contacts += 1;
//                 contacts.push_back(boost::tuple<int,int,double,int,int>(
//                                        a1->index,
//                                        a2->index,
//                                        dist,
//                                        subd1,
//                                        subd2));
//                 // try
//                 // {
//                 //     contacts.push_back(boost::tuple<int,int,double>(a1->index,a2->index,dist));

//                 // }
//                 // catch (const std::bad_alloc &contacts)
//                 // {
//                 //     std::cout << "Allocation failed for single contact: " << contacts.what() << std::endl;
//                 // }
//             }
//         }
//     }
//     return contacts;
// }

// SetContacts get_contacts_for_chain(vpAtoms chain1,
//                                    double cutoff,
//                                    MtIndexMap map,
//                                    int cid)
// {
//     // SINGLE, INTRACHAIN, w/MT

//     // std::cout << "Welcome to get_contacts_for_chain!" << std::endl;
//     // std::cout << "Cutoff: " << cutoff << std::endl;
//     // std::vector<boost::tuple<int,int,double>> contacts;
//     SetContacts contacts;
//     // std::cout << "chain1: " << chain1.size() << std::endl;
//     // std::cout << "chain2: " << chain2.size() << std::endl;


//     // int i1, i2;
//     // i1 = i2 = 0;
//     // int i = 0;

//     Vector a1v,a2v;
//     double dist = 0.0;

//     int subd1 = -1; // subdomains 1 & 2
//     int subd2 = -1;

//     // int total_contacts = 0;

//     // int ch1 = -1;
//     // int ch2 = -1;

//     for(auto a1: chain1)
//     {
//         // ch1 += 1;
//         // ch2 = -1;

//         // std::cout << i1 << " " << a1->index << std::endl;
//         // std::cout << a1->index << std::endl;
//         a1v.x = a1->x;
//         a1v.y = a1->y;
//         a1v.z = a1->z;

//         for(auto a2: chain1)
//         {
//             // ch2 += 1;

//             // Directive: Exclude an index being in contact with itself,
//             // or its +/- 2 neighbor.
//             // if((a2->index - 2 >= a1->index) and (a2->index + 2 <= a1->index))
//             // {
//             //     continue;
//             // }
//             if (a2->index - 2 <= a1->index)
//             {
//                 continue;
//             }


//             // std::cout << "\t" << i2 << "\t" << a2.index << std::endl;
//             // std::cout << "\t" << a2->index << std::endl;



//             a2v.x = a2->x;
//             a2v.y = a2->y;
//             a2v.z = a2->z;

//             dist = distance(a1v,a2v);


//             // a1->print_coords();
//             // a2->print_coords();


//             if(dist <= cutoff)
//             {
//                 // contacts.push_back(boost::tuple<int,int,int,double>(a1->index,a2->index,1,dist));
//                 // total_contacts += 1;

//                 // std::cout << "chaintype: " << m["chaintype"] << " \n"
//                 //           << "index: " << m["index"] << " \n"
//                 //           << "Nterm2: " << m["Nterm2"] << " \n"
//                 //           << "Mterm1: " << m["Mterm1"] << " \n"
//                 //           << "Mterm2: " << m["Mterm2"] << " \n"
//                 //           << "Cterm1: " << m["Cterm1"] << " \n"
//                 //           << "findex: " << m["findex"] << " \n"

//                 if((a1->index >= map[cid]["index"]) and (a1->index <= map[cid]["Nterm2"]))
//                 {
//                     subd1 = 0;
//                 }
//                 else if((a1->index >= map[cid]["Mterm1"]) and (a1->index <= map[cid]["Mterm2"]))
//                 {
//                     subd1 = 1;
//                 }
//                 else if((a1->index >= map[cid]["Cterm1"]) and (a1->index <= map[cid]["findex"]))
//                 {
//                     subd1 = 2;
//                 }

//                 if((a2->index >= map[cid]["index"]) and (a2->index <= map[cid]["Nterm2"]))
//                 {
//                     subd2 = 0;
//                 }
//                 else if((a2->index >= map[cid]["Mterm1"]) and (a2->index <= map[cid]["Mterm2"]))
//                 {
//                     subd2 = 1;
//                 }
//                 else if((a2->index >= map[cid]["Cterm1"]) and (a2->index <= map[cid]["findex"]))
//                 {
//                     subd2 = 2;
//                 }

//                 contacts.push_back(boost::tuple<int,int,double,int,int>(
//                                        a1->index,
//                                        a2->index,
//                                        dist,
//                                        subd1,
//                                        subd2));
//                 // try
//                 // {
//                 //     contacts.push_back(boost::tuple<int,int,double>(a1->index,a2->index,dist));

//                 // }
//                 // catch (const std::bad_alloc &contacts)
//                 // {
//                 //     std::cout << "Allocation failed for single contact: " << contacts.what() << std::endl;
//                 // }
//             }
//             // i2 += 1;

//             // }
//         }
//         // i1 += 1;
//         // i2 = 0;
//         // i += 1;
//     }

//     // std::cout << "contacts_size: " << contacts.size() << std::endl;
//     // std::cout << "contacts_counted: " << total_contacts << std::endl;
//     return contacts;

//     // if (total_contacts > 0)
//     // {
//     //     return contacts;
//     // }
//     // else
//     // {
//     //     // contacts.push_back(boost::tuple<int,int,double>(-1,-1,0));
//     //     return std::vector<boost::tuple<int,int,double>>(1,(-1,-1,0));
//     // }

// }


// void load_dcd_to_atoms(dcdhandle *dcd,Atom *aa);

// std::vector<boost::tuple<int,int,int,double>> get_contacts_for_chain_later(Atom *alater,
//                                                                            double cutoff,
//                                                                            double tolerance,
//                                                                            std::vector<boost::tuple
//                                                                            <int,int,int,double>> contacts)

SetContacts get_contacts_for_chain_later(vAtoms aa,
                                         double cutoff,
                                         double tolerance,
                                         SetContacts contacts)
{
// #pragma omp parallel
    // std::cout << "Welcome to Contacts - Later!" << std::endl;
    // std::vector<std::vector<boost::tuple<int,int,int,double>>> vec_contacts;
    // std::vector<boost::tuple<int,int,int,double>> cur_contacts;
    SetContacts cur_contacts;
    cur_contacts.reserve(contacts.size());
    double dist;
    double odist;
    double threshold;

    // Parallelize:
    // #pragma omp parallel
    // for(int i=0)
    // for(auto c: contacts)
    // #pragma omp parallel for
    // #pragma omp simd
    // {
    // #pragma omp for
    // #pragma omp parallel num_threads(4)
    // #pragma omp for ordered schedule(dynamic)

    // Failing due to race condition ...? probably not
    // https://software.intel.com/en-us/articles/getting-started-with-openmp

    // works
// #pragma omp parallel for private(dist,odist,threshold,p1,p2,cur_contacts)
        for(int i=0; i<contacts.size(); i++)
        {
            // onoff = -1; // will be set to 0,1 in all cases.
            dist = 0.0;
            odist = 0.0;

            // std::cout << boost::get<0>(c) << std::endl;
            // std::cout << boost::get<1>(c) << std::endl;
            // std::cout << boost::get<2>(c) << std::endl;
            // std::cout << "orig_dist: " << boost::get<3>(c) << std::endl;
            odist = boost::get<2>(contacts[i]);

            // std::cout << std::endl;
            // std::cout << aa[boost::get<0>(c)].x << std::endl;
            // std::cout << aa[boost::get<0>(c)].y << std::endl;
            // std::cout << aa[boost::get<0>(c)].z << std::endl;


            // No speed up. grrrrrrr!
            dist = matom_distance(aa[boost::get<0>(contacts[i])],
                                  aa[boost::get<1>(contacts[i])]);


            // p1.x = aa[boost::get<0>(contacts[i])].x;
            // p1.y = aa[boost::get<0>(contacts[i])].y;
            // p1.z = aa[boost::get<0>(contacts[i])].z;

            // std::cout << std::endl;
            // std::cout << aa[boost::get<1>(c)].x << std::endl;
            // std::cout << aa[boost::get<1>(c)].y << std::endl;
            // std::cout << aa[boost::get<1>(c)].z << std::endl;

            // p2.x = aa[boost::get<1>(contacts[i])].x;
            // p2.y = aa[boost::get<1>(contacts[i])].y;
            // p2.z = aa[boost::get<1>(contacts[i])].z;

            // dist = distance(p1,p2);

            // cutoff: 8
            // tolerance: 2
            if(odist > (cutoff - tolerance))
            {
                // 6 - 8
                threshold = odist + tolerance;
            }
            else
            {
                threshold = cutoff;
            }

            // if((dist < cutoff) or (dist < odist + tolerance))
            if(dist < threshold)
            {
                // index, index, dist, subdomain, subdomain, eh
                cur_contacts.push_back(
                    boost::tuple<int,int,double,int,int,double>
                    (boost::get<0>(contacts[i]),
                     boost::get<1>(contacts[i]),
                     dist,
                     boost::get<3>(contacts[i]),
                     boost::get<4>(contacts[i]),
                     boost::get<5>(contacts[i])));
            }
            // Contact:
            // tuple: <index,index, 0 or 1, current-distance>
            // cur_contacts.push_back(boost::tuple<int,int,int,double>(boost::get<0>(c),
            //                                                         boost::get<1>(c),
            //                                                         onoff,
            //                                                         dist));

            // cur_contacts.clear();
        }
        // } // pragma

        // break;

        return cur_contacts;
        // }
}

SetContacts get_contacts_for_chain_later(vAtoms aa,
                                         double cutoff,
                                         SetContacts contacts)
{
    // Hard cutoff:
    SetContacts cur_contacts;
    cur_contacts.reserve(contacts.size());
    // Vector p1, p2;
    double dist;
    double odist;

    for(int i=0; i<contacts.size(); i++)
    {
        // onoff = -1; // will be set to 0,1 in all cases.
        dist = 0.0;
        odist = 0.0;

        odist = boost::get<2>(contacts[i]);

        // No speed up. grrrrrrr!
        dist = matom_distance(aa[boost::get<0>(contacts[i])],
                              aa[boost::get<1>(contacts[i])]);


        // p1.x = aa[boost::get<0>(contacts[i])].x;
        // p1.y = aa[boost::get<0>(contacts[i])].y;
        // p1.z = aa[boost::get<0>(contacts[i])].z;
        // p2.x = aa[boost::get<1>(contacts[i])].x;
        // p2.y = aa[boost::get<1>(contacts[i])].y;
        // p2.z = aa[boost::get<1>(contacts[i])].z;
        // dist = distance(p1,p2);

        if(dist < cutoff)
        {
            // index, index, dist, subdomain, subdomain, eh
            cur_contacts.push_back(
                boost::tuple<int,int,double,int,int,double>
                (boost::get<0>(contacts[i]),
                 boost::get<1>(contacts[i]),
                 dist,
                 boost::get<3>(contacts[i]),
                 boost::get<4>(contacts[i]),
                 boost::get<5>(contacts[i])));
        }
    }

    return cur_contacts;
}

void output_global_contacts(SetGlobalContacts gc)
{
    std::cout << "Writing Global Contacts to file now." << std::endl;


    // FILE
    FILE * fp_contacts;
    fp_contacts = fopen("emol_mtcontacts.dat", "w+");
    // fprintf(fp_contacts,"\n");

    for(auto f: gc) // frame in global contact array
    {
        for(auto c: f) // chain (but actually dimer) in frame, 156
        {
            for(auto n: c) // 9 situations of 6 neighbors, 0,1,0-1; 0-2,3,4; 1-5,6,7
            {
                fprintf(fp_contacts,"%d ",n.size());
            }
        }
        fprintf(fp_contacts,"\n");
    }

    fclose(fp_contacts);
}

void output_global_contacts_by_subdomain(SetGlobalContacts gc)
{
    std::cout << "Writing Global Contacts by Subdomain to file now." << std::endl;

    int nn, nm, nc, mn, mm, mc, cn, cm, cc;
    nn = nm = nc = mn = mm = mc = cn = cm = cc = 0;

    // new sort by 3s, NMC (coming outward from Alpha or Beta)
    int N, M, C;
    N = M = C;

    // float ext, aext, bext;
    // float asn, aen, awn;
    // float asm, aem, awm;
    // float asc, aec, awc;
    // float bsn, ben, bwn;
    // float bsm, bem, bwm;
    // float bsc, bec, bwc;
    float fN, fM, fC;



    // FILE
    FILE * fp_contacts;
#ifdef TOPO
    fp_contacts = fopen("emol_mtcontacts_by_subdomain_top.dat", "w+");
#else
    fp_contacts = fopen("emol_mtcontacts_by_subdomain.dat", "w+");
#endif

    // || defined (DCD_WRITE_B) || defined (DCD_WRITE) || defined (DCD_WRITE_E)
    // TOPOmt2 = -DDCDREAD -DTOPO -DTOPO_read -DTOPO_mt_BEFORE \
    //     -DMTBUILDMAP -DTOPO_mt_SORT -DMTMAP2_DURING -DMTMAP2_AFTER


    // fprintf(fp_contacts,"\n");


    // FILE
    FILE * fp_contacts3;
#ifdef TOPO
    fp_contacts3 = fopen("emol_mtcontacts_by_subdomain3_top.dat", "w+");
#else
    fp_contacts3 = fopen("emol_mtcontacts_by_subdomain3.dat", "w+");
#endif

    fprintf(fp_contacts3,"#   0-NMC(3)");
    fprintf(fp_contacts3,"                1-NMC(7)");
    fprintf(fp_contacts3,"                0-1-NMC(11)");
    fprintf(fp_contacts3,"              0-2-NMC(15)");
    fprintf(fp_contacts3,"            0-3-NMC(19)");
    fprintf(fp_contacts3,"             0-4-NMC(23)");
    fprintf(fp_contacts3,"             1-5-NMC(27)");
    fprintf(fp_contacts3,"             1-6-NMC(31)");
    fprintf(fp_contacts3,"             1-7-NMC(35)\n");

    // FILE
    FILE * fp_contacts3n;
#ifdef TOPO
    fp_contacts3n = fopen("emol_mtcontacts_by_subdomain3n_top.dat", "w+");
#else
    fp_contacts3n = fopen("emol_mtcontacts_by_subdomain3n.dat", "w+");
#endif

    fprintf(fp_contacts3n,"#     0-NMC(3)");
    fprintf(fp_contacts3n,"              1-NMC(7)");
    fprintf(fp_contacts3n,"              0-1-NMC(11)");
    fprintf(fp_contacts3n,"           0-2-NMC(15)");
    fprintf(fp_contacts3n,"           0-3-NMC(19)");
    fprintf(fp_contacts3n,"           0-4-NMC(23)");
    fprintf(fp_contacts3n,"           1-5-NMC(27)");
    fprintf(fp_contacts3n,"           1-6-NMC(31)");
    fprintf(fp_contacts3n,"           1-7-NMC(35)\n");



    // fprintf(fp_contacts,"\n");


    int count_f;
    count_f = 0;

    for(auto f: gc) // frame in global contact array
    {
        count_f += 1;

        for(auto c: f) // chain (but actually dimer) in frame, 156
        {
            for(auto n: c) // 9 situations of 6 neighbors, 0,1,0-1; 0-2,3,4; 1-5,6,7
            {
                // fprintf(fp_contacts,"%d ",n.size());

                nn = nm = nc = mn = mm = mc = cn = cm = cc = 0;
                N = M = C = 0;

                for(auto contact: n)
                {
                    // std::cout << contact.get<3>() << " " << contact.get<4>() << std::endl;
                    if((contact.get<3>() == 0) and (contact.get<4>() == 0))
                    {
                        nn += 1;
                        N += 1;
                    }
                    else if((contact.get<3>() == 0) and (contact.get<4>() == 1))
                    {
                        nm += 1;
                        N += 1;
                    }
                    else if((contact.get<3>() == 0) and (contact.get<4>() == 2))
                    {
                        nc += 1;
                        N += 1;
                    }
                    else if((contact.get<3>() == 1) and (contact.get<4>() == 0))
                    {
                        mn += 1;
                        M += 1;
                    }
                    else if((contact.get<3>() == 1) and (contact.get<4>() == 1))
                    {
                        mm += 1;
                        M += 1;
                    }
                    else if((contact.get<3>() == 1) and (contact.get<4>() == 2))
                    {
                        mc += 1;
                        M += 1;
                    }
                    else if((contact.get<3>() == 2) and (contact.get<4>() == 0))
                    {
                        cn += 1;
                        C += 1;
                    }
                    else if((contact.get<3>() == 2) and (contact.get<4>() == 1))
                    {
                        cm += 1;
                        C += 1;
                    }
                    else if((contact.get<3>() == 2) and (contact.get<4>() == 2))
                    {
                        cc += 1;
                        C += 1;
                    }

                }
                fprintf(fp_contacts,"%5d %5d %5d %5d %5d %5d %5d %5d %5d %5d\n",
                        n.size(),
                        nn,nm,nc,
                        mn,mm,mc,
                        cn,cm,cc);


                // 21 + 1 (n.size()) = 22
                // ext = aext = bext = 0;
                // asn = aen = awn = 0;
                // asm = aem = awm = 0;
                // asc = aec = awc = 0;
                // bsn = ben = bwn = 0;
                // bsm = bem = bwm = 0;
                // bsc = bec = bwc = 0;
                // ext =


                // 3n-file
                fN = fM = fC = 0.0;

                if (n.size() != 0)
                {
                    fN = float(N) / n.size();
                    fM = float(M) / n.size();
                    fC = float(C) / n.size();
                }

                fprintf(fp_contacts3n,"%5d %4.2f %4.2f %4.2f  ",
                        n.size(),
                        fN,
                        fM,
                        fC);


                // 3-file
                fprintf(fp_contacts3,"%5d %5d %5d %5d ",
                        n.size(),
                        N,
                        M,
                        C);

            }
            // fprintf(fp_contacts,"\n");

            // 3n-file
            fprintf(fp_contacts3n,"\n");
            // 3-file
            fprintf(fp_contacts3,"\n");

            // fprintf(fp_contacts3,"%5d %4.2f %4.2f",
            //         n.size(),
            //         )

        }

        fprintf(fp_contacts,"# %d\n",count_f);
        fprintf(fp_contacts3,"# %d\n",count_f);
        fprintf(fp_contacts3n,"# %d\n",count_f);
    }
    fclose(fp_contacts);
    fclose(fp_contacts3);
    fclose(fp_contacts3n);
}

MtNeighbors get_map_of_mtneighbors(vIndexGroup isel_chain,vAtoms aa,DimerList dimers)
{
    // printf("Welcome to get_map_of_mtneighbors!\n");
    // std::cout << matrix.size() << std::endl;
    // std::vector<std::vector<int>> matrix(dimers.size(), std::vector<int>(8,-1));
    MtNeighbors matrix(dimers.size(), std::vector<int>(8,-1));

    // int ic;
    // ic = 0;

    Vector centroid;
    // centroid.print_Vector();
    std::vector<Vector> centroids;


    for(auto ch: isel_chain)
    {
        // std::cout << ch.size() << std::endl;
        centroid = get_centroid(ch,aa); // ch: list of indices, aa: allatoms_ref
        // centroid.print_Vector();
        centroids.push_back(centroid);

        // std::cout << "Return_Centroid: "
        //           << centroid[0]
        //           << " "
        //           << centroid[1]
        //           << " "
        //           << centroid[2]
        //           << std::endl;

        // ic += 1;
    }
    // std::cout << "# of centroids: " << centroids.size() << std::endl;
    // exit(0);

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

                // 90: 44 * 1.414 = ~62 + ...
                if (vdist_mag < 75.0)
                {
                    candidates.push_back(m2);
                }
            }
        }
        chain_candidates.push_back(candidates);
        candidates.clear();
    }
    // std::cout << "Candidates acquired." << std::endl;
    // exit(0);


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

    // MTMAP: the diagram
    //     West
    //     4  7
    //  2  0  1  5
    //     3  6
    //     East

    for(auto cc: chain_candidates)
    {
        // std::cout << "monomer: " << ican << std::endl;

        // ican == corresponding ALPHA from either an
        // actual ALPHA or BETA candidate
        // i.e. 72, 73 ==> 36 chain.
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
            // Axis Vector between Alpha and Beta monomers.
            // normalized.
            // magnitude.
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

            // std::cout << dimers[pdim].first
            //           << "-"
            //           << dimers[pdim].second
            //           << "   can: "
            //           << can
            //           << std::endl;


            // Vector between the candidate (can) and ..
            // the alpha monomer (ican)
            vdist = get_vector(centroids[ican],centroids[can]);
            // vdist = get_vector(centroids[can],centroids[ican]);
            vdist_n = normalize(vdist);
            vdist_mag = magnitude(vdist);
            // std::cout << can << " " << vdist_mag << " " << std::endl;



            dsin = get_sintheta(avec_n,vdist_n);
            dcos = get_costheta(avec_n,vdist_n);

            // std::cout << "sin: " << dsin << std::endl;
            // std::cout << "cos: " << dcos << std::endl;
            // std::cout << std::endl;

            // Find .. within protofilament. (longitudinal)
            // 54-55   can: 29  (8.8 deg)
            // sin: 0.153011
            // cos: 0.988224
            // if((dsin < 0.12) && (std::abs(dcos) > 0.88))
            if(std::abs(dcos) > 0.92)
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
            // End of protofilament (longitudinal).

            // find West/East Neighbors.
            // cos(69.94) = 0.343, if less than 0.343 => angle >= 70 deg. want greater
            // sin(65.5) = 0.91, if greater than 0.91 => angle >= 65.5 deg. want less
            else if ((dsin > 0.91) && (std::abs(dcos) < 0.343))
            {
                if(dcos < 0)
                {
                    // WEST:
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
                    // EAST:
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

void print_mt_map(MtNeighbors mt_matrix)
{
    // KEEP THIS.
    // Print Map of MT neigbors.
    for(auto mt_ch: mt_matrix)
    {
        for(auto mt_n: mt_ch)
        {
            std::cout << std::setw(4) << mt_n << " ";
        }
        std::cout << std::endl;
    }


    // KEEP (still the old way.)
    // Print Map of MT neighbors. (old way..?)
    // std::vector<std::vector<int>>::iterator itmap;
    // std::vector<int>::iterator itmap_n;
    // for(itmap = mt_matrix.begin(); itmap != mt_matrix.end(); itmap++)
    // {
    //     // std::cout << (*itmap).size() << std::endl;

    //     for(itmap_n = (*itmap).begin(); itmap_n != (*itmap).end(); itmap_n++)
    //     {
    //         std::cout << std::setw(4) << (*itmap_n) << " ";
    //     }
    //             std::cout << std::endl;
    // }

}

void print_global_contacts(SetGlobalContacts gc)
{
    // std::cout << "Writing Global Contacts to file now." << std::endl;

    int fn = 0;
    int cn = 0;
    int nn = 0;
    int gcn = 0; // global contact number.

    // FILE
    for(auto f: gc) // frame in global contact array
    {
        fn += 1;
        std::cout << "Frame: " << fn << std::endl;

        for(auto c: f) // chain (but actually dimer) in frame, 156
        {
            std::cout << "Chain: " << cn << std::endl;

            for(auto n: c) // 9 situations of 6 neighbors, 0,1,0-1; 0-2,3,4; 1-5,6,7
            {
                std::cout << "Neigbor: " << nn << " " << std::endl;

                for(auto c1: n)
                {
                    std::cout << c1.get<0>() << "-" << c1.get<1>() << " ";

                    gcn += 1;
                }
                std::cout << std::endl;
                nn += 1;
            }
            std::cout << std::endl;

            nn = 0;
            cn += 1;
        }
        // fprintf(fp_contacts,"\n");
    }
    // fclose(fp_contacts);
    std::cout << "The total number of contacts is: " << gcn << std::endl;
}

void print_global_contacts_count(SetGlobalContacts gc)
{
    // std::cout << "Writing Global Contacts to file now." << std::endl;

    int fn = 0;
    int cn = 0;
    int nn = 0;
    int gcn = 0; // global contact number.

    // FILE
    for(auto f: gc) // frame in global contact array
    {
        fn += 1;
        std::cout << "Frame: " << fn << std::endl;

        for(auto c: f) // chain (but actually dimer) in frame, 156
        {
            std::cout << "Chain: " << cn << std::endl;

            for(auto n: c) // 9 situations of 6 neighbors, 0,1,0-1; 0-2,3,4; 1-5,6,7
            {
                std::cout << "  n: " << nn << " has " << n.size() << " contacts."<< std::endl;
                for(auto c1: n)
                {
                    // std::cout << c1.get<0>() << "-" << c1.get<1>() << " ";
                    gcn += 1;
                }
                // std::cout << std::endl;
                nn += 1;
            }
            // std::cout << std::endl;

            nn = 0;
            cn += 1;
        }
        // fprintf(fp_contacts,"\n");
    }
    // fclose(fp_contacts);
    std::cout << "The total number of contacts is: " << gcn << std::endl;
}

void write_contacts_to_file(FILE *fp_topology,SetContacts contact_set)
{
    // TOPOLOGY

    // fprintf(fp_topology,"# contact set\n");
    for(auto c: contact_set)
    {
        fprintf(fp_topology,"%5d %5d    1  %9.5f   %8.5f\n",
                c.get<0>(),
                c.get<1>(),
                c.get<2>(),
                c.get<5>());
    }
}

SetContacts read_contacts_from_file(char filename[40])
{
    SetContacts contacts;
    std::string line;
    std::string line_after1;
    std::string line_after2;
    std::string fn(filename);

    std::string str_index1;
    std::string str_index2;
    std::string str_dist;
    std::string str_eh;

    int now_read;
    now_read = -1;

    int index1, index2;
    double dist, eh;
    index1 = index2 = -1;
    dist = eh = 0.0;

    std::cout << "Reading contacts from topology file: "
              << filename
              << std::endl;

    std::ifstream cfile(fn.c_str());

    std::size_t found1;
    std::size_t found2;

    // = str.find(str2);

    if(cfile.is_open())
    {
        std::cout << "file is now open." << std::endl;

        while(getline(cfile,line))
        {
            if (boost::algorithm::starts_with(line,"#") == 1)
            {
                continue;
            }
            else if (boost::algorithm::starts_with(line,";") == 1)
            {
                continue;
            }
            else if (boost::algorithm::starts_with(line,"[ native ]") == 1)
            {
                now_read += 1;
                continue;
            }


            if(now_read != -1)
            {
                now_read += 1;
                // std::cout << now_read << ":";
                // std::cout << line << std::endl;

                found1 = line.find(" ");

                // Gets the line that starts/ends with a space.
                if (found1 == std::string::npos)
                {
                    std::cout << "lines corresponding to contacts: " << now_read << std::endl;
                    break;
                }

                // Index1:
                found1 = line.find_first_not_of(" ");
                found2 = line.find_first_not_of("0123456789",found1+1);
                str_index1 = line.substr(found1,found2);

                // Index2:
                line_after1 = line.substr(found2+1,std::string::npos);
                found1 = line_after1.find_first_not_of(" ");
                found2 = line_after1.find_first_not_of("0123456789",found1+1);
                str_index2 = line_after1.substr(found1,found2);

                // dist: (with a skip for that weird "1" in the topology)
                line_after2 = line_after1.substr(found2+1,std::string::npos);
                found1 = line_after2.find_first_not_of(" ");
                line_after2 = line_after2.substr(found1+1,std::string::npos);
                found1 = line_after2.find_first_not_of(" ");
                found2 = line_after2.find_first_not_of("0123456789.",found1+1);
                str_dist = line_after2.substr(found1,found2);

                // eh:
                str_eh = line_after2.substr(found2+1,std::string::npos);

                // std::cout << found1 << " " << found2 << std::endl;
                // std::cout << line_after1 << std::endl;

                // if(now_read > 500)
                // {
                //     break;
                // }


                // // ---------------FORMER ----------
                // str_index1 = line.substr(1,6);
                // str_index2 = line.substr(9,6);
                // str_dist = line.substr(22,9);
                // str_eh = line.substr(34,8);

                // // std::cout << str_index1 << " "
                // //           << str_index2 << " "
                // //           << str_dist << " "
                // //           << str_eh << "\n"
                // //           << std::endl;
                // // ---------------FORMER ----------



                index1 = atoi(str_index1.c_str());
                index2 = atoi(str_index2.c_str());
                dist = std::stod(str_dist.c_str());
                eh = std::stod(str_eh.c_str());


                // std::cout << index1 << " "
                //           << index2 << " "
                //           << dist << " "
                //           << eh << "\n"
                //           << std::endl;

                contacts.push_back(boost::tuple<int,int,double,int,int,double>(
                                       index1,
                                       index2,
                                       dist,
                                       -1,
                                       -1,
                                       eh));
            }
        }
        cfile.close();
    }
    std::cout << "file is now closed." << std::endl;
    std::cout << "# of contacts from topology file:  " <<
        contacts.size() << std::endl;

    return contacts;
}

void print_set_contacts(SetContacts cn)
{
    for(auto c: cn)
    {
        std::cout << c.get<0>() << "-"
                  << c.get<1>() << "  ["
                  << c.get<2>() << "]   "
                  << c.get<3>() << "--"
                  << c.get<4>() << " "
                  << c.get<5>() << "\n";
    }
}

SetContacts set_eh_contacts(SetContacts cn,double eh)
{
    //       v-------Major fix. Now it captures tuple by reference.
    for(auto &c: cn)
    {
        // std::cout << c.get<5>() << " " << eh << std::endl;
        c.get<5>() = eh;
        // std::cout << c.get<5>() << std::endl;

    }

    // std::cout << "assignment complete.\n";
    // for(auto c: cn)
    // {
    //     std::cout << c.get<5>() << std::endl;
    // }

    return cn;
}

SetContacts set_mtsubdomain_flag_contacts(SetContacts cn,MtIndexMap mtmap)
{
    int index;
    index = -1;

    //       v-------Major fix. Now it captures tuple by reference.
    for(auto &c: cn)
    {
        index = c.get<0>();
        // mtentry["chaintype"] = betabool;
        // mtentry["index"] = low_index;
        // mtentry["findex"] = high_index;
        // mtentry["Nterm2"] = Nterm2;
        // mtentry["Mterm1"] = Mterm1;
        // mtentry["Mterm2"] = Mterm2;
        // mtentry["Cterm1"] = Cterm1;

        // std::cout << c.get<0>() << " " << c.get<1>() << std::endl;
        // std::cout << chain["index"] << " " << chain["findex"] << std::endl;

        for(auto chain: mtmap)
        {
            if((index >= chain["index"]) and
               (index <= chain["findex"]))
            {
                if(index <= chain["Nterm2"])
                {
                    c.get<3>() = 0;
                }
                else if(index <= chain["Mterm2"])
                {
                    c.get<3>() = 1;
                }
                else
                {
                    c.get<3>() = 2;
                }
                break;
            }
        }

        index = c.get<1>();

        for(auto chain: mtmap)
        {
            if((index >= chain["index"]) and
               (index <= chain["findex"]))
            {
                if(index <= chain["Nterm2"])
                {
                    c.get<4>() = 0;
                }
                else if(index <= chain["Mterm2"])
                {
                    c.get<4>() = 1;
                }
                else
                {
                    c.get<4>() = 2;
                }
                break;
            }
        }
    }

    // std::cout << "assignment complete.\n";
    // for(auto c: cn)
    // {
    //     std::cout << c.get<5>() << std::endl;
    // }

    return cn;
}

void output_framecontact_set(SetFrameContacts framecontact_set)
{
    std::cout << "Writing FrameContacts to file now." << std::endl;

    // FILE
    FILE * fp_contacts;
    fp_contacts = fopen("emol_mtcontacts_list.dat", "w+");
    fprintf(fp_contacts,"# emol_mtcontacts_list.dat\n");

    // for(auto f: framecontact_set) // frame in global contact array
    // {
        // for(auto c: f) // chain (but actually dimer) in frame, 156
        // {
            // for(auto cn: c) // 9 situations of 6 neighbors, 0,1,0-1; 0-2,3,4; 1-5,6,7
            // {
            //     fprintf(fp_contacts,"%d ",cn.size());
            // }
        // }
    //     fprintf(fp_contacts,"\n");
    // }
    fclose(fp_contacts);
}

void print_chain_index_boundaries(MtIndexMap mtmap)
{
    for(auto m: mtmap)
    {
        // std::cout << m.size();
        std::cout << " " << std::setw(6) << m["index"];
        std::cout << " " << std::setw(6) << m["findex"];
        std::cout << " " << std::setw(6) << m["findex"] - m["index"] + 1;
        std::cout << std::endl;
    }
}

// SetChains sort_contacts2(SetContacts cn,MtIndexMap mtmap)
SetChains sort_contacts2(SetContacts cn,MtIndexMap mtmap,MtNeighbors matrix)
{
    std::cout << "mtmap: " << mtmap.size() << std::endl;
    std::cout << "matrix: "  << matrix.size() << std::endl;

    // SetContacts cn0,cn1,cn2,cn3,cn4,cn5,cn6,cn7,cn8; // 9
    // Contact contact;
    SetContacts cn0;
    SetNeighbors neighbor_set;
    SetChains chain_set;
    int cin1,fin1,cin2,fin2,numchain,m1,m2,c1,c2;
    int cc;
    cc = -1;
    // cin, fin: chain index, final index
    //
    // m: monomer
    cin1 = cin2 = fin1 = fin2 = c1 = c2 = numchain = m1 = m2 = -1;

    std::vector<std::pair<int,int>> combos;
    combos.push_back(std::make_pair(0,0));
    combos.push_back(std::make_pair(1,1));
    combos.push_back(std::make_pair(0,1));
    combos.push_back(std::make_pair(0,2));
    combos.push_back(std::make_pair(0,3));
    combos.push_back(std::make_pair(0,4));
    combos.push_back(std::make_pair(1,5));
    combos.push_back(std::make_pair(1,6));
    combos.push_back(std::make_pair(1,7));
    // global (not present)  | SetGlobalcontacts
    // chains (1/2)          | SetChains
    // 9 neighbors           | SetNeighbors
    // contacts              | SetContacts

    // for(auto m: matrix)
    // {
    //     // 0-310, 156 total
    //     std::cout << " " << std::setw(3) << m[0];
    //     std::cout << " " << std::setw(3) << m[1];
    //     std::cout << " " << std::setw(3) << m[2];
    //     std::cout << " " << std::setw(3) << m[3];
    //     std::cout << " " << std::setw(3) << m[4];
    //     std::cout << " " << std::setw(3) << m[5];
    //     std::cout << " " << std::setw(3) << m[6];
    //     std::cout << " " << std::setw(3) << m[7];
    //     std::cout << std::endl;
    // }

    // for(auto m: mtmap)
    // {
    //     std::cout << m.size();
    //     std::cout << " " << m["index"];
    //     std::cout << " " << m["findex"];
    //     std::cout << std::endl;
    // }

    // n is all 8 neighbors --> 9 situations.
    // for(auto n: matrix)
    // 9 situations:
    // 0, 1, 0-1
    // 0: 2, 3, 4 (SEW)
    // 1: 5, 6, 7 (NEW)

    std::vector<int> n(8);

    for(int i=0; i<matrix.size(); i++)
    {
        // i: 0 - 155
        n = matrix[i];

        for(auto cb: combos)
        {
            cn0.clear();

            // m: monomer 15 in contact 40.
            // 14  15  -1  16  12  40  17  13
            // cb: 0-0, 1-1, 0-1, 0-2, 0-3, 0-4, 1-5,1-6,1-7
            m1 = n[cb.first];
            m2 = n[cb.second];

            if((m1 == -1) or (m2 == -1))
            {
                neighbor_set.push_back(cn0);
                continue;
            }
            else
            {
                cin1 = mtmap[m1]["index"];
                fin1 = mtmap[m1]["findex"];
                cin2 = mtmap[m2]["index"];
                fin2 = mtmap[m2]["findex"];

                // std::cout << cb.first << "-" << cb.second;
                // std::cout << "   " << n[cb.first] << "  " << n[cb.second];
                // // << std::endl;
                // std::cout << " | " << cin1 << " " << fin1
                //           << " " << cin2 << " " << fin2
                //           << " " << cn.size() << std::endl;

                // Begin:
                cc = -1;
                for(auto c: cn)
                {
                    cc += 1;
                    c1 = c.get<0>();
                    c2 = c.get<1>();

                    // PRINT HERE?
                    // std::cout << "c12: " << c1 << " " << c2 << std::endl;
                    // std::cout << "i12: " << cin1 << " " << fin1 << std::endl;
                    // std::cout << "n12: " << cin2 << " " << fin2 << std::endl;

                    if(((c1 >= cin1) and (c1 <= fin1)
                        and (c2 >= cin2) and (c2 <= fin2)))
                    {
                        // std::cout << c.get<0>() << " ";
                        // std::cout << c.get<1>() << " ";
                        // std::cout << c.get<2>() << " ";
                        // std::cout << c.get<3>() << " ";
                        // std::cout << c.get<4>() << " ";
                        // std::cout << c.get<5>() << " ";
                        // std::cout << std::endl;

                        // contact = c;
                        // cn0.push_back(contact);
                        cn0.push_back(c);
                        // cn.erase(cc);
                    }
                }
            }
            // std::cout << "cn-size: " << cn0.size() << std::endl;
            neighbor_set.push_back(cn0);

        } // next combo: 0, 1, 0-1, 0-234, 1-567;
        chain_set.push_back(neighbor_set);
        neighbor_set.clear();
    }

    // shows that all should be 9.
    // for(auto ch: chain_set)
    // {
    //     std::cout << " " << ch.size();
    // }
    // std::cout << std::endl;


    // for(auto ch: chain_set)
    // {
    //     for(auto n: ch)
    //     {
    //         for(auto ct: n)
    //         {
    //             std::cout << ct.get<0>() << "-" << ct.get<1>() << std::endl;
    //         }
    //     }
    // }

    std::cout << "Chain_size: " << chain_set.size() << std::endl;
    return chain_set;
}

// SetContacts sort_contacts(SetContacts cn,MtIndexMap mtmap,int c1, int c2)
// {
//     SetContacts contacts;

//     // std::cout << cn.size() << std::endl;
//     // std::cout << mtmap.size() << std::endl;
//     // std::cout << c1 << " " << c2 << std::endl;
//     // std::cout << mtmap[c1]["index"] << " " << mtmap[c1]["findex"] << std::endl;
//     // std::cout << mtmap[c2]["index"] << " " << mtmap[c2]["findex"] << std::endl;

//     for(auto c: cn)
//     {
//         // std::cout << c.get<0>() << " " << c.get<1>() << std::endl;

//         if((c.get<0>() >= mtmap[c1]["index"])
//            and (c.get<0>() <= mtmap[c1]["findex"])
//            and (c.get<1>() >= mtmap[c2]["index"])
//            and (c.get<1>() <= mtmap[c2]["findex"]))
//         {
//             contacts.push_back(c);
//         }
//     }
//     // exit(0);

//     return contacts;
// }
