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


/* ---------------------------------------------------------
   headers
   --------------------------------------------------------- */
#include "debug.h"
#include "contacts.hpp"
#include "md.h"
#include "system.hpp"
#include "microtubule.hpp"
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

// void get_contacts(Atom *a1,char *argv)
// {
//     debug("welcome to get_contacts. 1 (intra)\n");
// }

// void get_contacts_for_chain()
// std::vector<boost::tuple<int,int,int,double>> get_contacts_for_chain(std::vector <Atom> chain1,
//                                                                      std::vector <Atom> chain2,
//                                                                      float cutoff)
// std::vector<boost::tuple<int,int,double>> get_contacts_for_chain(std::vector <Atom> chain1,
//                                                                  std::vector <Atom> chain2,
//                                                                  float cutoff)
// std::vector<boost::tuple<int,int,double>> get_contacts_for_chain(std::vector <Atom> chain1,
//                                                                  std::vector <Atom> chain2,
//                                                                  float cutoff,
//                                                                  std::vector<boost::tuple<int,int,double>> contacts)
// void get_contacts_for_chain(std::vector <Atom> chain1,
//                             std::vector <Atom> chain2,
//                             float cutoff,
//                             std::vector<boost::tuple<int,int,double>> contacts)
SetContacts get_contacts_for_chain(std::vector <Atom> chain1,
                                   std::vector <Atom> chain2,
                                   double cutoff)
{
    // INTERCHAIN, Not MT

    // std::cout << "Welcome to get_contacts_for_chain!" << std::endl;
    // std::cout << "Cutoff: " << cutoff << std::endl;
    // std::vector<boost::tuple<int,int,double>> contacts;
    SetContacts contacts;
    // std::cout << "chain1: " << chain1.size() << std::endl;
    // std::cout << "chain2: " << chain2.size() << std::endl;


    // int i1, i2;
    // i1 = i2 = 0;
    // int i = 0;

    Vector a1v,a2v;
    double dist = 0.0;

    int subd1 = -1; // subdomains 1 & 2
    int subd2 = -1;

    // int total_contacts = 0;

    for(auto a1: chain1)
    {
        // std::cout << i1 << " " << a1.index << std::endl;
        // std::cout << a1.index << std::endl;
        a1v.x = a1.x;
        a1v.y = a1.y;
        a1v.z = a1.z;

        for(auto a2: chain2)
        {
            // if(a2.index <)

            // if(i2 <= i1)
            // {
            //     continue;
            // }


            // Directive: Exclude an index being in contact with itself,
            // or its +/- 2 neighbor.
            if((a2.index - 2 >= a1.index) and (a2.index + 2 <= a1.index))
               // or
               // (a2.index - 1 == a1.index))
            {
                continue;
            }


                // std::cout << "\t" << i2 << "\t" << a2.index << std::endl;
                // std::cout << "\t" << a2.index << std::endl;

                a2v.x = a2.x;
                a2v.y = a2.y;
                a2v.z = a2.z;

                dist = distance(a1v,a2v);


                // a1.print_coords();
                // a2.print_coords();


                if(dist <= cutoff)
                {
                    // contacts.push_back(boost::tuple<int,int,int,double>(a1.index,a2.index,1,dist));
                    // total_contacts += 1;

                    contacts.push_back(boost::tuple<int,int,double,int,int>(
                                           a1.index,
                                           a2.index,
                                           dist,
                                           subd1,
                                           subd2));
                    // try
                    // {
                    //     contacts.push_back(boost::tuple<int,int,double>(a1.index,a2.index,dist));

                    // }
                    // catch (const std::bad_alloc &contacts)
                    // {
                    //     std::cout << "Allocation failed for single contact: " << contacts.what() << std::endl;
                    // }
                }
                // i2 += 1;

            // }
        }
        // i1 += 1;
        // i2 = 0;
        // i += 1;
    }

    // std::cout << "contacts_size: " << contacts.size() << std::endl;
    // std::cout << "contacts_counted: " << total_contacts << std::endl;
    return contacts;

    // if (total_contacts > 0)
    // {
    //     return contacts;
    // }
    // else
    // {
    //     // contacts.push_back(boost::tuple<int,int,double>(-1,-1,0));
    //     return std::vector<boost::tuple<int,int,double>>(1,(-1,-1,0));
    // }



}

SetContacts get_contacts_for_chain(std::vector <Atom> chain1,
                                   std::vector <Atom> chain2,
                                   double cutoff,
                                   MtIndexMap map,
                                   int cid1,
                                   int cid2)
{
    // INTERCHAIN, w/MT

    // std::cout << "Welcome to get_contacts_for_chain! MT, inter-chain." << std::endl;
    // std::cout << "Cutoff: " << cutoff << std::endl;
    // std::vector<boost::tuple<int,int,double>> contacts;
    SetContacts contacts;
    // std::cout << "chain1: " << chain1.size() << std::endl;
    // std::cout << "chain2: " << chain2.size() << std::endl;


    // int i1, i2;
    // i1 = i2 = 0;
    // int i = 0;

    Vector a1v,a2v;
    double dist = 0.0;

    int subd1 = -1; // subdomains 1 & 2
    int subd2 = -1;

    // int total_contacts = 0;

    for(auto a1: chain1)
    {
        // std::cout << i1 << " " << a1.index << std::endl;
        // std::cout << a1.index << std::endl;
        a1v.x = a1.x;
        a1v.y = a1.y;
        a1v.z = a1.z;

        for(auto a2: chain2)
        {
            // if(a2.index <)

            // if(i2 <= i1)
            // {
            //     continue;
            // }


            // Directive: Exclude an index being in contact with itself,
            // or its +/- 2 neighbor.
            if((a2.index - 2 >= a1.index) and (a2.index + 2 <= a1.index))
                // or
                // (a2.index - 1 == a1.index))
            {
                continue;
            }


            // std::cout << "\t" << i2 << "\t" << a2.index << std::endl;
            // std::cout << "\t" << a2.index << std::endl;

            a2v.x = a2.x;
            a2v.y = a2.y;
            a2v.z = a2.z;

            dist = distance(a1v,a2v);


            // a1.print_coords();
            // a2.print_coords();


            if(dist <= cutoff)
            {
                // contacts.push_back(boost::tuple<int,int,int,double>(a1.index,a2.index,1,dist));
                // total_contacts += 1;

                if((a1.index >= map[cid1]["index"]) and (a1.index <= map[cid1]["Nterm2"]))
                {
                    subd1 = 0;
                }
                else if((a1.index >= map[cid1]["Mterm1"]) and (a1.index <= map[cid1]["Mterm2"]))
                {
                    subd1 = 1;
                }
                else if((a1.index >= map[cid1]["Cterm1"]) and (a1.index <= map[cid1]["findex"]))
                {
                    subd1 = 2;
                }

                if((a2.index >= map[cid2]["index"]) and (a2.index <= map[cid2]["Nterm2"]))
                {
                    subd2 = 0;
                }
                else if((a2.index >= map[cid2]["Mterm1"]) and (a2.index <= map[cid2]["Mterm2"]))
                {
                    subd2 = 1;
                }
                else if((a2.index >= map[cid2]["Cterm1"]) and (a2.index <= map[cid2]["findex"]))
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

                contacts.push_back(boost::tuple<int,int,double,int,int>
                                   (a1.index,
                                    a2.index,
                                    dist,
                                    subd1,
                                    subd2));
                // try
                // {
                //     contacts.push_back(boost::tuple<int,int,double>(a1.index,a2.index,dist));

                // }
                // catch (const std::bad_alloc &contacts)
                // {
                //     std::cout << "Allocation failed for single contact: " << contacts.what() << std::endl;
                // }
            }
            // i2 += 1;

            // }
        }
        // i1 += 1;
        // i2 = 0;
        // i += 1;
    }

    // std::cout << "contacts_size: " << contacts.size() << std::endl;
    // std::cout << "contacts_counted: " << total_contacts << std::endl;
    return contacts;

    // if (total_contacts > 0)
    // {
    //     return contacts;
    // }
    // else
    // {
    //     // contacts.push_back(boost::tuple<int,int,double>(-1,-1,0));
    //     return std::vector<boost::tuple<int,int,double>>(1,(-1,-1,0));
    // }



}


SetContacts get_contacts_for_chain(std::vector <Atom> chain1,
                                   double cutoff)
{
    // SINGLE, INTRACHAIN, Not MT

    // std::cout << "Welcome to get_contacts_for_chain! MT, intra." << std::endl;
    // std::cout << "Cutoff: " << cutoff << std::endl;
    // std::vector<boost::tuple<int,int,double>> contacts;
    SetContacts contacts;
    // std::cout << "chain1: " << chain1.size() << std::endl;
    // std::cout << "chain2: " << chain2.size() << std::endl;


    // int i1, i2;
    // i1 = i2 = 0;
    // int i = 0;

    Vector a1v,a2v;
    double dist = 0.0;

    int subd1 = -1; // subdomains 1 & 2
    int subd2 = -1;

    // int total_contacts = 0;

    for(auto a1: chain1)
    {
        // std::cout << i1 << " " << a1.index << std::endl;
        // std::cout << a1.index << std::endl;
        a1v.x = a1.x;
        a1v.y = a1.y;
        a1v.z = a1.z;

        for(auto a2: chain1)
        {
            // Directive: Exclude an index being in contact with itself,
            // or its +/- 2 neighbor.
            // if((a2.index - 2 >= a1.index) and (a2.index + 2 <= a1.index))
            // {
            //     continue;
            // }
            if (a2.index - 2 <= a1.index)
            {
                continue;
            }


                // std::cout << "\t" << i2 << "\t" << a2.index << std::endl;
                // std::cout << "\t" << a2.index << std::endl;

                a2v.x = a2.x;
                a2v.y = a2.y;
                a2v.z = a2.z;

                dist = distance(a1v,a2v);


                // a1.print_coords();
                // a2.print_coords();


                if(dist <= cutoff)
                {
                    // contacts.push_back(boost::tuple<int,int,int,double>(a1.index,a2.index,1,dist));
                    // total_contacts += 1;

                    contacts.push_back(boost::tuple<int,int,double,int,int>(
                                           a1.index,
                                           a2.index,
                                           dist,
                                           subd1,
                                           subd2));
                    // try
                    // {
                    //     contacts.push_back(boost::tuple<int,int,double>(a1.index,a2.index,dist));

                    // }
                    // catch (const std::bad_alloc &contacts)
                    // {
                    //     std::cout << "Allocation failed for single contact: " << contacts.what() << std::endl;
                    // }
                }
                // i2 += 1;

            // }
        }
        // i1 += 1;
        // i2 = 0;
        // i += 1;
    }

    // std::cout << "contacts_size: " << contacts.size() << std::endl;
    // std::cout << "contacts_counted: " << total_contacts << std::endl;
    return contacts;

    // if (total_contacts > 0)
    // {
    //     return contacts;
    // }
    // else
    // {
    //     // contacts.push_back(boost::tuple<int,int,double>(-1,-1,0));
    //     return std::vector<boost::tuple<int,int,double>>(1,(-1,-1,0));
    // }



}

SetContacts get_contacts_for_chain(std::vector <Atom> chain1,
                                   double cutoff,
                                   MtIndexMap map,
                                   int cid)
{
    // SINGLE, INTRACHAIN, w/MT

    // std::cout << "Welcome to get_contacts_for_chain!" << std::endl;
    // std::cout << "Cutoff: " << cutoff << std::endl;
    // std::vector<boost::tuple<int,int,double>> contacts;
    SetContacts contacts;
    // std::cout << "chain1: " << chain1.size() << std::endl;
    // std::cout << "chain2: " << chain2.size() << std::endl;


    // int i1, i2;
    // i1 = i2 = 0;
    // int i = 0;

    Vector a1v,a2v;
    double dist = 0.0;

    int subd1 = -1; // subdomains 1 & 2
    int subd2 = -1;

    // int total_contacts = 0;

    // int ch1 = -1;
    // int ch2 = -1;

    for(auto a1: chain1)
    {
        // ch1 += 1;
        // ch2 = -1;

        // std::cout << i1 << " " << a1.index << std::endl;
        // std::cout << a1.index << std::endl;
        a1v.x = a1.x;
        a1v.y = a1.y;
        a1v.z = a1.z;

        for(auto a2: chain1)
        {
            // ch2 += 1;

            // Directive: Exclude an index being in contact with itself,
            // or its +/- 2 neighbor.
            // if((a2.index - 2 >= a1.index) and (a2.index + 2 <= a1.index))
            // {
            //     continue;
            // }
            if (a2.index - 2 <= a1.index)
            {
                continue;
            }


            // std::cout << "\t" << i2 << "\t" << a2.index << std::endl;
            // std::cout << "\t" << a2.index << std::endl;



            a2v.x = a2.x;
            a2v.y = a2.y;
            a2v.z = a2.z;

            dist = distance(a1v,a2v);


            // a1.print_coords();
            // a2.print_coords();


            if(dist <= cutoff)
            {
                // contacts.push_back(boost::tuple<int,int,int,double>(a1.index,a2.index,1,dist));
                // total_contacts += 1;

                // std::cout << "chaintype: " << m["chaintype"] << " \n"
                //           << "index: " << m["index"] << " \n"
                //           << "Nterm2: " << m["Nterm2"] << " \n"
                //           << "Mterm1: " << m["Mterm1"] << " \n"
                //           << "Mterm2: " << m["Mterm2"] << " \n"
                //           << "Cterm1: " << m["Cterm1"] << " \n"
                //           << "findex: " << m["findex"] << " \n"

                if((a1.index >= map[cid]["index"]) and (a1.index <= map[cid]["Nterm2"]))
                {
                    subd1 = 0;
                }
                else if((a1.index >= map[cid]["Mterm1"]) and (a1.index <= map[cid]["Mterm2"]))
                {
                    subd1 = 1;
                }
                else if((a1.index >= map[cid]["Cterm1"]) and (a1.index <= map[cid]["findex"]))
                {
                    subd1 = 2;
                }

                if((a2.index >= map[cid]["index"]) and (a2.index <= map[cid]["Nterm2"]))
                {
                    subd2 = 0;
                }
                else if((a2.index >= map[cid]["Mterm1"]) and (a2.index <= map[cid]["Mterm2"]))
                {
                    subd2 = 1;
                }
                else if((a2.index >= map[cid]["Cterm1"]) and (a2.index <= map[cid]["findex"]))
                {
                    subd2 = 2;
                }

                contacts.push_back(boost::tuple<int,int,double,int,int>(
                                       a1.index,
                                       a2.index,
                                       dist,
                                       subd1,
                                       subd2));
                // try
                // {
                //     contacts.push_back(boost::tuple<int,int,double>(a1.index,a2.index,dist));

                // }
                // catch (const std::bad_alloc &contacts)
                // {
                //     std::cout << "Allocation failed for single contact: " << contacts.what() << std::endl;
                // }
            }
            // i2 += 1;

            // }
        }
        // i1 += 1;
        // i2 = 0;
        // i += 1;
    }

    // std::cout << "contacts_size: " << contacts.size() << std::endl;
    // std::cout << "contacts_counted: " << total_contacts << std::endl;
    return contacts;

    // if (total_contacts > 0)
    // {
    //     return contacts;
    // }
    // else
    // {
    //     // contacts.push_back(boost::tuple<int,int,double>(-1,-1,0));
    //     return std::vector<boost::tuple<int,int,double>>(1,(-1,-1,0));
    // }

}


// void load_dcd_to_atoms(dcdhandle *dcd,Atom *aa);

// std::vector<boost::tuple<int,int,int,double>> get_contacts_for_chain_later(Atom *alater,
//                                                                            double cutoff,
//                                                                            double tolerance,
//                                                                            std::vector<boost::tuple
//                                                                            <int,int,int,double>> contacts)

SetContacts get_contacts_for_chain_later(Atom *alater,
                                         double cutoff,
                                         double tolerance,
                                         SetContacts contacts)
{
    // std::cout << "Welcome to Contacts - Later!" << std::endl;
    // std::vector<std::vector<boost::tuple<int,int,int,double>>> vec_contacts;
    // std::vector<boost::tuple<int,int,int,double>> cur_contacts;
    SetContacts cur_contacts;


    // std::cout << std::get<0>(contacts[0]); << std::endl;

    Vector p1, p2;
    double dist;
    double odist;
    // int onoff; // onoff; 0/1.


    // Note to self.. fix the distance function to account for x,y,z
    // possibly with template or overloaded function.

    // for(auto a: contacts) // 156 in the 312 monomer case.
    // {
    // std::cout << "size: " << a.size() << std::endl;

    // std::cout << "evaluating " << contacts.size() << " contacts." << std::endl;

    for(auto c: contacts)
    {
        // onoff = -1; // will be set to 0,1 in all cases.
        dist = 0.0;
        odist = 0.0;

        // std::cout << boost::get<0>(c) << std::endl;
        // std::cout << boost::get<1>(c) << std::endl;
        // std::cout << boost::get<2>(c) << std::endl;
        // std::cout << "orig_dist: " << boost::get<3>(c) << std::endl;
        odist = boost::get<2>(c);

        // std::cout << std::endl;
        // std::cout << alater[boost::get<0>(c)].x << std::endl;
        // std::cout << alater[boost::get<0>(c)].y << std::endl;
        // std::cout << alater[boost::get<0>(c)].z << std::endl;

        p1.x = alater[boost::get<0>(c)].x;
        p1.y = alater[boost::get<0>(c)].y;
        p1.z = alater[boost::get<0>(c)].z;

        // std::cout << std::endl;
        // std::cout << alater[boost::get<1>(c)].x << std::endl;
        // std::cout << alater[boost::get<1>(c)].y << std::endl;
        // std::cout << alater[boost::get<1>(c)].z << std::endl;

        p2.x = alater[boost::get<1>(c)].x;
        p2.y = alater[boost::get<1>(c)].y;
        p2.z = alater[boost::get<1>(c)].z;


        dist = distance(p1,p2);

        // std::cout << "cur_dist: " << dist << std::endl;
        // std::cout << std::endl;

        if((dist < cutoff) or (dist < odist + tolerance))
        {
            // onoff = 1;
            // Contact:
            // tuple: <index,index, 0 or 1, current-distance>
            // cur_contacts.push_back(boost::tuple<int,int,int,double>(boost::get<0>(c),
            //                                                         boost::get<1>(c),
            //                                                         onoff,
            //                                                         dist));
            cur_contacts.push_back(boost::tuple<int,int,double,int,int>(
                                       boost::get<0>(c),
                                       boost::get<1>(c),
                                       dist,
                                       boost::get<3>(c),
                                       boost::get<4>(c)));

        }
        // else
        // {
        //     onoff = 0;
        // }


        // Contact:
        // tuple: <index,index, 0 or 1, current-distance>
        // cur_contacts.push_back(boost::tuple<int,int,int,double>(boost::get<0>(c),
        //                                                         boost::get<1>(c),
        //                                                         onoff,
        //                                                         dist));

        // cur_contacts.clear();
    }
    // break;

    return cur_contacts;
    // }
}

std::vector<boost::tuple<int,int,int,double>> output_contacts(std::vector<std::vector<boost::tuple
                                                              <int,int,int,double>>> contacts)
{
    std::cout << "Now printing contact file." << std::endl;
    std::cout << "1st-dimension: " << contacts.size() << std::endl;


    // FILE
    FILE * fp_contacts;
    fp_contacts = fopen("emol_contacts.dat","a+");

    for(auto cl: contacts)
    {
        // std::cout << cl.size() << std::endl;

        for(auto c: cl)
        {


        }
        fprintf(fp_contacts,"%d ",cl.size());
        // fprintf(fp_contacts,"%d ");
    }

    fclose(fp_contacts);
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


    // FILE
    FILE * fp_contacts;
    fp_contacts = fopen("emol_mtcontacts_by_subdomain.dat", "w+");
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

                for(auto contact: n)
                {
                    // std::cout << contact.get<3>() << " " << contact.get<4>() << std::endl;
                    if((contact.get<3>() == 0) and (contact.get<4>() == 0))
                    {
                        nn += 1;
                    }
                    else if((contact.get<3>() == 0) and (contact.get<4>() == 1))
                    {
                        nm += 1;
                    }
                    else if((contact.get<3>() == 0) and (contact.get<4>() == 2))
                    {
                        nc += 1;
                    }
                    else if((contact.get<3>() == 1) and (contact.get<4>() == 0))
                    {
                        mn += 1;
                    }
                    else if((contact.get<3>() == 1) and (contact.get<4>() == 1))
                    {
                        mm += 1;
                    }
                    else if((contact.get<3>() == 1) and (contact.get<4>() == 2))
                    {
                        mc += 1;
                    }
                    else if((contact.get<3>() == 2) and (contact.get<4>() == 0))
                    {
                        cn += 1;
                    }
                    else if((contact.get<3>() == 2) and (contact.get<4>() == 1))
                    {
                        cm += 1;
                    }
                    else if((contact.get<3>() == 2) and (contact.get<4>() == 2))
                    {
                        cc += 1;
                    }

                }
                fprintf(fp_contacts,"%5d %5d %5d %5d %5d %5d %5d %5d %5d %d\n",
                        n.size(),
                        nn,nm,nc,
                        mn,mm,mc,
                        cn,cm,cc);

            }
            // fprintf(fp_contacts,"\n");
        }
        fprintf(fp_contacts,"# %d\n",count_f);
    }
    fclose(fp_contacts);
}


SetGlobalContacts explore_global_contacts(SetGlobalContacts gc,
                                          MtIndexMap map,
                                          MtNeighbors mt_matrix)
{
    std::cout << "Exploring Global Contacts .." << std::endl;

    SetGlobalContacts gc3;
    SetChains chain_set3;
    SetNeighbors neighbor_set;

    // Contact Sets.
    SetContacts set_NN;
    SetContacts set_MM;
    SetContacts set_CC;

    SetContacts set_NM;
    SetContacts set_NC;

    SetContacts set_MN;
    SetContacts set_MC;

    SetContacts set_CN;
    SetContacts set_CM;


    int count_d = 0; // dimer
    int count_t = 0; // type: 0, 1, 0-1, 0-2, 0-3, 0-4, .. 9 tot.
    int count_n = 0;

    for(auto f: gc) // frame in global contact array
    {

        count_d = 0;
        for(auto d: f) // chain (but actually dimer) in frame, 156
        {


            count_n = 0;
            for(auto n: d) // 9 situations of 6 neighbors, 0,1,0-1; 0-2,3,4; 1-5,6,7
            {
                // std::cout << n.get<0> << std::endl;
                // std::cout << std::get<0>(n) << std::endl;

                std::cout << "the_alpha: " << count_d * 2 << std::endl;
                std::cout << "the_beta: " << count_d * 2 + 1 << std::endl;
                std::cout << "count_n: " << count_n << std::endl;
                std::cout << "neighbor .. " << mt_matrix[count_d][count_n] << std::endl;

                count_t = 0;
                for(auto t: n)
                {
                    // std::cout << std::get<0>(t) << std::endl;
                    std::cout << t.get<0>() << "  "
                              << t.get<1>() << "  "
                              << t.get<2>() << "  "
                              << t.get<3>() << "  "
                              << t.get<4>() << "  "
                              << std::endl;

                    // Assign N, M, C.
                    // 9 ifs?
                    // 0 .. 8

                    // 0: alpha-alpha
                    // 1: beta-beta
                    // 2: alpha-beta
                    // 345: alpha-234
                    // 678: beta-567





                    if(count_t >= 24)
                    {
                        break;
                    }
                    count_t += 1;
                }

                count_n += 1;
            }

            if (count_d > 12)
            {
                break;
            }

            count_d += 1;
        }
    }

    return gc3;
}


MtNeighbors get_map_of_mtneighbors(std::vector<std::vector <Atom>> chain_ref,
                                   DimerList dimers)
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

        // ic += 1;
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
