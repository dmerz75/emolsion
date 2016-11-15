// system.cpp

/* ---------------------------------------------------------
   standard libraries
   --------------------------------------------------------- */
// #include <stdio.h>
// #include <stdlib.h> // strtod?, stod
// #include <assert.h>
#include <string> // string. for std::string
// #include <cctype>
// #include <algorithm> // remove_if, count
// #include <iostream>
// #include <fstream>
// #include <ctime>
// #include <list>        // std::list
// #include <vector>
// #include <iterator> // istream_iterator


/* ---------------------------------------------------------
   headers
   --------------------------------------------------------- */
#include "debug.h"
#include "system.hpp"



/* ---------------------------------------------------------
   functions
   --------------------------------------------------------- */
int system_select(Atom *aa,char const *criterion,int total)
{

    // std::cout << "making selection .." << endl;
    std::string selection(criterion);
    // std::cout << criterion << endl;
    // std::cout << selection << endl;

    // printf("num_atoms: %d\n",aa[0].num_atoms);
    int num_all = aa[0].num_atoms;
    int num;
    num = 0;


    std::string str1("chain "); // with space, or it will need development...
    std::size_t found1 = selection.find(str1);


    if (found1 != std::string::npos)
    {
        // std::cout << "got it!" << endl;
        std::string sel_chain = selection.replace(found1,str1.length(),"");
        // std::cout << "chain: " << sel_chain << endl;

        for(int i=0; i < num_all; i++)
        {
            // std::cout << i << endl;
            // std::cout << aa[i].chain << endl;

            // DEBUG: had to remove whitespace.
            // std::cout << "chain: " << aa[i].chain.length() << endl;
            // std::cout << "sel: " << sel_chain.length() << endl;


            if(aa[i].chain.compare(sel_chain) == 0)
            {
                num += 1;
            }

            // if(aa[i].chain.compare(sel_chain) == 0)
            // {
            //     num += 1;
            // }
        }
        total = num;
        return total;
    }
    // else
    // {
    //     std::cout << "no chain!" << endl;
    // }



    std::string str2("resid ");
    std::size_t found2 = selection.find(str2);

    if (found2 != std::string::npos) // found "resid"
    {
        // MUST HAVE "resid .."
        std::string sel_resid0 = selection.replace(found2,str2.length(),"");
        std::string strto("to");
        std::size_t foundto = sel_resid0.find("to");

        // std::cout << "sel_resid0:" << sel_resid0 << endl;

        if(foundto != std::string::npos) // found "to"
        {
            // MUST HAVE "resid 0 to 5"
            std::string sel_resid1 = sel_resid0.substr(0,foundto);
            std::string sel_resid2 = sel_resid0.substr(foundto+strto.length(),sel_resid0.length());
            int r1,r2;
            r1 = r2 = -1;

            r1 = stoi(sel_resid1);
            r2 = stoi(sel_resid2);

            // std::cout << "sel_resid1:" << sel_resid1 << '\t' << r1 << endl;
            // std::cout << "sel_resid2:" << sel_resid2 << '\t' << r2 << endl;


            for(int i=0; i < num_all; i++)
            {
                if((aa[i].resid >= r1) and (aa[i].resid <= r2))
                {
                    num += 1;
                }
            }
        }
        else
        {
            // THIS IS CASE: resid 14
            int r1;
            r1 = -1;
            r1 = stoi(sel_resid0);
            // std::cout << "sel_resid0:" << sel_resid0 << '\t' << r1 << endl;

            for(int i=0; i < num_all; i++)
            {

                if(aa[i].resid == r1)
                {
                    num += 1;
                }
            }
        }
        total = num;
        return total;
    }



    std::string str3("index ");
    std::size_t found3 = selection.find(str3);

    if (found3 != std::string::npos) // found "index"
    {
        // MUST HAVE "index .."
        std::string sel_index0 = selection.replace(found3,str3.length(),"");
        std::string strto("to");
        std::size_t foundto = sel_index0.find("to");

        // std::cout << "sel_index0:" << sel_index0 << endl;

        if(foundto != std::string::npos) // found "to"
        {
            // MUST HAVE "index 0 to 5"
            std::string sel_index1 = sel_index0.substr(0,foundto);
            std::string sel_index2 = sel_index0.substr(foundto+strto.length(),sel_index0.length());
            int r1,r2;
            r1 = r2 = -1;

            r1 = stoi(sel_index1);
            r2 = stoi(sel_index2);

            // std::cout << "sel_index1:" << sel_index1 << '\t' << r1 << endl;
            // std::cout << "sel_index2:" << sel_index2 << '\t' << r2 << endl;


            for(int i=0; i < num_all; i++)
            {
                if((aa[i].index >= r1) and (aa[i].index <= r2))
                {
                    num += 1;
                }
            }
        }
        else
        {
            // THIS IS CASE: index 14
            int r1;
            r1 = -1;
            r1 = stoi(sel_index0);
            // std::cout << "sel_index0:" << sel_index0 << '\t' << r1 << endl;

            for(int i=0; i < num_all; i++)
            {

                if(aa[i].index == r1)
                {
                    num += 1;
                }
            }
        }
    }


    total = num;
    return total;
}

void system_select(Atom *aa,char const *criterion,int total,Atom *asel)
{
    /*
      Round 2. Get the selection.

    */

    // std::cout << "making selection .." << endl;
    std::string selection(criterion);
    // std::cout << criterion << endl;
    // std::cout << selection << endl;

    // printf("num_atoms: %d\n",aa[0].num_atoms);
    int num_all = aa[0].num_atoms;
    int num;
    num = 0;


    std::string str1("chain "); // with space, or it will need development...
    std::size_t found1 = selection.find(str1);


    if (found1 != std::string::npos)
    {
        // std::cout << "got it!" << endl;
        std::string sel_chain = selection.replace(found1,str1.length(),"");
        // std::cout << "chain: " << sel_chain << endl;

        for(int i=0; i < num_all; i++)
        {
            // std::cout << i << endl;
            // std::cout << aa[i].chain << endl;

            // DEBUG: had to remove whitespace.
            // std::cout << "chain: " << aa[i].chain.length() << endl;
            // std::cout << "sel: " << sel_chain.length() << endl;


            if(aa[i].chain.compare(sel_chain) == 0)
            {
                asel[num] = aa[i];
                num += 1;
            }

            // if(aa[i].chain.compare(sel_chain) == 0)
            // {
            //     num += 1;
            // }
        }
        total = num;
        // return total;
    }
    // else
    // {
    //     std::cout << "no chain!" << endl;
    // }



    std::string str2("resid ");
    std::size_t found2 = selection.find(str2);

    if (found2 != std::string::npos) // found "resid"
    {
        // MUST HAVE "resid .."
        std::string sel_resid0 = selection.replace(found2,str2.length(),"");
        std::string strto("to");
        std::size_t foundto = sel_resid0.find("to");

        // std::cout << "sel_resid0:" << sel_resid0 << endl;

        if(foundto != std::string::npos) // found "to"
        {
            // MUST HAVE "resid 0 to 5"
            std::string sel_resid1 = sel_resid0.substr(0,foundto);
            std::string sel_resid2 = sel_resid0.substr(foundto+strto.length(),sel_resid0.length());
            int r1,r2;
            r1 = r2 = -1;

            r1 = stoi(sel_resid1);
            r2 = stoi(sel_resid2);

            // std::cout << "sel_resid1:" << sel_resid1 << '\t' << r1 << endl;
            // std::cout << "sel_resid2:" << sel_resid2 << '\t' << r2 << endl;


            for(int i=0; i < num_all; i++)
            {
                if((aa[i].resid >= r1) and (aa[i].resid <= r2))
                {
                    asel[num] = aa[i];
                    num += 1;
                }
            }
        }
        else
        {
            // THIS IS CASE: resid 14
            int r1;
            r1 = -1;
            r1 = stoi(sel_resid0);
            // std::cout << "sel_resid0:" << sel_resid0 << '\t' << r1 << endl;

            for(int i=0; i < num_all; i++)
            {

                if(aa[i].resid == r1)
                {
                    asel[num] = aa[i];
                    num += 1;
                }
            }
        }
        total = num;
        // return total;
    }



    std::string str3("index ");
    std::size_t found3 = selection.find(str3);

    if (found3 != std::string::npos) // found "index"
    {
        // MUST HAVE "index .."
        std::string sel_index0 = selection.replace(found3,str3.length(),"");
        std::string strto("to");
        std::size_t foundto = sel_index0.find("to");

        // std::cout << "sel_index0:" << sel_index0 << endl;

        if(foundto != std::string::npos) // found "to"
        {
            // MUST HAVE "index 0 to 5"
            std::string sel_index1 = sel_index0.substr(0,foundto);
            std::string sel_index2 = sel_index0.substr(foundto+strto.length(),sel_index0.length());
            int r1,r2;
            r1 = r2 = -1;

            r1 = stoi(sel_index1);
            r2 = stoi(sel_index2);

            // std::cout << "sel_index1:" << sel_index1 << '\t' << r1 << endl;
            // std::cout << "sel_index2:" << sel_index2 << '\t' << r2 << endl;


            for(int i=0; i < num_all; i++)
            {
                if((aa[i].index >= r1) and (aa[i].index <= r2))
                {
                    asel[num] = aa[i];
                    num += 1;
                }
            }
        }
        else
        {
            // THIS IS CASE: index 14
            int r1;
            r1 = -1;
            r1 = stoi(sel_index0);
            // std::cout << "sel_index0:" << sel_index0 << '\t' << r1 << endl;

            for(int i=0; i < num_all; i++)
            {

                if(aa[i].index == r1)
                {
                    asel[num] = aa[i];
                    num += 1;
                }
            }
        }
    }


    total = num;
    // return total;
}

void get_minmax(System sys)
{
    debug("you made it!\n");

    // float minx, miny, minz;
    // float maxx, maxy, maxz;

    // minx = miny = minz = 9999.0;
    // maxx = maxy = maxz = -9999.0;


    // for(int i=0; i<sys.num_segment; i++)
    // {
    //     // printf("%d\n",i);

    //     for(int j=0; j<sys.segment[i].num_residue; j++)
    //     {
    //         // printf("%d ",j);

    //         for(int k=0; k<sys.segment[i].residue[j].num_atom; k++)
    //         {
    //             // sys.segment[i].residue[j].atom[k].print_info();
    //             // sys.segment[i].residue[j].atom[k].print_coords();

    //             if (sys.segment[i].residue[j].atom[k].x < minx)
    //             {
    //                 minx = sys.segment[i].residue[j].atom[k].x;
    //             }
    //             if (sys.segment[i].residue[j].atom[k].y < miny)
    //             {
    //                 miny = sys.segment[i].residue[j].atom[k].y;
    //             }
    //             if (sys.segment[i].residue[j].atom[k].z < minz)
    //             {
    //                 minz = sys.segment[i].residue[j].atom[k].z;
    //             }

    //             if (sys.segment[i].residue[j].atom[k].x > maxx)
    //             {
    //                 maxx = sys.segment[i].residue[j].atom[k].x;
    //             }
    //             if (sys.segment[i].residue[j].atom[k].y > maxy)
    //             {
    //                 maxy = sys.segment[i].residue[j].atom[k].y;
    //             }
    //             if (sys.segment[i].residue[j].atom[k].z > maxz)
    //             {
    //                 maxz = sys.segment[i].residue[j].atom[k].z;
    //             }



    //         }
    //         // printf("\n");
    //     }
    //     // printf("\n");
    // }

    // sys.minx = minx;
    // sys.miny = miny;
    // sys.minz = minz;

    // sys.maxx = maxx;
    // sys.maxy = maxy;
    // sys.maxz = maxz;
}
