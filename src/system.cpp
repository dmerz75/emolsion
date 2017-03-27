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
// vpAtoms select(vAtoms aa,char const *criterion)
IndexGroup select(vAtoms aa,char const *criterion)
{
    // std::cout << "Selecting .." << std::endl;
    std::string selection(criterion);
    // std::cout << criterion << std::endl;

    // vpAtoms asel;
    IndexGroup isel; // Indices selected.

    int num;
    num = 0;

    // CHAIN
    std::string str1("chain "); // with space, or it will need development...
    std::size_t found1 = selection.find(str1);
    if (found1 != std::string::npos)
    {
        std::string sel_chain = selection.replace(found1,str1.length(),"");
        // std::cout << "chain: " << sel_chain << std::endl;

        for(int i=0; i < aa.size(); i++)
        {
            if(aa[i].chain.compare(sel_chain) == 0)
            {
                // asel[num] = aa[i];
                // asel.push_back(&aa[i]);
                isel.push_back(aa[i].index);
                num += 1;
            }
        }
    }

    // RESID
    std::string str2("resid ");
    std::size_t found2 = selection.find(str2);
    if (found2 != std::string::npos) // found "resid"
    {
        // MUST HAVE "resid .."
        std::string sel_resid0 = selection.replace(found2,str2.length(),"");
        std::string strto("to");
        std::size_t foundto = sel_resid0.find("to");
        // std::cout << "sel_resid0:" << sel_resid0 << std::endl;

        if(foundto != std::string::npos) // found "to"
        {
            // MUST HAVE "resid 0 to 5"
            std::string sel_resid1 = sel_resid0.substr(0,foundto);
            std::string sel_resid2 = sel_resid0.substr(foundto+strto.length(),sel_resid0.length());
            int r1,r2;
            r1 = r2 = -1;
            r1 = stoi(sel_resid1);
            r2 = stoi(sel_resid2);
            // std::cout << "sel_resid1:" << sel_resid1 << '\t' << r1 << std::endl;
            // std::cout << "sel_resid2:" << sel_resid2 << '\t' << r2 << std::endl;

            for(int i=0; i < aa.size(); i++)
            {
                if((aa[i].resid >= r1) and (aa[i].resid <= r2))
                {
                    // asel[num] = aa[i];
                    // asel.push_back(&aa[i]);
                    isel.push_back(aa[i].index);
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
            // std::cout << "sel_resid0:" << sel_resid0 << '\t' << r1 << std::endl;

            for(int i=0; i < aa.size(); i++)
            {

                if(aa[i].resid == r1)
                {
                    // asel[num] = aa[i];
                    // asel.push_back(&aa[i]);
                    isel.push_back(aa[i].index);
                    num += 1;
                }
            }
        }
    }

    // INDEX
    std::string str3("index ");
    std::size_t found3 = selection.find(str3);

    if (found3 != std::string::npos) // found "index"
    {
        // MUST HAVE "index .."
        std::string sel_index0 = selection.replace(found3,str3.length(),"");
        std::string strto("to");
        std::size_t foundto = sel_index0.find("to");
        // std::cout << "sel_index0:" << sel_index0 << std::endl;

        if(foundto != std::string::npos) // found "to"
        {
            // MUST HAVE "index 0 to 5"
            std::string sel_index1 = sel_index0.substr(0,foundto);
            std::string sel_index2 = sel_index0.substr(foundto+strto.length(),sel_index0.length());
            int r1,r2;
            r1 = r2 = -1;
            r1 = stoi(sel_index1);
            r2 = stoi(sel_index2);
            // std::cout << "sel_index1:" << sel_index1 << '\t' << r1 << std::endl;
            // std::cout << "sel_index2:" << sel_index2 << '\t' << r2 << std::endl;

            for(int i=0; i < aa.size(); i++)
            {
                if((aa[i].index >= r1) and (aa[i].index <= r2))
                {
                    // asel[num] = aa[i];
                    // asel.push_back(&aa[i]);
                    isel.push_back(aa[i].index);
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
            // std::cout << "sel_index0:" << sel_index0 << '\t' << r1 << std::endl;

            for(int i=0; i < aa.size(); i++)
            {
                if(aa[i].index == r1)
                {
                    // asel[num] = aa[i];
                    // asel.push_back(&aa[i]);
                    isel.push_back(aa[i].index);
                    num += 1;
                }
            }
        }
    }

    // ALL
    std::string str4("all");
    std::size_t found4 = selection.find(str4);
    if (found4 != std::string::npos) // found "all"
    {
        for(int i=0; i < aa.size(); i++)
        {
            // asel[num] = aa[i];
            // asel.push_back(&aa[i]);
            isel.push_back(aa[i].index);
            num += 1;
        }
    }


    // CHAINID
    std::string str5("chainid ");
    std::size_t found5 = selection.find(str5);
    if (found5 != std::string::npos) // found "chainid"
    {
        // MUST HAVE "chainid .."
        std::string sel_chainid0 = selection.replace(found5,str5.length(),"");
        std::string strto("to");
        std::size_t foundto = sel_chainid0.find("to");
        // std::cout << "sel_chainid0:" << sel_chainid0 << std::endl;

        if(foundto != std::string::npos) // found "to"
        {
            // MUST HAVE "chainid 0 to 5"
            std::string sel_chainid1 = sel_chainid0.substr(0,foundto);
            std::string sel_chainid2 = sel_chainid0.substr(foundto+strto.length(),sel_chainid0.length());
            int cid1,cid2;
            cid1 = cid2 = -1;

            cid1 = stoi(sel_chainid1);
            cid2 = stoi(sel_chainid2);
            // std::cout << "sel_chainid1:" << sel_chainid1 << '\t' << cid1 << std::endl;
            // std::cout << "sel_chainid2:" << sel_chainid2 << '\t' << cid2 << std::endl;


            for(int i=0; i < aa.size(); i++)
            {
                if((aa[i].chainid >= cid1) and (aa[i].chainid <= cid2))
                {
                    // asel[num] = aa[i];
                    // asel.push_back(&aa[i]);
                    isel.push_back(aa[i].index);
                    num += 1;
                }
            }
        }
        else
        {
            // THIS IS CASE: chainid 14
            int cid1;
            cid1 = -1;
            cid1 = stoi(sel_chainid0);
            // std::cout << "sel_chainid0:" << sel_chainid0 << '\t' << cid1 << std::endl;

            for(int i=0; i < aa.size(); i++)
            {

                if(aa[i].chainid == cid1)
                {
                    // asel[num] = aa[i];
                    // asel.push_back(&aa[i]);
                    isel.push_back(aa[i].index);
                    num += 1;
                }
            }
        }
    }

    // Atomtype
    std::string str6("atomtype ");
    std::size_t found6 = selection.find(str6);
    if (found6 != std::string::npos) // found "atomtype"
    {
        // MUST HAVE "atomtype .."
        std::string sel_atomtype0 = selection.replace(found6,str6.length(),"");

        for(int i=0; i < aa.size(); i++)
        {
            // std::cout << aa[i].atomtype;
            // std::cout << " " << sel_atomtype0 << std::endl;

            if(aa[i].atomtype.compare(sel_atomtype0) == 0)
            {
                // asel[num] = aa[i];
                // asel.push_back(&aa[i]);
                isel.push_back(aa[i].index);
                num += 1;
            }
        }
    }

    // return asel;
    return isel;
}

int system_select(Atom *aa,char const *criterion,int total)
{
    // 1.

    // std::cout << "making selection .." << std::endl;
    std::string selection(criterion);
    // std::cout << criterion << std::endl;
    // std::cout << selection << std::endl;

    // printf("num_atoms: %d\n",aa[0].num_atoms);
    int num_all = aa[0].num_atoms;
    int num;
    num = 0;


    std::string str1("chain "); // with space, or it will need development...
    std::size_t found1 = selection.find(str1);


    if (found1 != std::string::npos)
    {
        // std::cout << "got it!" << std::endl;
        std::string sel_chain = selection.replace(found1,str1.length(),"");
        // std::cout << "chain: " << sel_chain << std::endl;

        for(int i=0; i < num_all; i++)
        {
            // std::cout << i << std::endl;
            // std::cout << aa[i].chain << std::endl;

            // DEBUG: had to remove whitespace.
            // std::cout << "chain: " << aa[i].chain.length() << std::endl;
            // std::cout << "sel: " << sel_chain.length() << std::endl;


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
    //     std::cout << "no chain!" << std::endl;
    // }



    std::string str2("resid ");
    std::size_t found2 = selection.find(str2);

    if (found2 != std::string::npos) // found "resid"
    {
        // MUST HAVE "resid .."
        std::string sel_resid0 = selection.replace(found2,str2.length(),"");
        std::string strto("to");
        std::size_t foundto = sel_resid0.find("to");

        // std::cout << "sel_resid0:" << sel_resid0 << std::endl;

        if(foundto != std::string::npos) // found "to"
        {
            // MUST HAVE "resid 0 to 5"
            std::string sel_resid1 = sel_resid0.substr(0,foundto);
            std::string sel_resid2 = sel_resid0.substr(foundto+strto.length(),sel_resid0.length());
            int r1,r2;
            r1 = r2 = -1;

            r1 = stoi(sel_resid1);
            r2 = stoi(sel_resid2);

            // std::cout << "sel_resid1:" << sel_resid1 << '\t' << r1 << std::endl;
            // std::cout << "sel_resid2:" << sel_resid2 << '\t' << r2 << std::endl;


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
            // std::cout << "sel_resid0:" << sel_resid0 << '\t' << r1 << std::endl;

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

        // std::cout << "sel_index0:" << sel_index0 << std::endl;

        if(foundto != std::string::npos) // found "to"
        {
            // MUST HAVE "index 0 to 5"
            std::string sel_index1 = sel_index0.substr(0,foundto);
            std::string sel_index2 = sel_index0.substr(foundto+strto.length(),sel_index0.length());
            int r1,r2;
            r1 = r2 = -1;

            r1 = stoi(sel_index1);
            r2 = stoi(sel_index2);

            // std::cout << "sel_index1:" << sel_index1 << '\t' << r1 << std::endl;
            // std::cout << "sel_index2:" << sel_index2 << '\t' << r2 << std::endl;


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
            // std::cout << "sel_index0:" << sel_index0 << '\t' << r1 << std::endl;

            for(int i=0; i < num_all; i++)
            {

                if(aa[i].index == r1)
                {
                    num += 1;
                }
            }
        }
    }


    std::string str4("all");
    std::size_t found4 = selection.find(str4);

    if (found4 != std::string::npos) // found "all"
    {
        total = num_all;
    }


    std::string str5("chainid ");
    std::size_t found5 = selection.find(str5);

    if (found5 != std::string::npos) // found "chainid"
    {
        // MUST HAVE "chainid .."
        std::string sel_chainid0 = selection.replace(found5,str5.length(),"");
        std::string strto("to");
        std::size_t foundto = sel_chainid0.find("to");

        // std::cout << "sel_chainid0:" << sel_chainid0 << std::endl;

        if(foundto != std::string::npos) // found "to"
        {
            // MUST HAVE "chainid 0 to 5"
            std::string sel_chainid1 = sel_chainid0.substr(0,foundto);
            std::string sel_chainid2 = sel_chainid0.substr(foundto+strto.length(),sel_chainid0.length());
            int cid1,cid2;
            cid1 = cid2 = -1;

            cid1 = stoi(sel_chainid1);
            cid2 = stoi(sel_chainid2);

            // std::cout << "sel_chainid1:" << sel_chainid1 << '\t' << cid1 << std::endl;
            // std::cout << "sel_chainid2:" << sel_chainid2 << '\t' << cid2 << std::endl;


            for(int i=0; i < num_all; i++)
            {
                if((aa[i].chainid >= cid1) and (aa[i].chainid <= cid2))
                {
                    num += 1;
                }
            }
        }
        else
        {
            // THIS IS CASE: chainid 14
            int cid1;
            cid1 = -1;
            cid1 = stoi(sel_chainid0);
            // std::cout << "sel_chainid0:" << sel_chainid0 << '\t' << cid1 << std::endl;

            for(int i=0; i < num_all; i++)
            {

                if(aa[i].chainid == cid1)
                {
                    num += 1;
                }
            }
        }
        total = num;
        return total;
    }


    total = num;
    return total;
}

void system_select(Atom *aa,char const *criterion,int total,Atom *asel)
{
    // Round 2. Get the selection.

    // asel =

    // std::cout << "making selection .." << std::endl;
    std::string selection(criterion);
    // std::cout << criterion << std::endl;
    // std::cout << selection << std::endl;

    // printf("num_atoms: %d\n",aa[0].num_atoms);
    int num_all = aa[0].num_atoms;
    int num;
    num = 0;


    std::string str1("chain "); // with space, or it will need development...
    std::size_t found1 = selection.find(str1);


    if (found1 != std::string::npos)
    {
        // std::cout << "got it!" << std::endl;
        std::string sel_chain = selection.replace(found1,str1.length(),"");
        // std::cout << "chain: " << sel_chain << std::endl;

        for(int i=0; i < num_all; i++)
        {
            // std::cout << i << std::endl;
            // std::cout << aa[i].chain << std::endl;

            // DEBUG: had to remove whitespace.
            // std::cout << "chain: " << aa[i].chain.length() << std::endl;
            // std::cout << "sel: " << sel_chain.length() << std::endl;


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
    //     std::cout << "no chain!" << std::endl;
    // }



    std::string str2("resid ");
    std::size_t found2 = selection.find(str2);

    if (found2 != std::string::npos) // found "resid"
    {
        // MUST HAVE "resid .."
        std::string sel_resid0 = selection.replace(found2,str2.length(),"");
        std::string strto("to");
        std::size_t foundto = sel_resid0.find("to");

        // std::cout << "sel_resid0:" << sel_resid0 << std::endl;

        if(foundto != std::string::npos) // found "to"
        {
            // MUST HAVE "resid 0 to 5"
            std::string sel_resid1 = sel_resid0.substr(0,foundto);
            std::string sel_resid2 = sel_resid0.substr(foundto+strto.length(),sel_resid0.length());
            int r1,r2;
            r1 = r2 = -1;

            r1 = stoi(sel_resid1);
            r2 = stoi(sel_resid2);

            // std::cout << "sel_resid1:" << sel_resid1 << '\t' << r1 << std::endl;
            // std::cout << "sel_resid2:" << sel_resid2 << '\t' << r2 << std::endl;


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
            // std::cout << "sel_resid0:" << sel_resid0 << '\t' << r1 << std::endl;

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

        // std::cout << "sel_index0:" << sel_index0 << std::endl;

        if(foundto != std::string::npos) // found "to"
        {
            // MUST HAVE "index 0 to 5"
            std::string sel_index1 = sel_index0.substr(0,foundto);
            std::string sel_index2 = sel_index0.substr(foundto+strto.length(),sel_index0.length());
            int r1,r2;
            r1 = r2 = -1;

            r1 = stoi(sel_index1);
            r2 = stoi(sel_index2);

            // std::cout << "sel_index1:" << sel_index1 << '\t' << r1 << std::endl;
            // std::cout << "sel_index2:" << sel_index2 << '\t' << r2 << std::endl;


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
            // std::cout << "sel_index0:" << sel_index0 << '\t' << r1 << std::endl;

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

    std::string str4("all");
    std::size_t found4 = selection.find(str4);

    if (found4 != std::string::npos) // found "all"
    {
        for(int i=0; i < num_all; i++)
        {
            asel[num] = aa[i];
            num += 1;
        }
    }

    // total = num;
    // return total;



    std::string str5("chainid ");
    std::size_t found5 = selection.find(str5);

    if (found5 != std::string::npos) // found "chainid"
    {
        // MUST HAVE "chainid .."
        std::string sel_chainid0 = selection.replace(found5,str5.length(),"");
        std::string strto("to");
        std::size_t foundto = sel_chainid0.find("to");

        // std::cout << "sel_chainid0:" << sel_chainid0 << std::endl;

        if(foundto != std::string::npos) // found "to"
        {
            // MUST HAVE "chainid 0 to 5"
            std::string sel_chainid1 = sel_chainid0.substr(0,foundto);
            std::string sel_chainid2 = sel_chainid0.substr(foundto+strto.length(),sel_chainid0.length());
            int cid1,cid2;
            cid1 = cid2 = -1;

            cid1 = stoi(sel_chainid1);
            cid2 = stoi(sel_chainid2);

            // std::cout << "sel_chainid1:" << sel_chainid1 << '\t' << cid1 << std::endl;
            // std::cout << "sel_chainid2:" << sel_chainid2 << '\t' << cid2 << std::endl;


            for(int i=0; i < num_all; i++)
            {
                if((aa[i].chainid >= cid1) and (aa[i].chainid <= cid2))
                {
                    asel[num] = aa[i];
                    num += 1;
                }
            }
        }
        else
        {
            // THIS IS CASE: chainid 14
            int cid1;
            cid1 = -1;
            cid1 = stoi(sel_chainid0);
            // std::cout << "sel_chainid0:" << sel_chainid0 << '\t' << cid1 << std::endl;

            for(int i=0; i < num_all; i++)
            {

                if(aa[i].chainid == cid1)
                {
                    asel[num] = aa[i];
                    num += 1;
                }
            }
        }
        total = num;
        // return total;
    }


    total = num;
    // return total;
}

int system_select_atomtype(Atom *aa,char const *criterion,int total,Atom *asel)
{
    // Get the selection.
    // criterion: N, C, CA, O
    // the length of the search, CA=2, HA1=3.
    // std::cout << "the criterion is: " << criterion << std::endl;
    int crit_length = 0;
    crit_length = strlen(criterion);
    // std::cout << "criterion_length: " << crit_length << std::endl;

    std::string selection(criterion);
    int num_all = aa[0].num_atoms;
    int num;
    num = 0;


    // std::string str1("chain "); // with space, or it will need development...
    // std::size_t found1 = selection.find(str1);
    // std::size_t found1 = selection.find();
    std::string atomtype;

    if(strcmp(criterion,"backbone") == 0)
    {
        // std::cout << "Found backbone." << std::endl;
        for(int i=0; i < num_all; i++)
        {
            atomtype = "C";
            if(atomtype.length() == aa[i].atomtype.length())
            {
                std::size_t found1 = aa[i].atomtype.find(atomtype);
                if(found1 != std::string::npos)
                {
                    asel[num] = aa[i];
                    num += 1;
                }
            }
            atomtype = "CA";
            if(atomtype.length() == aa[i].atomtype.length())
            {
                std::size_t found1 = aa[i].atomtype.find(atomtype);
                if(found1 != std::string::npos)
                {
                    asel[num] = aa[i];
                    num += 1;
                }
            }
            atomtype = "O";
            if(atomtype.length() == aa[i].atomtype.length())
            {
                std::size_t found1 = aa[i].atomtype.find(atomtype);
                if(found1 != std::string::npos)
                {
                    asel[num] = aa[i];
                    num += 1;
                }
            }
            atomtype = "N";
            if(atomtype.length() == aa[i].atomtype.length())
            {
                std::size_t found1 = aa[i].atomtype.find(atomtype);
                if(found1 != std::string::npos)
                {
                    asel[num] = aa[i];
                    num += 1;
                }
            }
        }
    }
    else
    {
        for(int i=0; i < num_all; i++)
        {
            // if(i < 20)
            // {
            //     std::cout << "count: " << i << " ";
            //     std::cout << "index: " << aa[i].index << " ";
            //     std::cout << "atomtype: " << aa[i].atomtype << std::endl;
            // }
            // if(i >= 20)
            // {
            //     break;
            // }
            // std::cout << "atomtype_length: " << aa[i].atomtype.length() << std::endl;

            if(aa[i].atomtype.length() != crit_length)
            {
                continue;
            }
            else
            {
                std::size_t found1 = aa[i].atomtype.find(criterion);
                if(found1 != std::string::npos)
                {
                    asel[num] = aa[i];
                    num += 1;
                }
            }
        }
    }

    total = num;
    return num;
    // std::cout << "the number of atoms by atomtype was: " << total << std::endl;
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


// std::vector<Atom> system_select(Atom *aa,char const *criterion,int total)
// std::vector<Atom> system_select(Atom *aa,char const *criterion)
// {
//     /*
//       Round 3. Get the selection.
//     */

//     std::vector<Atom> aabin;
//         // std::vector<Atom*> chain1;


//     std::cout << "making selection .." << std::endl;
//     std::string selection(criterion);
//     // std::cout << criterion << std::endl;
//     std::cout << selection << std::endl;

//     // printf("num_atoms: %d\n",aa[0].num_atoms);
//     int num_all = aa[0].num_atoms;
//     int num, total;
//     num = total = 0;


//     std::string str1("chain "); // with space, or it will need development...
//     std::size_t found1 = selection.find(str1);


//     if (found1 != std::string::npos)
//     {
//         // std::cout << "got it!" << std::endl;
//         std::string sel_chain = selection.replace(found1,str1.length(),"");
//         // std::cout << "chain: " << sel_chain << std::endl;

//         for(int i=0; i < num_all; i++)
//         {
//             // std::cout << i << std::endl;
//             // std::cout << aa[i].chain << std::endl;

//             // DEBUG: had to remove whitespace.
//             // std::cout << "chain: " << aa[i].chain.length() << std::endl;
//             // std::cout << "sel: " << sel_chain.length() << std::endl;


//             if(aa[i].chain.compare(sel_chain) == 0)
//             {
//                 // asel[num] = aa[i];
//                 Atom asel = aa[i];
//                 num += 1;
//             }

//             // if(aa[i].chain.compare(sel_chain) == 0)
//             // {
//             //     num += 1;
//             // }
//         }
//         total = num;
//         // return total;
//     }
//     // else
//     // {
//     //     std::cout << "no chain!" << std::endl;
//     // }



//     std::string str2("resid ");
//     std::size_t found2 = selection.find(str2);

//     if (found2 != std::string::npos) // found "resid"
//     {
//         // MUST HAVE "resid .."
//         std::string sel_resid0 = selection.replace(found2,str2.length(),"");
//         std::string strto("to");
//         std::size_t foundto = sel_resid0.find("to");

//         // std::cout << "sel_resid0:" << sel_resid0 << std::endl;

//         if(foundto != std::string::npos) // found "to"
//         {
//             // MUST HAVE "resid 0 to 5"
//             std::string sel_resid1 = sel_resid0.substr(0,foundto);
//             std::string sel_resid2 = sel_resid0.substr(foundto+strto.length(),sel_resid0.length());
//             int r1,r2;
//             r1 = r2 = -1;

//             r1 = stoi(sel_resid1);
//             r2 = stoi(sel_resid2);

//             // std::cout << "sel_resid1:" << sel_resid1 << '\t' << r1 << std::endl;
//             // std::cout << "sel_resid2:" << sel_resid2 << '\t' << r2 << std::endl;


//             for(int i=0; i < num_all; i++)
//             {
//                 if((aa[i].resid >= r1) and (aa[i].resid <= r2))
//                 {
//                     // asel[num] = aa[i];
//                     Atom asel = aa[i];
//                     num += 1;
//                 }
//             }
//         }
//         else
//         {
//             // THIS IS CASE: resid 14
//             int r1;
//             r1 = -1;
//             r1 = stoi(sel_resid0);
//             // std::cout << "sel_resid0:" << sel_resid0 << '\t' << r1 << std::endl;

//             for(int i=0; i < num_all; i++)
//             {

//                 if(aa[i].resid == r1)
//                 {
//                     // asel[num] = aa[i];
//                     Atom asel = aa[i];
//                     num += 1;
//                 }
//             }
//         }
//         total = num;
//         // return total;
//     }



//     std::string str3("index ");
//     std::size_t found3 = selection.find(str3);

//     if (found3 != std::string::npos) // found "index"
//     {
//         // MUST HAVE "index .."
//         std::string sel_index0 = selection.replace(found3,str3.length(),"");
//         std::string strto("to");
//         std::size_t foundto = sel_index0.find("to");

//         // std::cout << "sel_index0:" << sel_index0 << std::endl;

//         if(foundto != std::string::npos) // found "to"
//         {
//             // MUST HAVE "index 0 to 5"
//             std::string sel_index1 = sel_index0.substr(0,foundto);
//             std::string sel_index2 = sel_index0.substr(foundto+strto.length(),sel_index0.length());
//             int r1,r2;
//             r1 = r2 = -1;

//             r1 = stoi(sel_index1);
//             r2 = stoi(sel_index2);

//             // std::cout << "sel_index1:" << sel_index1 << '\t' << r1 << std::endl;
//             // std::cout << "sel_index2:" << sel_index2 << '\t' << r2 << std::endl;


//             for(int i=0; i < num_all; i++)
//             {
//                 if((aa[i].index >= r1) and (aa[i].index <= r2))
//                 {
//                     // asel[num] = aa[i];
//                     Atom asel = aa[i];
//                     num += 1;
//                 }
//             }
//         }
//         else
//         {
//             // THIS IS CASE: index 14
//             int r1;
//             r1 = -1;
//             r1 = stoi(sel_index0);
//             // std::cout << "sel_index0:" << sel_index0 << '\t' << r1 << std::endl;

//             for(int i=0; i < num_all; i++)
//             {

//                 if(aa[i].index == r1)
//                 {
//                     // asel[num] = aa[i];
//                     Atom asel = aa[i];
//                     num += 1;
//                 }
//             }
//         }
//     }

//     std::string str4("all");
//     std::size_t found4 = selection.find(str4);

//     if (found4 != std::string::npos) // found "all"
//     {
//         for(int i=0; i < num_all; i++)
//         {
//             // asel[num] = aa[i];
//             Atom asel = aa[i];
//             num += 1;
//         }
//     }

//     total = num;
//     // return asel;
//     // return total;
// }

Vector get_centroid(IndexGroup ig,vAtoms aa)
{
    // OVERLOADED
    // std::cout << "Welcome to get_centroid!" << std::endl;
    // std::cout << a.size() << std::endl;

    Vector v;
    float x,y,z;
    x = y = z = 0.0;

    // for(auto a: aa)
    // {
    //     x += a->x;
    //     y += a->y;
    //     z += a->z;
    // }

    for(auto i:ig)
    {
        x += aa[i].x;
        y += aa[i].y;
        z += aa[i].z;
    }

    // v.x = x / aa.size();
    // v.y = y / aa.size();
    // v.z = z / aa.size();

    v.x = x / ig.size();
    v.y = y / ig.size();
    v.z = z / ig.size();

    return v;
}
Vector get_centroid(vAtoms aa)
{
    // OVERLOADED
    // std::cout << "Welcome to get_centroid!" << std::endl;
    // std::cout << a.size() << std::endl;

    Vector v;
    float x,y,z;
    x = y = z = 0.0;

    for(auto a: aa)
    {
        x += a.x;
        y += a.y;
        z += a.z;
    }

    v.x = x / aa.size();
    v.y = y / aa.size();
    v.z = z / aa.size();

    return v;
}

// vAtoms set_chainid(vAtoms aa)
vAtoms set_chainid(vAtoms aa)
{
    // std::cout << "Getting Segments/chains." << std::endl;
    // std::cout << "Num_atoms: " << aa.size() << std::endl;

    int chainid = 0;

    std::string chain;
    chain = aa[0].chain; // Get "A"

    vAtoms a_ret; // return
    a_ret.reserve(aa.size());

    for(int i=0; i<aa.size(); i++)
    {
        // printf("chain: %s\n",aa[i].chain.c_str());
        // printf("atom: %d  chain: %s  resid: %d  atoms-resid: %d\n",
        //        aa[i].index,aa[i].chain.c_str(),aa[i].resid,
        //        aa[i].num_atoms_res);

        Atom a = aa[i];
        if(aa[i].chain.compare(chain) != 0)
        {
            chainid += 1;
        }

        // Set chainid.
        a.chainid = chainid;
        // std::cout << "chainid: " << a.chainid << std::endl;

        // Reset chain to next in line, "A" -> "B"
        chain = aa[i].chain;
        // Push it back:
        a_ret.push_back(a);
    }

    return a_ret;
    // return a_ret;
    // exit(0);
}

vIndexGroup sort_segment_chain(vAtoms aa)
{

    // int chainid = 0;
    int chainid = aa[0].chainid;
    // vvpAtoms segment;
    // vpAtoms pa;
    // Atom *a;

    IndexGroup segment;
    vIndexGroup segments;


    std::cout << "Sorting segments & chains. "
              << "First chainid: " << chainid
              << std::endl;

    // for(auto a: aa)
    // {
    //     std::cout << "Chainid: " << a.chainid << std::endl;
    // }


    for(int i=0; i<aa.size(); i++)
    {
        // a = *aa[i];

        if(chainid == aa[i].chainid)
        {
            // std::cout << "Pushing back atom." << std::endl;
            // pa.push_back(&aa[i]);
            // pa.push_back(a);
            segment.push_back(aa[i].index);
        }
        else
        {
            // std::cout << "Pushing back segment. "
            //           << chainid
            //           << std::endl;
            // std::cout << "Size: " << pa.size() << std::endl;
            segments.push_back(segment);
            segment.clear();
            chainid = aa[i].chainid;

            // std::cout << "current_point_size: " << pa.size() << std::endl;
            // pa.push_back(&aa[i]);
            // pa.push_back(a);
            segment.push_back(aa[i].index);

            // std::cout << "current_point_size: " << pa.size() << std::endl;
        }
    }

    // std::cout << "Pushing back segment. "
    //           << chainid
    //           << std::endl;
    // std::cout << "Size: " << pa.size() << std::endl;
    // segment.push_back(pa);
    segments.push_back(segment);
    // pa.clear();
    segment.clear();

    std::cout << "Number of chains: " << segments.size() << std::endl;
    return segments;
}
