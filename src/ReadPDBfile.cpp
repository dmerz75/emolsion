// ReadPDBfile.cpp

/* ---------------------------------------------------------
   standard libraries
   --------------------------------------------------------- */
// #include<debug.h>
// #include <stdio.h>
// #include <stdlib.h> // strtod?, stod
// #include <assert.h>
#include <string.h> // string.
// #include <cctype>
// #include <algorithm> // remove_if, count
#include <iostream>
#include <fstream> // ifstream
#include <boost/algorithm/string.hpp>

#include <system.hpp>

// #include <ctime>
// #include <list>        // std::list
// #include <vector>
// #include <iterator> // istream_iterator


/* ---------------------------------------------------------
   headers
   --------------------------------------------------------- */
#include "debug.h"
#include "ReadPDBfile.hpp"


using namespace std;

/* ---------------------------------------------------------
   functions
   --------------------------------------------------------- */
int ReadPDBfile(char filename[40],int total_atoms,Atom *aa)
{
    // char filename
    // string fn
    // string line
    // ifstream pdbfile

    string fn(filename);
    cout << "Reading: " << fn.c_str() << endl;


    string line;
    ifstream pdbfile(fn.c_str());


    // While reading, count the atoms.
    int num_atoms;
    num_atoms = 0;
    // ,total_atoms;
    // total_atoms = -1;


    if (pdbfile.is_open())
    {
        cout << "file is now open." << endl;

        if (total_atoms == -1)
        {
            while(getline(pdbfile,line))
            {
                if (boost::algorithm::starts_with(line,"ATOM") == 1)
                {
                    num_atoms += 1;
                    // Atom a;
                }
            }
            pdbfile.close();
            total_atoms = num_atoms;
            cout << "file is now closed." << endl;

            return total_atoms;
        }
        else // total_atoms != -1.
        {
            while(getline(pdbfile,line))
            {
                if (boost::algorithm::starts_with(line,"ATOM") == 1)
                {

                    // Atom a;

                    /* std::string str_index = line.substr(5,6);      // int */
                    /* std::string str_chain = line.substr(21,1); */
                    /* std::string str_resid = line.substr(23,3);     // int */
                    string str_x = line.substr(30,8);         // float
                    string str_y = line.substr(38,8);         // float
                    string str_z = line.substr(46,8);         // float
                    aa[num_atoms].x = atof(str_x.c_str());
                    aa[num_atoms].y = atof(str_y.c_str());
                    aa[num_atoms].z = atof(str_z.c_str());

                    // printf("index: %d\n",num_atoms);
                    aa[num_atoms].index = num_atoms;
                    aa[num_atoms].atomtype = line.substr(13,3);// CA,CG1,OD1
                    aa[num_atoms].restype = line.substr(17,3); // GLY,ILE

                    try
                    {
                        aa[num_atoms].general_atomtype = line.substr(77,1); // C,N,O,F
                    }
                    catch (const std::out_of_range & e)
                    {
                        // std::cout<< e.what(); // message to console: vector::_M_range_check
                    }



                    // For Residue:
                    aa[num_atoms].resid = atoi(line.substr(23,3).c_str()); // 4,603 for 2KHO

                    // For Chain:
                    aa[num_atoms].chain = line.substr(21,1); // like "A"


                    num_atoms += 1;
                }
            }
            pdbfile.close();
            total_atoms = num_atoms;
            cout << "file is now closed." << endl;
        }
    }
    else
    {
        cout << "Unable to open file";
    }

    cout << "The number of atoms found was: " << num_atoms << endl;
    return total_atoms;
}
