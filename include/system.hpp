// foo.h
#ifndef _system_hpp_
#define _system_hpp_


// libraries:
// #include <stdio.h>
/* #include <stdlib.h> */
#include <string>
#include <iostream>
// #include<vector>


// headers:
#include "debug.h"
// #include "radius.hpp"
// #include "pdbfile.hpp"

// Definitions:
/* #define BUFFERSIZE 900 */

using namespace std;

/* ---------------------------------------------------------
   Class declarations
   --------------------------------------------------------- */

/* ---------------------------------------------------------
   System: (Chain,Residue,Atom)
   --------------------------------------------------------- */
class System
{
public:
    System();
    ~System();
    void print_prop();


    int num_chains;
    int num_residues;
    int num_atoms;


    float minx, miny, minz;
    float maxx, maxy, maxz;

private:


};
inline System::System()
{
    // cout << "System construction commencing." << endl;
    // debug("System construction commencing.\n");
    num_chains = -1;
    num_residues = -1;
    num_atoms = -1;

    minx = miny = minz = -999999.0;
    maxx = maxy = maxz =  999999.0;
}
inline System::~System()
{
    // cout << "System destruction commencing." << endl;
    // debug("System destruction commencing.\n");

}
inline void System::print_prop()
{
    cout << "SYSTEM:"
         << "\n\tChains: " << num_chains
         << "\n\tResidues: " << num_residues
         << "\n\tAtoms: " << num_atoms
         << endl;
    cout << "SIZE:"
         << "\n\tX: " << minx << " " << maxx
         << "\n\tY: " << miny << " " << maxy
         << "\n\tZ: " << minz << " " << maxz
         << endl;
}

/* ---------------------------------------------------------
   Chain
   --------------------------------------------------------- */
class Chain: public System
{
public:
    Chain();
    ~Chain();
    void print_prop();

    string chain; // by pdb (pdbfile.cpp, 2nd pass)
    int num_residues;

private:
    int chainid;

};
inline Chain::Chain()
{
    // cout << "Chain construction commencing." << endl;
    // debug("Chain construction commencing.\n");
    chainid = -1;
    num_residues = -1;
    chain = "zz";
}
inline Chain::~Chain()
{
    // cout << "Chain destruction commencing." << endl;
    // debug("Chain destruction commencing.\n");

}
inline void Chain::print_prop()
{
    cout << "Chain: "
         << "\n\t" << chainid << " has " << num_residues
         << "." << endl;

}

/* ---------------------------------------------------------
   Residue
   --------------------------------------------------------- */
class Residue: public Chain
{
public:
    Residue();
    ~Residue();
    void print_prop();

    int id_global;
    int id_local;
    int resid;
    int num_atoms_res;
    string restype;   // HIS, GLU, ILE, LEU ..

private:

};
inline Residue::Residue()
{
    // cout << "Residue construction commencing." << endl;
    // debug("Residue construction commencing.\n");
    id_global = -1;
    id_local = -1;
    resid = -1;
    num_atoms_res = -1;
}
inline Residue::~Residue()
{
    // cout << "Residue destruction commencing." << endl;
    // debug("Residue destruction commencing.\n");

}
inline void Residue::print_prop()
{
    cout << "Residue: "
         << "\n\tglobal,local,resid: "
         << id_global << " " << id_local << " " << resid << endl;

}


/* ---------------------------------------------------------
   Atom
   --------------------------------------------------------- */
class Atom: public Residue
{
public:
    Atom();
    ~Atom();
    void print_coords();

    float x, y, z;
    int index;


    float radius;     // VDW, Richards, Gerstein.

    string atomtype;  // N, CA, C, O, CB, ..

    string general_atomtype; // 77,1

private:
    /* std::string str_index = line.substr(5,6);      // int */
    /* std::string str_atomtype = line.substr(13,3); */
    /* std::string str_restype = line.substr(17,3); */
    /* std::string str_chain = line.substr(21,1); */
    /* std::string str_resid = line.substr(23,3);     // int */
    /* std::string str_x = line.substr(30,8);         // float */
    /* std::string str_y = line.substr(38,8);         // float */
    /* std::string str_z = line.substr(46,8);         // float */

};
inline Atom::Atom()
{
    // cout << "Atom construction commencing." << endl;
    // debug("Atom construction commencing.\n");
    x = 0.00000;
    y = 0.00000;
    z = 0.00000;
    radius = 0.0000;
    index = -1;

}
inline Atom::~Atom()
{
    // cout << "Atom destruction commencing." << endl;
    // debug("Atom destruction commencing.\n");

}
inline void Atom::print_coords()
{
    printf("ATOM %6d  %3s %3s          %8.3f%8.3f%8.3f",
           index,atomtype.c_str(),
           restype.c_str(),x,y,z);
    printf("                       %s\n",general_atomtype.c_str());

}


// inline void Atom::print_info()
// {
//     // std::cout << restype << ": " << atomtype << " " << general_atomtype << std::endl;
//     printf("%s: %s (%s)\n",restype.c_str(),atomtype.c_str(),general_atomtype.c_str());
// }


/* ---------------------------------------------------------
   function declarations
   --------------------------------------------------------- */
// void ReadPDBfile(PDBfile *pdbfile,char filename[40]);
// int select(Atom *aa,char criterion[40],int num);
int system_select(Atom *aa,char const *criterion,int total);
void system_select(Atom *aa,char const *criterion,int total,Atom *asel);


// int select(Atom aa,Atom asel,string param,string criterion,int num);
void get_minmax(System sys);


#endif
