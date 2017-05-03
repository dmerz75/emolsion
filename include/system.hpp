// foo.h
#ifndef _system_hpp_
#define _system_hpp_



/* ---------------------------------------------------------
   libraries:
   --------------------------------------------------------- */
// #include <stdio.h>
/* #include <stdlib.h> */
#include <string>
#include <iostream>
#include<vector>



/* ---------------------------------------------------------
   headers:
   --------------------------------------------------------- */
#include "debug.h"
#include "md.h"
// #include "radius.hpp"
// #include "pdbfile.hpp"


/* ---------------------------------------------------------
   Definitions, Typedefs
   --------------------------------------------------------- */


/* ---------------------------------------------------------
   Class declarations
   --------------------------------------------------------- */
   // System: (Chain,Residue,Atom)

class System
{
public:
    System();
    ~System();
    void print_Prop();

    int num_chains;
    int num_residues;
    int num_atoms;


    float minx, miny, minz;
    float maxx, maxy, maxz;

private:


};
inline System::System()
{
    // std::cout << "System construction commencing." << std::endl;
    // debug("System construction commencing.\n");
    num_chains = -1;
    num_residues = -1;
    num_atoms = -1;

    minx = miny = minz = -999999.0;
    maxx = maxy = maxz =  999999.0;
}
inline System::~System()
{
    // std::cout << "System destruction commencing." << std::endl;
    // debug("System destruction commencing.\n");

}
inline void System::print_Prop()
{
    std::cout << "SYSTEM:"
              << "\n\tChains: " << num_chains
              << "\n\tResidues: " << num_residues
              << "\n\tAtoms: " << num_atoms
              << std::endl;
    std::cout << "SIZE:"
              << "\n\tX: " << minx << " " << maxx
              << "\n\tY: " << miny << " " << maxy
              << "\n\tZ: " << minz << " " << maxz
              << std::endl;
}

/* ---------------------------------------------------------
   Chain
   --------------------------------------------------------- */
class Chain: public System
{
public:
    Chain();
    Chain(const Chain &obj);
    ~Chain();
    void print_Prop();

    std::string chain; // by pdb (pdbfile.cpp, 2nd pass)
    int num_residues;
    // int num_atoms;
    int chainid;

private:

};
inline Chain::Chain()
{
    // std::cout << "Chain construction commencing." << std::endl;
    // debug("Chain construction commencing.\n");
    chainid = -1;
    num_residues = -1;
    chain = "zz";
}
inline Chain::Chain(const Chain &obj)
{
    // Copy constructor.
    chain = obj.chain;
    num_residues = obj.num_residues;
    chainid = obj.chainid;
}
inline Chain::~Chain()
{
    // std::cout << "Chain destruction commencing." << std::endl;
    // debug("Chain destruction commencing.\n");

}
inline void Chain::print_Prop()
{
    std::cout << "Chain: "
              << "\n\t" << chainid << " has " << num_residues
              << "." << std::endl;

}

/* ---------------------------------------------------------
   Residue
   --------------------------------------------------------- */
class Residue: public Chain
{
public:
    Residue();
    Residue(const Residue &obj);
    ~Residue();
    void print_Prop();

    int id_global;
    int id_local;
    int resid;
    int num_atoms_res;
    std::string restype;   // HIS, GLU, ILE, LEU ..

private:

};
inline Residue::Residue()
{
    // std::cout << "Residue construction commencing." << std::endl;
    // debug("Residue construction commencing.\n");
    id_global = -1;
    id_local = -1;
    resid = -1;
    num_atoms_res = -1;
}
inline Residue::Residue(const Residue &obj)
{
    // Copy Constructor.
    std::cout << "Allocating with copy constructor.(R)" << std::endl;
    // ptr = new int;
    // *ptr = *obj.ptr; // copy the value
    id_global = obj.id_global;
    id_local = obj.id_local;
    resid = obj.resid;
    num_atoms_res = obj.num_atoms_res;
    restype = obj.restype;
}
inline Residue::~Residue()
{
    // std::cout << "Residue destruction commencing." << std::endl;
    // debug("Residue destruction commencing.\n");

}
inline void Residue::print_Prop()
{
    std::cout << "Residue: "
              << "\n\tglobal,local,resid: "
              << id_global << " " << id_local << " " << resid << std::endl;

}


/* ---------------------------------------------------------
   Atom
   --------------------------------------------------------- */
class Atom: public Residue
{
public:
    Atom();
    Atom(const Atom &obj);
    // Atom(const Atom &a) : Atom(a);


    ~Atom();
    void print_Coords();
    void print_Atom();

    float x, y, z;
    int index;

    float radius;     // VDW, Richards, Gerstein.
    std::string atomtype;  // N, CA, C, O, CB, ..
    std::string general_atomtype; // 77,1

    std::string chain; // "A"
    int chainid;
    int resid;

private:
    /* std::string str_index = line.substr(5,6);      // int */
    /* std::string str_atomtype = line.substr(13,3); */
    /* std::string str_restype = line.substr(17,3); */
    /* std::string str_chain = line.substr(21,1); */
    /* std::string str_resid = line.substr(23,3);     // int */
    /* std::string str_x = line.substr(30,8);         // float */
    /* std::string str_y = line.substr(38,8);         // float */
    /* std::string str_z = line.substr(46,8);         // float */
    int *ptr;
};
inline Atom::Atom()
{
    // std::cout << "Atom construction commencing." << std::endl;
    // debug("Atom construction commencing.\n");
    x = 0.00000;
    y = 0.00000;
    z = 0.00000;
    chain = "zz";
    chainid = -1;
    resid = -1;
    radius = 0.0000;
    index = -1;

}
inline Atom::Atom(const Atom &obj)
{
    // Copy Constructor.
    // std::cout << "Allocating with copy constructor." << std::endl;


    // ptr = new int;
    // *ptr = *obj.ptr; // copy the value
    x = obj.x;
    y = obj.y;
    z = obj.z;
    radius = obj.radius;
    index = obj.index;

    atomtype = obj.atomtype;
    general_atomtype = obj.general_atomtype;

    chain = obj.chain;
    chainid = obj.chainid;
    resid = obj.resid;
}
inline Atom::~Atom()
{
    // std::cout << "Atom destruction commencing." << std::endl;
    // debug("Atom destruction commencing.\n");
    // delete ptr;
    // delete obj;

}
inline void Atom::print_Atom()
{
    std::cout << "index: "
              << index << "   "
              << "x: "
              << x << "  "
              << "y: "
              << y << "  "
              << "z: "
              << z << "  "
              << std::endl;
    std::cout << "chain: " << chain
              << " resid: "  << resid
              << std::endl;
}
inline void Atom::print_Coords()
{
    printf("ATOM %6d  %3s %3s          %8.3f%8.3f%8.3f",
           index,atomtype.c_str(),
           restype.c_str(),x,y,z);
    printf("                       %s\n",general_atomtype.c_str());
}

// inline void Atom::print_Coords()
// {
//     printf("ATOM %6d  %3s %3s          %8.3f%8.3f%8.3f",
//            index,atomtype.c_str(),
//            restype.c_str(),x,y,z);
//     printf("                       %s\n",general_atomtype.c_str());

// }
// inline void Atom::print_info()
// {
//     // std::cout << restype << ": " << atomtype << " " << general_atomtype << std::endl;
//     printf("%s: %s (%s)\n",restype.c_str(),atomtype.c_str(),general_atomtype.c_str());
// }


/* ---------------------------------------------------------
   function declarations
   --------------------------------------------------------- */
// int system_select_atomtype(Atom *aa,char const *criterion,int total);
int system_select_atomtype(Atom *aa,char const *criterion,int total,Atom *asel);
int system_select(Atom *aa,char const *criterion,int total);
void system_select(Atom *aa,char const *criterion,int total,Atom *asel);


// vAtoms, vvAtoms | Vector of Atoms, Vector of Vector of Atoms
typedef std::vector<Atom> vAtoms;
typedef std::vector<std::vector<Atom>> vvAtoms;

// vpAtoms,
typedef std::vector<Atom *> vpAtoms;
typedef std::vector<std::vector<Atom *>> vvpAtoms;

typedef std::vector<int> IndexGroup;
typedef std::vector<std::vector<int>> vIndexGroup;

// SELECT:
// vpAtoms select(vAtoms aa,char const *criterion);
IndexGroup select(vAtoms aa,char const *criterion);

vAtoms set_chainid(vAtoms aa);
// vvpAtoms sort_segment_chain(vAtoms aa);
Vector get_centroid(IndexGroup ig,vAtoms aa);


void get_minmax(System sys);
vIndexGroup sort_segment_chain(vAtoms aa);

// typedef std::vector<
// typedef std::vector<Atom> Dihedral(4);
// Atoms aa_backbone;

#endif
