# // 1: c++, c: contacts
# // 2: name:
# // .contacts
#ifndef _contacts_hpp_
#define _contacts_hpp_

/* ---------------------------------------------------------
   libraries:
   --------------------------------------------------------- */
// #include <stdio.h>
// #include <stdlib.h> // strtod?, stod
// #include <assert.h>
// #include <string.h> // string
// #include <cctype.h>
/* #include <algorithm> // remove_if, count */
/* #include <iostream> */
/* #include <fstream> */
/* #include <ctime> */
// #include <list>        // std::list
// #include <vector>
// #include <iterator> // istream_iterator
#include "system.hpp"
#include "dcd.h"
// #include "dcdio.h" // inside of dcd.h
#include "boost/tuple/tuple.hpp"
#include "microtubule.hpp"

/* ---------------------------------------------------------
   headers:
   --------------------------------------------------------- */
#include "debug.h"
/* #include "system.hpp" */
// #include "atom.hpp"



/* ---------------------------------------------------------
   Definitions:
   --------------------------------------------------------- */
/* #define BUFFERSIZE 900 */



/* ---------------------------------------------------------
   Classes:
   --------------------------------------------------------- */
// header_class
typedef boost::tuple<int,int,double> Contact;
typedef std::vector<Contact> SetContacts;
typedef std::vector<SetContacts> SetNeighbors;
typedef std::vector<SetNeighbors> SetChains;
typedef std::vector<SetChains> SetGlobalContacts;


/* ---------------------------------------------------------
   Contacts
   --------------------------------------------------------- */
// class Contacts
// {
// public:
//     Contacts();
//     Contacts(const Contacts &obj);
//     ~Contacts();
//     void print_all_contacts();

//     // int size;
//     int index;

//     std::vector<std::vector<boost::tuple<int,int,double>>> frame_contacts;
//     std::vector<boost::tuple<int,int,double>> chain_contacts;
//     std::vector<boost::tuple<int,int,double>> initial_contacts;

// private:

// };
// inline Contacts::Contacts()
// {
//     // std::cout << "Contacts construction commencing." << std::endl;
//     // debug("Contacts construction commencing.\n");
//     // chainid = -1;
//     // num_residues = -1;
//     // chain = "zz";
// }
// inline Contacts::Contacts(const Contacts &obj)
// {
//     // Copy constructor.
//     index = obj.index;
//     // chain = obj.chain;
//     // num_residues = obj.num_residues;
//     // chainid = obj.chainid;
// }
// inline Contacts::~Contacts()
// {
//     // std::cout << "Contacts destruction commencing." << std::endl;
//     // debug("Contacts destruction commencing.\n");
//     // chain_contacts.clear();
//     frame_contacts.clear();

// }
// inline void Contacts::print_all_contacts()
// {
//     std::cout << "Contacts: " << std::endl;

//     // for(auto c: contact)
//     // {
//     //     std::cout << c[0] << " "
//     //               << c[1] << " "
//     //               << c[2] << std::endl;
//     // }
// }

// class Residue: public Chain
// {
// public:
//     Residue();
//     Residue(const Residue &obj);
//     ~Residue();
//     void print_prop();

//     int id_global;
//     int id_local;
//     int resid;
//     int num_atoms_res;
//     std::string restype;   // HIS, GLU, ILE, LEU ..

// private:

// };
// inline Residue::Residue()
// {
//     // std::cout << "Residue construction commencing." << std::endl;
//     // debug("Residue construction commencing.\n");
//     id_global = -1;
//     id_local = -1;
//     resid = -1;
//     num_atoms_res = -1;
// }




/* ---------------------------------------------------------
   function declarations:
   --------------------------------------------------------- */
// void ReadPDBfile(PDBfile *pdbfile,char filename[40]);
// void get_contacts(Atom *a1,char *argv,int num_atoms);
void get_contacts(Atom *a1,Atom *a2,char dcdfilename[40],int num_atoms);
// void get_map_of_mtneighbors(std::vector<Atom*> chain_ref,std::vector<std::vector<int>> matrix);
// void get_map_of_mtneighbors(std::vector<std::vector <Atom>> chain_ref,std::vector<std::vector<int>> matrix,
//                             std::vector<std::pair<int,int>> dimers);
std::vector<std::vector<int>> get_map_of_mtneighbors(std::vector<std::vector <Atom>> chain_ref,
                                                     std::vector<std::pair<int,int>> dimers);

// void get_contacts_for_chain();
// std::vector<boost::tuple<int,int,int>> get_contacts_for_chain(std::vector<std::vector <Atom>> chain_ref);
// std::vector<boost::tuple<int,int,int,double>> get_contacts_for_chain(std::vector <Atom> chain1,
//                                                                      std::vector <Atom> chain2,
//                                                                      float cutoff);
// std::vector<boost::tuple<int,int,double>> get_contacts_for_chain(std::vector <Atom> chain1,
//                                                                  std::vector <Atom> chain2,
//                                                                  float cutoff);
// std::vector<boost::tuple<int,int,double>> get_contacts_for_chain(std::vector <Atom> chain1,
//                                                                  std::vector <Atom> chain2,
//                                                                  float cutoff,
//                                   std::vector<boost::tuple<int,int,double>> chain_contacts);
// void get_contacts_for_chain(std::vector <Atom> chain1,
//                             std::vector <Atom> chain2,
//                             float cutoff,
//                             std::vector<boost::tuple<int,int,double>> contacts);
SetContacts get_contacts_for_chain(std::vector <Atom> chain1,
                                   std::vector <Atom> chain2,
                                   double cutoff);
SetContacts get_contacts_for_chain(std::vector <Atom> chain1,
                                   double cutoff);



// std::vector<boost::tuple<int,int,int,double>> get_contacts_for_chain_later(Atom *aalater,
//                                                                            double cutoff,
//                                                                            double tolerance,
//                                                                            std::vector<boost::tuple
//                                                                            <int,int,int,double>> contacts);


// SetContacts get_contacts_for_chain_later(Atom *alater,
//                                          double cutoff,
//                                          double tolerance,
//                                          SetContacts contacts);
SetContacts get_contacts_for_chain_later(Atom *alater,
                                         double cutoff,
                                         double tolerance,
                                         SetContacts contacts);





// std::vector<boost::tuple<int,int,int,double>> output_contacts(std::vector<boost::tuple
//                                                               <int,int,int,double>> contacts);
std::vector<boost::tuple<int,int,int,double>> output_contacts(std::vector<std::vector<boost::tuple
                                                              <int,int,int,double>>> contacts);
void output_global_contacts(SetGlobalContacts gc);
// void explore_global_contacts(SetGlobalContacts gc);
SetGlobalContacts explore_global_contacts(SetGlobalContacts gc,MtIndexMap map);


// chain_contact = get_contacts_for_chain_later(aa_later,8.0,2.0,chain_contacts_0[0]);


// std::vector<boost::tuple<int,int,int,double>> get_contacts_for_chain(std::vector <Atom> chain1,
//                                                                      std::vector <Atom> chain2,
//                                                                      float cutoff);
    // ,
                                                              // float extra);


#endif
