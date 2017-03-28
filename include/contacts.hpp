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
// #include "dcdio.h" // inside of dcd.h
#include "system.hpp"
#include "dcd.h"
#include "microtubule.hpp"
#include "boost/tuple/tuple.hpp"

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


/* ---------------------------------------------------------
   Typdefs:
   --------------------------------------------------------- */
// for the contact arrarys, in time, by chain, by neighbor.
// Contacts, Set of Contacts, 9 Neighbors, Set of Chains/dimers ~ 156,
// Lastly, by frame
typedef boost::tuple<int,int,double,int,int> Contact;
typedef std::vector<Contact> SetContacts;
typedef std::vector<SetContacts> SetNeighbors;
typedef std::vector<SetNeighbors> SetChains;
typedef std::vector<SetChains> SetGlobalContacts;


/* ---------------------------------------------------------
   function declarations:
   --------------------------------------------------------- */
void get_contacts(Atom *a1,Atom *a2,char dcdfilename[40],int num_atoms);

// MtNeighbors get_map_of_mtneighbors(std::vector<std::vector <Atom>> chain_ref,
//                                    DimerList dimers);
// MtNeighbors get_map_of_mtneighbors(vvpAtoms chains,
//                                    DimerList dimers);
MtNeighbors get_map_of_mtneighbors(vIndexGroup isel_chain,vAtoms aa,DimerList dimers);
void print_mt_map(MtNeighbors mt_matrix);

// SetContacts get_contacts_for_chain(std::vector <Atom> chain1,
//                                    std::vector <Atom> chain2,
//                                    double cutoff);
// SetContacts get_contacts_for_chain(std::vector <Atom> chain1,
//                                    double cutoff);
// SetContacts get_contacts_for_chain(std::vector <Atom> chain1,
//                                    double cutoff,
//                                    MtIndexMap map,
//                                    int cid);
// SetContacts get_contacts_for_chain(std::vector <Atom> chain1,
//                                    std::vector <Atom> chain2,
//                                    double cutoff,
//                                    MtIndexMap map,
//                                    int cid1,
//                                    int cid2);


// SM
// vAtoms chain2,
// SM - inter
SetContacts get_contacts_for_chain(vAtoms aa,
                                   double cutoff,
                                   IndexGroup ig1,
                                   IndexGroup ig2);

// MT - inter
SetContacts get_contacts_for_chain(vAtoms aa,
                                   double cutoff,
                                   MtIndexMap map,
                                   IndexGroup ig1,
                                   IndexGroup ig2,
                                   int cid1,int cid2);
// contact_set = get_contacts_for_chain(allatoms_ref,8.0,isel_chain[c[0]],
//                                      c[0]);

// MT - intra
SetContacts get_contacts_for_chain(vAtoms aa,
                                   double cutoff,
                                   MtIndexMap map,
                                   IndexGroup ig1,
                                   int cid1);

// SetContacts get_contacts_for_chain(vAtoms aa,
//                                    double cutoff,
//                                    IndexGroup ig1);
// // MT
// SetContacts get_contacts_for_chain(vAtoms chain1,
//                                    double cutoff,
//                                    MtIndexMap map,
//                                    int cid);
// SetContacts get_contacts_for_chain(vAtoms chain1,
//                                    double cutoff,
//                                    MtIndexMap map,
//                                    int cid1,
//                                    int cid2);

// vAtoms chain2,


// SetContacts get_contacts_for_chain_later(Atom *alater,
//                                          double cutoff,
//                                          double tolerance,
//                                          SetContacts contacts);
SetContacts get_contacts_for_chain_later(vAtoms allatoms,
                                         double cutoff,
                                         double tolerance,
                                         SetContacts contacts);

void print_global_contacts(SetGlobalContacts gc);
void print_global_contacts_count(SetGlobalContacts gc);
void output_global_contacts(SetGlobalContacts gc);
void output_global_contacts_by_subdomain(SetGlobalContacts gc);

#endif
