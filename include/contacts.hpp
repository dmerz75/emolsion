# // 1: c++, c: contacts
# // 2: name:
# // .contacts
#ifndef _contacts_hpp_
#define _contacts_hpp_

/* ---------------------------------------------------------
   libraries & headers
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
#include "debug.h"


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
//
// i1 index 1, i2 index 2
// dist: distance between atoms
// mt1: sector of a tubulin monomer for index 1 (contact)
// mt2: sector of a tubulin monomer for index 2 (contacted)
// eh: energy scaling of the Lennard-Jones potential
//
//                    i1  i2  dist  mt1 mt2 eh
typedef boost::tuple<int,int,double,int,int,double> Contact;
typedef std::vector<Contact> SetContacts;
typedef std::vector<SetContacts> SetNeighbors;
typedef std::vector<SetContacts> SetFrameContacts;
typedef std::vector<SetNeighbors> SetChains;
typedef std::vector<SetChains> SetGlobalContacts;


/* ---------------------------------------------------------
   function declarations:
   --------------------------------------------------------- */
void get_contacts(Atom *a1,Atom *a2,char dcdfilename[40],int num_atoms);

// SM
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
// MT - intra
SetContacts get_contacts_for_chain(vAtoms aa,
                                   double cutoff,
                                   MtIndexMap map,
                                   IndexGroup ig1,
                                   int cid1);


// SM or MT (time later)
SetContacts get_contacts_for_chain_later(vAtoms allatoms,
                                         double cutoff,
                                         double tolerance,
                                         SetContacts contacts);

MtNeighbors get_map_of_mtneighbors(vIndexGroup isel_chain,vAtoms aa,DimerList dimers);
void print_mt_map(MtNeighbors mt_matrix);

// code development
void print_set_contacts(SetContacts cn);
void print_global_contacts(SetGlobalContacts gc);
void print_global_contacts_count(SetGlobalContacts gc);

SetContacts set_eh_contacts(SetContacts cn,double eh);

// for data analysis
void output_global_contacts(SetGlobalContacts gc);
void output_global_contacts_by_subdomain(SetGlobalContacts gc);
void output_framecontact_set(SetFrameContacts framecontact_set);

// for topology creation
void write_contacts_to_file(FILE *fp_topology,SetContacts contact_set);
SetContacts read_contacts_from_file(char filename[40]);
SetContacts sort_contacts(SetContacts cn,MtIndexMap mtmap,int c1, int c2);

#endif
