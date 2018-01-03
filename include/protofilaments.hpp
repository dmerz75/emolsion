# // 1: c++, c: protofilaments
# // 2: name:
# // .protofilaments
#ifndef _protofilaments_
#define _protofilaments_

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
/* #include <vector> */
// #include <iterator> // istream_iterator



/* ---------------------------------------------------------
   headers:
   --------------------------------------------------------- */
#include "debug.h"
#include "system.hpp"
#include <stdio.h>
#include <sstream>
#include "boost/tuple/tuple.hpp"
// #include "atom.hpp"
// #include "microtubule.hpp"



/* ---------------------------------------------------------
   Definitions:
   --------------------------------------------------------- */
/* #define BUFFERSIZE 900 */

typedef std::pair<int,int> Dimer;
typedef std::vector<Dimer> Protofilament;
typedef std::vector<Protofilament> Protofilaments;


/* ---------------------------------------------------------
   Classes:
   --------------------------------------------------------- */
// header_class

class SystemPF {

public:

    // Constructor
    SystemPF();

    // Destructor
    ~SystemPF(){
    };

    int num_protofilaments; // defaults 13.
    Protofilaments protofilaments; // initially, 0-1, .. 26-27,

    std::vector<boost::tuple<Vector,Vector,Vector>> vAxis;
    // first alpha-beta axis/vector;
    // vAxis, centroids_1A, centroids_1B;
    // std::vector<Vector> centroids_1A;
    // std::vector<Vector> centroids_1B;

    IndexGroup vindexgroup_A;   // all alpha,beta indices.
    IndexGroup vindexgroup_B;


    void print_PF();
    void build_initial_dimers();
    void get_first_ab_axes(vAtoms aa);
    void identify_chains_on_pf(vAtoms aa,vIndexGroup isel_chain);


private:

};

inline SystemPF::SystemPF()
{
    num_protofilaments = 13; // use default.. ?
}
inline void SystemPF::print_PF()
{
    std::cout << "This SystemPF has " << num_protofilaments
              << " protofilament(s)."
              << std::endl;

    std::cout << "They are: " << std::endl;

    // 1 dimer deep in the protofilament stack.
    for(auto pf:protofilaments)
    {
        std::cout << pf[0].first << ", " << pf[0].second << std::endl;
    }

}

inline void SystemPF::build_initial_dimers()
{

    Dimer dimer;
    Protofilament pf;
    // Protofilaments lst_allpf;

    int a,b;
    // typedef std::pair<int,int> Dimer;
    // typedef std::vector<Dimer> Protofilament;
    // typedef std::vector<Protofilament> Protofilaments;

    std::cout << "Building initial dimers." << std::endl;


    // typically 13 protofilaments
    for (int i=0; i < num_protofilaments; i++)
    {
        a = 2 * i;
        b = 2 * i + 1;
        // std::cout << "protofilament: " << i << std::endl;
        // std::cout << a << b << std::endl;
        dimer = std::make_pair(a,b);
        pf.push_back(dimer);
        protofilaments.push_back(pf);
        pf.clear();
        // dimer.clear();
    }
    // return lst_allpf;
}
inline void SystemPF::get_first_ab_axes(vAtoms aa)
{
    std::cout << "Getting first ab axes." << std::endl;

    std::ostringstream select_a; // str: "chainid 0"
    std::ostringstream select_b;
    IndexGroup indexgroup_a;   // (1st) alpha,beta
    IndexGroup indexgroup_b;   // indices, and selection.
    Vector centroid_a;
    Vector centroid_b;

    Vector axis;

    for(auto pf: protofilaments)
    {
        std::cout << pf[0].first << ", " << pf[0].second << std::endl;

        select_a.str(""); // clear string.
        select_b.str(""); // clear string.

        select_a << "chainid " << pf[0].first; // build string: chained 0
        select_b << "chainid " << pf[0].second;
        // std::cout << "alpha:" << sel_a.str();
        // std::cout << "  beta:" << sel_b.str() << std::endl;
        indexgroup_a = select(aa,select_a.str().c_str());
        indexgroup_b = select(aa,select_b.str().c_str());
        std::cout << "Selection_a: " << indexgroup_a.size() << std::endl;
        std::cout << "Selection_b: " << indexgroup_b.size() << std::endl;

        centroid_a = get_centroid(indexgroup_a,aa);
        centroid_b = get_centroid(indexgroup_b,aa);
        // centroid_a.print_Vector();
        // centroid_b.print_Vector();

        axis = get_vector(centroid_b,centroid_a);
        axis.print_Vector();

        vAxis.push_back(boost::make_tuple(axis,centroid_a,centroid_b));
    }
}
inline void SystemPF::identify_chains_on_pf(vAtoms aa,vIndexGroup isel_chain)
{
    std::cout << "identify chains on the same protofilament." << std::endl;

    Vector centroid;
    Vector centroida, centroidb;
    Vector vector_ab;
    Vector vector_ac;
    int i = 0;
    int counter = 0;
    // int index = -1;
    std::ostringstream chainid;

    Vector vector_anew;

    Vector nvec_ab, nvec_anew;
    double dcos, dsin, magnitude_ac;
    dcos = dsin = 0.0;

    for(auto c: isel_chain)
    {

        std::cout << c.size() << std::endl;
        // index = c[0];
        std::cout << "chainid: " << aa[c[0]].chainid << std::endl;
        // std::cout << "chainid: " << c[0] << std::endl;
        centroid = get_centroid(c,aa);
        centroid.print_Vector();

        vector_ab = vAxis[i].get<0>();
        centroida = vAxis[i].get<1>();
        centroidb = vAxis[i].get<2>();

        vector_ac = get_vector(centroida,centroid);
        magnitude_ac = magnitude(vector_ac);
        if(magnitude_ac < 70.0)
        {
            continue;
        }


        vector_anew = get_vector(centroida,centroid);
        nvec_ab = normalize(vector_ab);
        nvec_anew = normalize(vector_anew);


        dcos = get_costheta(nvec_ab,nvec_anew);
        if(dcos < 0.9)
        {
            continue;
        }
        counter += 1;
        dsin = get_sintheta(nvec_ab,nvec_anew);

        nvec_ab.print_Vector();
        nvec_anew.print_Vector();

        std::cout << "costheta: " << dcos << "   >>  sintheta: " << dsin << std::endl;
        std::cout << "counter: " << counter;
        std::cout << " chainid: " << aa[c[0]].chainid << std::endl;


        i += 1;
    }

}

/* ---------------------------------------------------------
   function declarations:
   --------------------------------------------------------- */
// void ReadPDBfile(PDBfile *pdbfile,char filename[40]);
// Protofilaments determine_num_protofilaments(vAtoms aa);
// Protofilaments get_full_protofilament(vAtoms aa);

#endif
