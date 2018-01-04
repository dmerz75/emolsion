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
#include <math.h>
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
    std::vector<int> chainids;
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
    void get_unique_chainids(vAtoms aa, vIndexGroup isel_chains);
    void print_Chainids();
    void get_first_ab_axes(vAtoms aa);
    void identify_chains_on_pf(vAtoms aa,vIndexGroup isel_chain);
    void get_bending_angle(vAtoms aa,vIndexGroup isel_chain);

private:

};
inline SystemPF::SystemPF()
{
    // constructor:
    num_protofilaments = 13; // use default.. ?
}
inline void SystemPF::print_Chainids()
{
    for(auto c: chainids)
    {
        std::cout << c
                  << " ";
    }
    std::cout << std::endl;

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

        for(int i=0; i < pf.size(); i++)
        {

            std::cout << pf[i].first << ", " << pf[i].second << ", ";
        }
        std::cout << std::endl;
    }


}
inline void SystemPF::get_unique_chainids(vAtoms aa, vIndexGroup chains)
{


    for(int i=0; i < chains.size(); i++)
    {
        // std::cout << aa[chains[i][0]].chainid << std::endl;
        chainids.push_back(aa[chains[i][0]].chainid);
    }
    // exit(0);
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
        // std::cout << pf[0].first << ", " << pf[0].second << std::endl;

        select_a.str(""); // clear string.
        select_b.str(""); // clear string.

        select_a << "chainid " << pf[0].first; // build string: chained 0
        select_b << "chainid " << pf[0].second;
        // std::cout << "alpha:" << sel_a.str();
        // std::cout << "  beta:" << sel_b.str() << std::endl;
        indexgroup_a = select(aa,select_a.str().c_str());
        indexgroup_b = select(aa,select_b.str().c_str());
        // std::cout << "Selection_a: " << indexgroup_a.size() << std::endl;
        // std::cout << "Selection_b: " << indexgroup_b.size() << std::endl;

        centroid_a = get_centroid(indexgroup_a,aa);
        centroid_b = get_centroid(indexgroup_b,aa);
        // centroid_a.print_Vector();
        // centroid_b.print_Vector();

        axis = get_vector(centroid_a,centroid_b);
        // axis.print_Vector();

        vAxis.push_back(boost::make_tuple(axis,centroid_a,centroid_b));
    }
}
inline void SystemPF::identify_chains_on_pf(vAtoms aa,vIndexGroup isel_chain)
{
    std::cout << "identify chains on the same protofilament." << std::endl;

    std::vector<int> unused_chainids;
    unused_chainids = chainids;

    Vector centroid;
    Vector centroida, centroidb;

    int i, a, b, csel;
    i = 0;
    a = b = csel = -1;

    int counter = 0;
    std::ostringstream chainid;

    Vector vector_ab, vector_ac;
    Vector nvec_ab, nvec_ac;

    double dcos, dsin, magnitude_ac;
    dcos = dsin = 0.0;

    Dimer dimer;


    // remove used chainids:
    // for(int c=0; c < unused_chainids.size; c++)
    for(auto pf: protofilaments)
    {
        // pf[0].first
        // pf[0].second
        unused_chainids.erase(std::find(unused_chainids.begin(),
                                        unused_chainids.end(),
                                        pf[0].first));
        unused_chainids.erase(std::find(unused_chainids.begin(),
                                        unused_chainids.end(),
                                        pf[0].second));
    }
    std::cout << "Remaining unused chainids: " << unused_chainids.size() << std::endl;


    for(int np=0; np < num_protofilaments - 1; np++)
    {
        std::cout << "dimer iteration: " << np << std::endl;

        i = 0; // reset the counting through protofilaments.
        for(auto pf: protofilaments) // of pair/tuple. (0,1),(2,3)..
        {
            // std::cout << pf[0].first << ", " << pf[0].second << std::endl;

            // USE THIS
            // std::cout << i << " --which protofilament-- "
            //           << pf.back().first
            //           << ", "
            //           << pf.back().second
            //           << std::endl;

            vector_ab = vAxis[i].get<0>();
            centroida = vAxis[i].get<1>();
            centroidb = vAxis[i].get<2>();

            // for(int c=pf.back().first; c < isel_chain.size(); c++)
            // for(auto c: isel_chain)
            // for(int c=0; c < unused_chainids.size(); c++)
            while(!unused_chainids.empty())
            {
                // std::cout << "chainid: " << aa[c[0]].chainid
                //           << " (" << c.size() << ")" << std::endl;


                // csel = unused_chainids[0];
                csel = unused_chainids[0];
                // std::cout << "csel(1): " << csel << std::endl;
                // std::cout << "i: " << i << std::endl;
                // std::cout << "ab: " << a << ", " << b << std::endl;



                // if(i > isel_chain.size())
                // {
                //     break;
                // }
                if(counter + num_protofilaments * 2 + 1 >= isel_chain.size())
                {
                    break;
                }


                centroid = get_centroid(isel_chain[csel],aa);
                // // centroid.print_Vector();
                vector_ac = get_vector(centroida,centroid);
                magnitude_ac = magnitude(vector_ac);

                // if((magnitude_ac < 70.0) || (magnitude_ac > 100.0))
                // if(magnitude_ac > 110.0)
                // {
                //     continue;
                // }

                // normalize, cosine, sine of a-b, a-c
                nvec_ab = normalize(vector_ab);
                nvec_ac = normalize(vector_ac);
                dcos = get_costheta(nvec_ac,nvec_ab);
                dsin = get_sintheta(nvec_ac,nvec_ab);

                // nvec_ab.print_Vector();
                // nvec_ac.print_Vector();
                // std::cout << "magnitude: " << magnitude_ac << std::endl;
                // std::cout << "costheta: " << dcos << "   " << dsin << " :sintheta." << std::endl;
                // std::cout << "counter: " << counter << std::endl;

                // Not.
                // Remove a,b:
                if((dcos > 0.970) && (dsin < 0.06))
                {
                    counter += 1;

                    if(isel_chain[csel].size() > 431)
                    {
                        a = csel;
                        unused_chainids.erase(std::find(unused_chainids.begin(),
                                                        unused_chainids.end(),
                                                        a));
                    }
                    else
                    {
                        b = csel;
                        unused_chainids.erase(std::find(unused_chainids.begin(),
                                                        unused_chainids.end(),
                                                        b));
                    }
                }
                else
                {
                    i += 1;
                    continue;
                }


                if((a == -1) || (b == -1))
                {
                    continue;
                }

                // std::cout << "add_ab: " << a << ", " << b << std::endl;
                // std::cout << "unused_size: " << unused_chainids.size() << std::endl;
                dimer = std::make_pair(a,b);
                protofilaments[i].push_back(dimer);
                // counter = 0;
                a = b = -1;
                break;
            }

            i += 1; // 0-12 for the 13 protofilaments.
            // break;
        } // the protofilaments.

    } // np for loop
}
inline void SystemPF::get_bending_angle(vAtoms aa,vIndexGroup isel_chain)
{
    std::cout << "Getting Bending Angle." << std::endl;

    // FILE
    FILE * fp_bending_angle;
    fp_bending_angle = fopen("emol_mtpf_bending_angle.dat", "a+");
    fprintf(fp_bending_angle,"#\n");
    // fprintf(fp_bending_angle,"\n");



    // std::ostringstream select_a; // str: "chainid 0"
    // std::ostringstream select_b;
    // IndexGroup indexgroup_a;   // (1st) alpha,beta
    // IndexGroup indexgroup_b;   // indices, and selection.
    // Vector centroid_a;
    // Vector centroid_b;
    // Vector axis;

    int num_angles; // angles (dimers * 0.5 + 1), for 13: should be 7
    int start;
    int pf1, pf2, pf3, pf4;
    Vector cen1, cen2, cen3, cen4;
    Vector v12, v34, n12, n34;
    double rad_ang, deg_ang;

    // 0-1;
    // 1-2;
    // 2-3; *
    // 3-4; *
    // 4-5; *
    // 5-6; *
    // 6-7; *
    // 7-8; *
    // 8-9; *
    // 9-10;
    // 10-11;

    for(auto pf: protofilaments)
    {
        // std::cout << pf[0].first << ", " << pf[0].second << std::endl;
        // std::cout << pf.size() << std::endl;
        num_angles = pf.size() * 0.5 + 1;
        start = (pf.size() - num_angles) * 0.5; // for the 13 case, starts at 2-3.

        for(int s=start; s < start + num_angles; s++)
        {
            pf1 = s;
            pf2 = s + 1;
            pf3 = s + 2;
            pf4 = s + 3;

            cen1 = get_centroid(isel_chain[pf1],aa);
            cen2 = get_centroid(isel_chain[pf2],aa);
            cen3 = get_centroid(isel_chain[pf3],aa);
            cen4 = get_centroid(isel_chain[pf4],aa);

            v12 = get_vector(cen1,cen2);
            v34 = get_vector(cen3,cen4);

            n12 = normalize(v12);
            n34 = normalize(v34);

            rad_ang = get_costheta(n12,n34);
            deg_ang = acos(rad_ang);

            fprintf(fp_bending_angle,"%4.1f ",deg_ang);
        }

        fprintf(fp_bending_angle,"\n");

        // Dimers:
        // for(auto d: pf)
        // {
        //     std::cout << "dimer: " << d.first
        //               << " " << d.second
        //               << std::endl;
        // }
    }

    fclose(fp_bending_angle);

}

/* ---------------------------------------------------------
   function declarations:
   --------------------------------------------------------- */
// void ReadPDBfile(PDBfile *pdbfile,char filename[40]);
// Protofilaments determine_num_protofilaments(vAtoms aa);
// Protofilaments get_full_protofilament(vAtoms aa);

#endif
