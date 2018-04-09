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
    void get_bending_angle1(vAtoms aa,vIndexGroup isel_chain);
    void get_bending_angle2(vAtoms aa,vIndexGroup isel_chain);


    void get_distance_centroid(vAtoms aa,vIndexGroup isel_chain);
    void get_beta_angle(vAtoms aa,vIndexGroup isel_chain);

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

    int c = 0;

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
            std::cout << i << " --which protofilament-- "
                      << pf.back().first
                      << ", "
                      << pf.back().second
                      << std::endl;

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

                // csel: integer of the chainid.
                csel = unused_chainids[0];
                // std::cout << "csel(pf): " << csel << std::endl;
                // std::cout << "i: " << i << std::endl;
                // std::cout << "ab: " << a << ", " << b << std::endl;

                // if(i > isel_chain.size())
                // {
                //     break;
                // }
                if(i >= num_protofilaments)
                {
                    break;
                }

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
                // Hard cosine and sine.
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

            // std::cout << "num_pf: " << num_protofilaments << std::endl;
            // if(i >= num_protofilaments)
            // {
            //     break;
            // }

        } // the protofilaments.

    } // np for loop
}
inline void SystemPF::get_bending_angle1(vAtoms aa,vIndexGroup isel_chain)
{
    // Description:
    // The bending angle is described as the 2-3 monomer-monomer vector
    // dotted with the 4-5 monomer-monomer vector. Take the acos.
    // Angles should increase from 0 degrees to 30++.
    // std::cout << "Getting Bending Angle." << std::endl;

    // FILE
    FILE * fp_bending_angle;
    fp_bending_angle = fopen("emol_mtpfbending_angle.dat", "a+");
    fprintf(fp_bending_angle,"#\n");

    // Variables:
    int num_angles; // angles (dimers * 0.5 + 1), for 13: should be 7
    int start;
    int pf1, pf2, pf3, pf4;

    Vector cen1, cen2, cen3, cen4;
    Vector v12, v23, v34;
    Vector n12, n23, n34;
    Vector cross13, cross14, cross24;
    Vector n13, n24;

    double cos_ang, acos_ang, sign;
    double m23;

    // Vector cendist;
    // double mag_cendist;
    // double sin_ang, asin_ang;
    // std::cout << "Num_Protofilaments: " << num_protofilaments << std::endl;
    // std::cout << "Num_Dimers: " << protofilaments[0].size() << std::endl;

    if(protofilaments[0].size() == 12)
    {

    }
    else if(protofilaments[0].size() == 8)
    {

    }
    else
    {
        fprintf(stderr,"[WARN] %s:%d: errno: %s\n",__FILE__,__LINE__,
                "Protofilaments' bending angles were not computed.");
        return;
    }
    num_angles = protofilaments[0].size() * 0.5 + 1;
    start = (protofilaments[0].size() - num_angles) * 0.5; // 12->2, 8->1 (pos)
    // 12:  7  (position 0 .. 1 .. 2)  --> 4-5
    //  8:  5  (position 0 .. 1)       --> 2-3
    // Dimers:   12, 8
    // 1.  0-1
    // 2.  2-3       _ 2-3  v 4-5
    // 3.  4-5    _  _ 4-5  v 6-7
    // 4.  6-7    _  _ 6-7  v 8-9
    // 5.  8-9    _  _ 8-9  v 10-11
    // 6.  10-11  _  _ 10-11 v 12-13
    // 7.  12-13  _
    // 8.  14-15  _
    // 9.  16-17  _
    // 10. 18-19
    // 11. 20-21
    // 12. 22-23

    for(auto pf: protofilaments)
    {
        // std::cout << pf[0].first << ", " << pf[0].second << std::endl;
        // std::cout << pf.size() << std::endl;
        // num_angles = pf.size() * 0.5 + 1;
        // start = (pf.size() - num_angles) * 0.5; // for the 13 case, starts at 2-3.
        for(int s=start; s < start + num_angles; s++)
        {
            // Identify 4 monomers.
            pf1 = pf[s].first;
            pf2 = pf[s].second;
            pf3 = pf[s + 1].first;
            pf4 = pf[s + 1].second;
            // std::cout << "Monomers: " << pf1
            //           << " " << pf2
            //           << " " << pf3
            //           << " " << pf4
            //           << std::endl;
            // 20, 21
            //     Monomers: 72 73 98 99
            //     Monomers: 98 99 124 125
            //     Monomers: 124 125 150 151
            //     Monomers: 150 151 176 177
            //     Monomers: 176 177 202 203
            //     Monomers: 202 203 228 229
            //     Monomers: 228 229 254 255
            //     22, 23
            //     Monomers: 74 75 100 101
            //     Monomers: 100 101 126 127
            //     Monomers: 126 127 152 153
            //     Monomers: 152 153 178 179
            //     Monomers: 178 179 204 205
            //     Monomers: 204 205 230 231
            //     Monomers: 230 231 256 257
            //     24, 25
            //     Monomers: 76 77 102 103
            //     Monomers: 102 103 128 129
            //     Monomers: 128 129 154 155
            //     Monomers: 154 155 180 181
            //     Monomers: 180 181 206 207
            //     Monomers: 206 207 232 233
            //     Monomers: 232 233 258 259

            // 4 Centroids.
            cen1 = get_centroid(isel_chain[pf1],aa);
            cen2 = get_centroid(isel_chain[pf2],aa);
            cen3 = get_centroid(isel_chain[pf3],aa);
            cen4 = get_centroid(isel_chain[pf4],aa);

            // 3 Vectors.
            // v12 = get_vector(cen2,cen1); // !REVERSED
            v12 = get_vector(cen1,cen2);
            v23 = get_vector(cen2,cen3);
            v34 = get_vector(cen3,cen4);

            // Middle vector length. Dimer-Dimer separation.
            m23 = magnitude(v23);

            // 3 normalized vectors.
            n12 = normalize(v12);
            n23 = normalize(v23);
            n34 = normalize(v34);
            // Can remove the normalization later, after the cross.
            // cross v23 with v12
            // cross v23 with v34


            // Normal to the Dihedral Planes.
            cross13 = cross_product(v23,v12);
            cross24 = cross_product(v23,v34);
            n13 = normalize(cross13);
            n24 = normalize(cross24);

            // Law of Sines Rule:  a x b / |a| |b|
            // The defining cross product, for 2 vectors passing by each other.
            // Right hand rule: thumb out, thumb in.
            cross14 = cross_product(n13,n24);
            sign = dot_product(cross14,n23);


            // Next step for dihedral angle calculation:
            // Get the orthogonal unit vectors.
            cos_ang = get_costheta(n12,n34);
            acos_ang = (acos(cos_ang) * 180.0) / M_PI;

            // ALTERNATE:
            // http://x3dna.org/highlights/how-to-calculate-torsion-angle
            // sign:
            // cross14 = cross_product(n12,n34);
            // sign = dot_product(cross14,n23);

            // std::cout << "The sign: " << sign << std::endl;
            // if (sign < 0)
            // {
            //     acos_ang = -1 * acos_ang;
            //     // acos_ang = -1 * acos_ang + 180.0;
            //     // std::cout << "Negative!" << std::endl;
            // }


            // else
            // {
                // acos_ang = acos_ang - 180.0;
            // }

            std::cout << "Angle(cos): " << cos_ang << std::endl;
            std::cout << "Angle(acos): " << acos_ang << std::endl;

            // sin_ang = get_sintheta(n12,n34);
            // asin_ang = asin(sin_ang) / M_PI * 180.0;
            // std::cout << "Angle(sin): " << sin_ang << std::endl;
            // std::cout << "Angle(asin): " << asin_ang << std::endl;

            // sin_ang = get_sintheta(n34,n12);
            // asin_ang = asin(sin_ang) / M_PI * 180.0;
            // std::cout << "Angle(sin): " << sin_ang << std::endl;
            // std::cout << "Angle(asin): " << asin_ang << std::endl;

            // v12.print_Vector();
            // v34.print_Vector();

            // if(sin_ang < 0)
            // {
            //     acos_ang = acos_ang * -1;
            //     fprintf(stderr,"[WARN] %s:%d: errno: %s\n",__FILE__,__LINE__,
            //             "A negative arcsin found!");
            //     std::cout << "Angle(sin): " << sin_ang << std::endl;
            //     std::cout << "Angle(asin): " << asin_ang << std::endl;
            // }

            fprintf(fp_bending_angle,"%6.1f ",acos_ang);
            fprintf(fp_bending_angle,"%6.1f ",m23);
        }

        fprintf(fp_bending_angle,"\n");
        // exit(0);

        // Dimers:
        // for(auto d: pf)
        // {
        //     std::cout << "dimer: " << d.first
        //               << " " << d.second
        //               << std::endl;
        // }
    }

    fclose(fp_bending_angle);
    // exit(0);

}


inline void SystemPF::get_bending_angle2(vAtoms aa,vIndexGroup isel_chain)
{
    // Description:
    // The bending angle is described as the 2-3 monomer-monomer vector
    // dotted with the 4-5 monomer-monomer vector. Take the acos.
    // Angles should increase from 0 degrees to 30++.
    // std::cout << "Getting Bending Angle." << std::endl;

    // FILE
    FILE * fp_bending_angle;
    fp_bending_angle = fopen("emol_mtpfbending_angle.dat", "a+");
    fprintf(fp_bending_angle,"#\n");

    // Variables:
    int num_angles; // angles (dimers * 0.5 + 1), for 13: should be 7
    int start;
    int pf1, pf2, pf3, pf4;

    Vector cen1, cen2, cen3, cen4;
    Vector v12, v23, v34;
    Vector n12, n23, n34;


    // Vector cross13, cross14, cross24;
    // Vector n13, n24;
    // double cos_ang, acos_ang, sign;

    double m23;
    Vector c12, c23, on12, on23, n1b2;
    double dihedral_rad, dihedral_deg, yv, xv;

    // std::cout << "Num_Protofilaments: " << num_protofilaments << std::endl;
    // std::cout << "Num_Dimers: " << protofilaments[0].size() << std::endl;
    if(protofilaments[0].size() == 12)
    {

    }
    else if(protofilaments[0].size() == 8)
    {

    }
    else
    {
        fprintf(stderr,"[WARN] %s:%d: errno: %s\n",__FILE__,__LINE__,
                "Protofilaments' bending angles were not computed.");
        return;
    }
    num_angles = protofilaments[0].size() * 0.5 + 1;
    start = (protofilaments[0].size() - num_angles) * 0.5; // 12->2, 8->1 (pos)



    for(auto pf: protofilaments)
    {

        for(int s=start; s < start + num_angles; s++)
        {
            // Identify 4 monomers.
            pf1 = pf[s].first;
            pf2 = pf[s].second;
            pf3 = pf[s + 1].first;
            pf4 = pf[s + 1].second;

            // 4 Centroids.
            cen1 = get_centroid(isel_chain[pf1],aa);
            cen2 = get_centroid(isel_chain[pf2],aa);
            cen3 = get_centroid(isel_chain[pf3],aa);
            cen4 = get_centroid(isel_chain[pf4],aa);

            // 3 Vectors.
            v12 = get_vector(cen1,cen2);
            v23 = get_vector(cen2,cen3);
            v34 = get_vector(cen3,cen4);

            // 3 normalized vectors.
            n12 = normalize(v12);
            n23 = normalize(v23); // b2
            n34 = normalize(v34);

            // Middle vector length. Dimer-Dimer separation.
            m23 = magnitude(v23);


            // Normal to the Dihedral Planes. (orthonormal).
            // https://math.stackexchange.com/questions/47059/
            // how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates
            c12 = cross_product(v12,v23);
            on12 = normalize(c12);        // n1

            c23 = cross_product(v23,v34);
            on23 = normalize(c23);        // n2

            // on12, n23,
            n1b2 = cross_product(on12,n23); // m1

            xv = dot_product(n12,n23);
            yv = dot_product(n1b2,on23);

            dihedral_rad = atan2(yv,xv);
            dihedral_deg = dihedral_rad * 180.0 / M_PI;
            // std::cout << "Angle: " << dihedral_deg << std::endl;

            fprintf(fp_bending_angle,"%6.1f ",dihedral_deg);
            fprintf(fp_bending_angle,"%6.1f ",m23);

        }

        fprintf(fp_bending_angle,"\n");

    }

    fclose(fp_bending_angle);

}


inline void SystemPF::get_distance_centroid(vAtoms aa,vIndexGroup isel_chain)
{
    // Description:


    // FILE
    FILE * fp_dist_centroid;
    fp_dist_centroid = fopen("emol_mtpfdist_centroid.dat","a+");
    fprintf(fp_dist_centroid,"#\n");

    // Variables:
    // int num_angles; // angles (dimers * 0.5 + 1), for 13: should be 7
    // int start;
    int pf1, pf2, pf3, pf4;
    Vector cen1, cen2, cen3, cen4;
    Vector v23;

    // Vector v12, v23, v34;
    // Vector n12, n23, n34;

    double m23;
    // Vector c12, c23, on12, on23, n1b2;
    // double dihedral_rad, dihedral_deg, yv, xv;

    // std::cout << "Num_Protofilaments: " << num_protofilaments << std::endl;
    // std::cout << "Num_Dimers: " << protofilaments[0].size() << std::endl;
    // std::cout << num_protofilaments << std::endl;
    // std::cout << protofilaments[0].size() << std::endl;

    if(protofilaments[0].size() == 12)
    {

    }
    else if(protofilaments[0].size() == 8)
    {

    }
    else
    {
        fprintf(stderr,"[WARN] %s:%d: errno: %s\n",__FILE__,__LINE__,
                "Protofilaments' bending angles were not computed.");
        return;
    }
    // num_angles = protofilaments[0].size() * 0.5 + 1;
    // start = (protofilaments[0].size() - num_angles) * 0.5; // 12->2, 8->1 (pos)




    for(auto pf: protofilaments)
    {

        for(int s=0; s < protofilaments[0].size()-1; s++)
        {
        //     // Identify 4 monomers.
            pf1 = pf[s].first;
            pf2 = pf[s].second;
            pf3 = pf[s + 1].first;
            pf4 = pf[s + 1].second;

            // std::cout << pf1 << " " << pf2 << " " << pf3 << " "<< pf4 << std::endl;

        //     // 4 Centroids.
        //     // cen1 = get_centroid(isel_chain[pf1],aa);
            cen2 = get_centroid(isel_chain[pf2],aa);
            cen3 = get_centroid(isel_chain[pf3],aa);
        //     // cen4 = get_centroid(isel_chain[pf4],aa);

        //     // 3 Vectors.
        //     // v12 = get_vector(cen1,cen2);
            v23 = get_vector(cen2,cen3);
        //     // v34 = get_vector(cen3,cen4);

        //     // 3 normalized vectors.
        //     // n12 = normalize(v12);
        //     // n23 = normalize(v23); // b2
        //     // n34 = normalize(v34);

        //     // Middle vector length. Dimer-Dimer separation.
            m23 = magnitude(v23);

        //     // Normal to the Dihedral Planes. (orthonormal).
        //     // https://math.stackexchange.com/questions/47059/
        //     // how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates
        //     // c12 = cross_product(v12,v23);
        //     // on12 = normalize(c12);        // n1

        //     // c23 = cross_product(v23,v34);
        //     // on23 = normalize(c23);        // n2

        //     // // on12, n23,
        //     // n1b2 = cross_product(on12,n23); // m1

        //     // xv = dot_product(n12,n23);
        //     // yv = dot_product(n1b2,on23);

        //     // dihedral_rad = atan2(yv,xv);
        //     // dihedral_deg = dihedral_rad * 180.0 / M_PI;
        //     // std::cout << "Angle: " << dihedral_deg << std::endl;

            fprintf(fp_dist_centroid,"%3d %3d ",pf2,pf3);
            fprintf(fp_dist_centroid,"%6.1f ",m23);
            fprintf(fp_dist_centroid,"%6.1f %6.1f %6.1f ",cen2.x,cen2.y,cen2.z);
            fprintf(fp_dist_centroid,"%6.1f %6.1f %6.1f ",cen3.x,cen3.y,cen3.z);
            fprintf(fp_dist_centroid,"\n");

        //     fprintf(fp_dist_centroid,"%6.1f %6.1f %6.1f ",cen3.x,cen3.y,cen3.z);
        }

        // fprintf(fp_dist_centroid,"\n");

    }

    fclose(fp_dist_centroid);

}


inline void SystemPF::get_beta_angle(vAtoms aa,vIndexGroup isel_chain)
{
    // Description:


    // FILE
    FILE * fp_beta_angle;
    fp_beta_angle = fopen("emol_mtpf_beta_angle.dat","a+");
    fprintf(fp_beta_angle,"#\n");

    // Variables:
    // int num_angles; // angles (dimers * 0.5 + 1), for 13: should be 7
    // int start;
    int pf1, pf2, pf3, pfa;
    Vector cen1, cen2, cen3, cen4;
    Vector n12,n23;
    Vector v12,v23,va;
    double m,vdot,acos_vdot;


    // std::cout << "Num_Protofilaments: " << num_protofilaments << std::endl;
    // std::cout << "Num_Dimers: " << protofilaments[0].size() << std::endl;
    // std::cout << num_protofilaments << std::endl;
    // std::cout << protofilaments[0].size() << std::endl;

    if(protofilaments[0].size() == 12)
    {

    }
    else if(protofilaments[0].size() == 8)
    {

    }
    else
    {
        fprintf(stderr,"[WARN] %s:%d: errno: %s\n",__FILE__,__LINE__,
                "Protofilaments' bending angles were not computed.");
        return;
    }
    // num_angles = protofilaments[0].size() * 0.5 + 1;
    // start = (protofilaments[0].size() - num_angles) * 0.5; // 12->2, 8->1 (pos)




    for(auto pf: protofilaments)
    {

        // for(int s=1; s < protofilaments[0].size()-2; s++)
        for(int s=1; s < pf.size()-1; s++)
        {
            // Identify 4 monomers.
            pf1 = pf[s-1].second; // beta -1
            pf2 = pf[s].second;   // beta
            pf3 = pf[s + 1].second; // beta+1

            pfa = pf[s + 1].first;    // alpha


            // 4 Centroids.
            cen1 = get_centroid(isel_chain[pf1],aa);
            cen2 = get_centroid(isel_chain[pf2],aa);
            cen3 = get_centroid(isel_chain[pf3],aa);
            cen4 = get_centroid(isel_chain[pfa],aa);

            // 3 Vectors.
            v12 = get_vector(cen2,cen1);
            v23 = get_vector(cen3,cen2);
            va  = get_vector(cen4,cen2); // alpha
            m   = magnitude(va);

            // 3 normalized vectors.
            n12 = normalize(v12);
            n23 = normalize(v23); // b2


            vdot = dot_product(n12,n23); // scalar.
            acos_vdot = (acos(vdot)) * 180.0 / M_PI;




            fprintf(fp_beta_angle,"%7.2f %7.2f ",m,acos_vdot);
            // fprintf(fp_beta_angle,"%6.1f %6.1f %6.1f ",cen2.x,cen2.y,cen2.z);
            // fprintf(fp_beta_angle,"%6.1f %6.1f %6.1f ",cen3.x,cen3.y,cen3.z);
            // fprintf(fp_beta_angle,"\n");

        //     fprintf(fp_beta_angle,"%6.1f %6.1f %6.1f ",cen3.x,cen3.y,cen3.z);
        }

        // fprintf(fp_beta_angle,"#");
        fprintf(fp_beta_angle,"\n");

    }

    fclose(fp_beta_angle);

}


#endif
