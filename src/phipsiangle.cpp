// phipsiangle.cpp

/* ---------------------------------------------------------
   standard libraries
   --------------------------------------------------------- */
// #include <stdio.h>
// #include <stdlib.h> // strtod?, stod
// #include <assert.h>
// #include <string.h> // string.
// #include <cctype>
// #include <algorithm> // remove_if, count
// #include <iostream>
// #include <fstream>
// #include <ctime>
// #include <list>        // std::list
// #include <iterator> // istream_iterator
#include <cmath>
// #include <vector>
// #include "boost/tuple/tuple.hpp"


/* ---------------------------------------------------------
   headers
   --------------------------------------------------------- */
#include "debug.h"
#include "md.h"
#include "system.hpp"
// #include "microtubule.hpp"
// #include "dcd.h"
// #include "dcdio.h"

/* ---------------------------------------------------------
   functions
   --------------------------------------------------------- */
void compute_phipsi(Dihedral dh)
{
    std::cout << "Computing phi / psi." << std::endl;

    // for(auto a: bb)
    // {
    //     std::cout << a.x << std::endl;
    // }

    for(auto d: dh)
    {

        std::cout << d[0].atomtype;
        std::cout << ": " << d[0].x   << " " << d[0].y   << " " << d[0].z;
        std::cout << std::endl;
        std::cout << d[1].atomtype;
        std::cout << ": " << d[1].x << " " << d[1].y << " " << d[1].z;
        std::cout << std::endl;
        std::cout << d[2].atomtype;
        std::cout << ": " << d[2].x << " " << d[2].y << " " << d[2].z;
        std::cout << std::endl;
        std::cout << d[3].atomtype;
        // std::cout << ": " << d[3].x << " " << d[3].y << " " << d[3].z;
        // std::cout << std::endl;
        // std::cout << std::endl;

        std::cout << "Max-x:" << Max(d[0].x,d[1].x) << std::endl;
        std::cout << "Max-y:" << Max(d[0].y,d[1].y) << std::endl;
    }
    // N: -0.905 -0.978 -0.56
    // CA: 0.01 -0.05 0.08
    // C: 1.463 -0.49 -0.078
    // O: 2.366 0.329 -0.305 -- not needed.

    // Atom phi1, phi2, phi3;
    // Atom psi1, psi2, psi3;

    // 5 positions from residues: r1-CO, r2-N, r2-CA, r2-CO, r3-N
    // phi: r1-CO -- r2-N -- r2-CA -- r2-CO   | h1, h2, h3, h4
    // psi: r2-N  -- r2-CA-- r2-CO -- r3-N    | h2, h3, h4, h5
    Vector h1, h2, h3, h4, h5;

    Vector phi1, phi2; // 2 vectors to form the phi plane
    // Vector vec1, vec2; // 2 vectors for the middle plane
    Vector psi1, psi2; // 2 vectors to form the psi plane

    Vector ortho_phi, ortho_psi, ortho_mid;  // orthogonal vectors.
    Vector nphi, npsi, nmid;                 // normalized, orthogonal vector to plane.
    double cosphi, cospsi;                   // cosine(angles)
    double phi,psi;                          // angles.

    for(int i=1; i<dh.size()-1; i++)
    {
        // phi1 = dh[i-1][2]; // resid 1, C
        // phi2 = dh[i][0];   // resid 2, N
        // phi3 = dh[i][1];   // resid 2, CA
        h1.create_Vector(dh[i-1][2].x,
                         dh[i-1][2].y,
                         dh[i-1][2].z);
        h2.create_Vector(dh[i][0].x,
                         dh[i][0].y,
                         dh[i][0].z);
        h3.create_Vector(dh[i][1].x,
                         dh[i][1].y,
                         dh[i][1].z);

        // psi1 = dh[i][1];   // resid 2, CA
        // psi2 = dh[i][2];   // resid 2, C
        // psi3 = dh[i+1][0]; // resid 3, N
        // s1.create_Vector(dh[i][1].x,
        //                  dh[i][1].y,
        //                  dh[i][1].z);
        h4.create_Vector(dh[i][2].x,
                         dh[i][2].y,
                         dh[i][2].z);
        h5.create_Vector(dh[i+1][0].x,
                         dh[i+1][0].y,
                         dh[i+1][0].z);


        // get_vector: 2nd - 1st
        phi1 = get_vector(h1,h2);
        phi2 = get_vector(h2,h3);

        psi1 = get_vector(h3,h4);
        psi2 = get_vector(h4,h5);

        // Get 3 planes:
        ortho_phi = cross_product(phi1,phi2);
        ortho_mid = cross_product(psi1,phi2);
        ortho_psi = cross_product(psi1,psi2);

        // Normalize:
        nphi = normalize(ortho_phi);
        npsi = normalize(ortho_psi);
        nmid = normalize(ortho_mid);

        // Get cos-angles:
        cosphi = dot_product(nphi,nmid);
        cospsi = dot_product(npsi,nmid);

        // Get angles:
        phi = acos(cosphi);
        psi = acos(cospsi);

        // std::cout << phi1.atomtype;
        // std::cout << ": " << phi1.x   << " " << phi1.y   << " " << phi1.z;
        // std::cout << std::endl;
        // std::cout << phi2.atomtype;
        // std::cout << ": " << phi2.x   << " " << phi2.y   << " " << phi2.z;
        // std::cout << std::endl;
        // std::cout << phi3.atomtype;
        // std::cout << ": " << phi3.x   << " " << phi3.y   << " " << phi3.z;
        // std::cout << std::endl;




    }




    // for(int i=0; i<bb.size(); i+=4)
    // {
    //     std::cout << i << std::endl;
    //     std::cout << "0: " << bb[i].x   << " " << bb[i].y   << " " << bb[i].z;
    //     std::cout << std::endl;
    //     std::cout << "1: " << bb[i+1].x << " " << bb[i+1].y << " " << bb[i+1].z;
    //     std::cout << std::endl;
    //     std::cout << "2: " << bb[i+2].x << " " << bb[i+2].y << " " << bb[i+2].z;
    //     std::cout << std::endl;
    //     std::cout << "3: " << bb[i+3].x << " " << bb[i+3].y << " " << bb[i+3].z;
    //     std::cout << std::endl;
    // }


}
