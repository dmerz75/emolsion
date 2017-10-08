/*this is the file angle_var.c*/

#include <stdio.h>
#include <malloc.h>
#include <math.h>

static FILE *fp;

/*this function determines the Phi and Psi angles between 2 amino acids.
  It builds the normal to the plane of the first amino acid = the vector n, where
  the plane of the 1st amino acid is determined by the cartesian coordinates of
  the C_{i}(x1,y1,z1), CA_{i}(x2,y2,z2) and N_{i}(x3,y3,z3) atoms. Then it
  builds the normal from C_{i+1}(x0,y0,z0) to the plane. If the parameter t>0
  then C_{i+1} is under the plane of the 1st amino acid, while for t<0 it is
  above the plane.*/

void angle_Phi(x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,z0,Phi)
    double x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,z0,*Phi;
{
    double p1x,p1y,p1z,p2x,p2y,p2z,dx,dy,dz;
    double n1x,n2x,n1y,n2y,n1z,n2z,n1,n2, cosinus;
    double px,py,pz,cs;

    p1x = x0-x2;
    p1y = y0-y2;
    p1z = z0-z2;

    dx = x3-x2;
    dy = y3-y2;
    dz = z3-z2;

    p2x = x3-x1;
    p2y = y3-y1;
    p2z = z3-z1;

    n1x = dy*p1z - dz*p1y;
    n1y = dz*p1x - dx*p1z;
    n1z = dx*p1y - dy*p1x;

    n1 = sqrt(n1x*n1x + n1y*n1y + n1z*n1z);

    n2x = dy*p2z - dz*p2y;
    n2y = dz*p2x - dx*p2z;
    n2z = dx*p2y - dy*p2x;

    n2 = sqrt(n2x*n2x + n2y*n2y + n2z*n2z);

    cosinus = (n1x*n2x + n1y*n2y + n1z*n2z)/(n1*n2);

    *Phi = acos(cosinus)*90./acos(0.);

    px = -n1y*n2z + n1z*n2y;
    py = -n1z*n2x + n1x*n2z;
    pz = -n1x*n2y + n1y*n2x;

    cs = px*dx + py*dy + pz*dz;

    if(cs < 0.)
        *Phi = 360. - *Phi;
}


void angle_Psi(x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,z0,Psi)
    double x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,z0,*Psi;
{
    double p1x,p1y,p1z,p2x,p2y,p2z,dx,dy,dz;
    double n1x,n2x,n1y,n2y,n1z,n2z,n1,n2, cosinus;
    double px,py,pz,cs;

    p1x = x0-x2;
    p1y = y0-y2;
    p1z = z0-z2;

    dx = x3-x2;
    dy = y3-y2;
    dz = z3-z2;

    p2x = x3-x1;
    p2y = y3-y1;
    p2z = z3-z1;

    n1x = dy*p1z - dz*p1y;
    n1y = dz*p1x - dx*p1z;
    n1z = dx*p1y - dy*p1x;

    n1 = sqrt(n1x*n1x + n1y*n1y + n1z*n1z);

    n2x = dy*p2z - dz*p2y;
    n2y = dz*p2x - dx*p2z;
    n2z = dx*p2y - dy*p2x;

    n2 = sqrt(n2x*n2x + n2y*n2y + n2z*n2z);

    cosinus = (n1x*n2x + n1y*n2y + n1z*n2z)/(n1*n2);

    *Psi = acos(cosinus)*90./acos(0.);

    px = n2y*n1z - n2z*n1y;
    py = n2z*n1x - n2x*n1z;
    pz = n2x*n1y - n2y*n1x;

    cs = (px*dx + py*dy + pz*dz);

    if(cs < 0.)
        *Psi = 360. - *Psi;

}
