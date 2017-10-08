/*this is the file main_angle_tot.c*/

/* The command line is: run <datafile> (i.e. the PDB structure file)*/
/*compile with gcc -O3 main_angle_tot.c angle.c -lm -Wall -o run*/

/*this program first reads the structure and the secondary structure of a
  given protein directly from the PDB file. Then it determines the Phi and the
  Psi angles for each point of the backbone.*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <malloc.h>
#include <math.h>
#include "define.h"

typedef struct Aminoacid{
    char Name[4], SS; /*amino acid name and secondary structure*/
    int  ind,angle;
    /*angle < 10 is good and angle = 10 is bad (acc. to Ramachandran's plot)*/
}Amino;

static Amino   AminoA[L_max];

static double  *D1, *D2, *D3; /*coord of each atom in the protein*/
static char    **Type, **Code; /*Type=type of atom (e.g. CA or OXT)*/
                               /*Code=3 letter code of the amino acid*/

static char    **Number1_h,**Number2_h,**Index;
static char    **Number1_s,**Number2_s;

static void    fatal(char *);
static void    fatalf(char *,char *);
static FILE    *ckopen(char *,char *);
static char    *ckalloc(int);

static FILE    *fp;

int main(int argc, char *argv[])
{
    char   buffer[200], *s, **symbol, filename[20];
    int    count,count_h,count_s,i,j,k,pos,N; /*N = # residues in protein*/
    int    n_alpha,n_beta,mism_alpha,mism_beta;
    int    mism_a[N_class+1],mism_b[N_class+1],nr_a[N_class+1],nr_b[N_class+1];
    int    start_h[H_max],end_h[H_max],start_s[S_max],end_s[S_max],Ind_H,Ind_S;
    int    seq_pos[Atom_max], Ind_pos, string_to_integer();
    double *xA, *yA, *zA; /*coord of the CA atoms with test prot.*/
    double *xN, *yN, *zN; /*coord of the N atoms with test prot.*/
    double *xC, *yC, *zC; /*coord of the C atoms with test prot.*/
    double Phi[L_max], Psi[L_max];
    void   free_vector(), angle_Phi(), angle_Psi(), identify();

    /*corresp. between number and 3-letter symbol for a.a.*/
    symbol = (char **)ckalloc((N_class+1)*sizeof(char *));

    for(i=1;i<=N_class;i++) symbol[i] = (char *)ckalloc(4*sizeof(char));

    symbol[1]  = "CYS";
    symbol[2]  = "PHE";
    symbol[3]  = "LEU";
    symbol[4]  = "TRP";
    symbol[5]  = "VAL";
    symbol[6]  = "ILE";
    symbol[7]  = "MET";
    symbol[8]  = "HIS";
    symbol[9]  = "TYR";
    symbol[10] = "ALA";
    symbol[11] = "GLY";
    symbol[12] = "PRO";
    symbol[13] = "ASN";
    symbol[14] = "THR";
    symbol[15] = "SER";
    symbol[16] = "ARG";
    symbol[17] = "GLN";
    symbol[18] = "ASP";
    symbol[19] = "LYS";
    symbol[20] = "GLU";

    /*read the 3D structure file with the coord. of the test prot.*/
    /*read the helices*/
    count_h = 0;
    fp = ckopen(argv[1],"r");
    while(fgets(buffer,200,fp) && (buffer[0] != '\n')){
        if(!strncmp(buffer,"HELIX",5)){
            count_h++;
        }
    }

    if(count_h != 0){
        Number1_h = (char **)ckalloc(count_h*sizeof(char *));
        Number2_h = (char **)ckalloc(count_h*sizeof(char *));
    }

    rewind(fp);

    if(count_h != 0){

        i = 0;
        while(fgets(buffer,200,fp) && (buffer[0] != '\n')){
            if(!strncmp(buffer,"HELIX",5)){

                s = buffer;

                s += 22;

                Number1_h[i] = (char *)ckalloc(4*sizeof(char));
                for(pos=0;pos<3;Number1_h[i][pos++]=*s++);
                Number1_h[i][pos] = '\0';

                s += 9;

                Number2_h[i] = (char *)ckalloc(4*sizeof(char));
                for(pos=0;pos<3;Number2_h[i][pos++]=*s++);
                Number2_h[i][pos] = '\0';

                i++;
            }

        }
    }
    fclose(fp);

    /*read the sheets*/
    count_s = 0;
    fp = ckopen(argv[1],"r");

    while(fgets(buffer,200,fp) && (buffer[0] != '\n')){
        if(!strncmp(buffer,"SHEET",5)){
            count_s++;
        }
    }

    if(count_s != 0){
        Number1_s = (char **)ckalloc(count_s*sizeof(char *));
        Number2_s = (char **)ckalloc(count_s*sizeof(char *));
    }

    rewind(fp);

    if(count_s != 0){

        i = 0;
        while(fgets(buffer,200,fp) && (buffer[0] != '\n')){
            if(!strncmp(buffer,"SHEET",5)){

                s = buffer;

                s += 23;

                Number1_s[i] = (char *)ckalloc(4*sizeof(char));
                for(pos=0;pos<3;Number1_s[i][pos++]=*s++);
                Number1_s[i][pos] = '\0';

                s += 8;

                Number2_s[i] = (char *)ckalloc(4*sizeof(char));
                for(pos=0;pos<3;Number2_s[i][pos++]=*s++);
                Number2_s[i][pos] = '\0';

                i++;
            }

        }
    }
    fclose(fp);

    /*read the positions of the amino acids in the sequence*/
    count = 0;
    fp = ckopen(argv[1],"r");

    while(fgets(buffer,200,fp) && (buffer[0] != '\n')){
        if(!strncmp(buffer,"ATOM",4)){
            count++;
        }
    }

    Index     = (char **)ckalloc((count+1)*sizeof(char *));
    D1        = (double *)ckalloc((count+1)*sizeof(double));
    D2        = (double *)ckalloc((count+1)*sizeof(double));
    D3        = (double *)ckalloc((count+1)*sizeof(double));
    Type      = (char **)ckalloc((count+1)*sizeof(char *));
    Code      = (char **)ckalloc((count+1)*sizeof(char *));

    rewind(fp);

    j = 0;
    while(j<count){

        fgets(buffer,200,fp);

        if(!strncmp(buffer,"ATOM",4)){

            s = buffer;

            s += 13;

            pos = 0;
            Type[j+1] = (char *)ckalloc(4*sizeof(char));
            while(*s != ' ') Type[j+1][pos++] = *s++;
            Type[j+1][pos] = '\0';

            while(*s==' ') s++;

            Code[j+1] = (char *)ckalloc(4*sizeof(char));
            for(pos=0;pos<3;Code[j+1][pos++]=*s++);
            Code[j+1][3] = '\0';

            s += 3;

            Index[j+1] = (char *)ckalloc(5*sizeof(char));
            for(pos=0;pos<3;Index[j+1][pos++]=*s++);
            Index[j+1][pos] = '\0';

            s += 5;

            sscanf(s, "%lf %lf %lf", &D1[j+1], &D2[j+1], &D3[j+1]);

            j++;

        }

    }
    fclose(fp);

    /*------------end read the structure of the protein---------------------*/

    printf("for protein %s\n",argv[1]);

    /*----check what was read-------*/
    fp=fopen("check","a");

    fprintf(fp,"\n\nfor protein %s\n",argv[1]);

    Ind_H = 0;
    for(i=0;i<count_h;i++){
        if(Number1_h[i] != NULL){
            fprintf(fp,"%5d %s\n", i, Number1_h[i]);
            Ind_H++;
            start_h[Ind_H] = string_to_integer(Number1_h[i]);
            fprintf(fp,"%5d start_h = %5d\n", Ind_H, start_h[Ind_H]);
        }
        if(Number2_h[i] != NULL){
            fprintf(fp,"%5d %s\n", i, Number2_h[i]);
            end_h[Ind_H] = string_to_integer(Number2_h[i]);
            fprintf(fp,"%5d end_h = %5d\n", Ind_H, end_h[Ind_H]);
        }
    }

    Ind_S = 0;
    for(i=0;i<count_s;i++){
        if(Number1_s[i] != NULL){
            fprintf(fp,"%5d %s\n", i, Number1_s[i]);
            Ind_S++;
            start_s[Ind_S] = string_to_integer(Number1_s[i]);
            fprintf(fp,"%5d start_s = %5d\n", Ind_S, start_s[Ind_S]);
        }

        if(Number2_s[i] != NULL){
            fprintf(fp,"%5d %s\n", i, Number2_s[i]);
            end_s[Ind_S] = string_to_integer(Number2_s[i]);
            fprintf(fp,"%5d end_s = %5d\n", Ind_S, end_s[Ind_S]);
        }
    }

    Ind_pos = 0;
    for(i=1;i<=count;i++){
        if(Index[i] != NULL){
            fprintf(fp,"%5d %s\n", i, Index[i]);
            if(i==1){
                Ind_pos++;
                seq_pos[Ind_pos] = string_to_integer(Index[i]);
                fprintf(fp,"%5d seq_pos = %5d\n", Ind_pos, seq_pos[Ind_pos]);
            }

            /*monitor the position only for the whole a.a.*/
            else{
                if(!strncmp(Index[i],Index[i-1],3));
                else{
                    Ind_pos++;
                    seq_pos[Ind_pos] = string_to_integer(Index[i]);
                    fprintf(fp,"%5d seq_pos = %5d\n", Ind_pos, seq_pos[Ind_pos]);
                }
            }
        }
    }

    for(i=1;i<=Ind_pos;i++){
        AminoA[i].SS = 'T';
    }

    for(i=1;i<=Ind_pos;i++){
        if(Ind_H != 0){
            for(j=1;j<=Ind_H;j++){
                if((seq_pos[i] > start_h[j]) && (seq_pos[i] < end_h[j])){
                    /*if this helix is longer than 6 residues*/
                    if((end_h[j]- start_h[j] + 1) >= 6)
                        AminoA[i].SS = 'H';
                }
            }
        }

        if(Ind_S != 0){
            for(j=1;j<=Ind_S;j++){
                if((seq_pos[i] >= start_s[j]) && (seq_pos[i] <= end_s[j])){
                    AminoA[i].SS = 'E';
                }
            }
        }

        fprintf(fp,"%5d %5d %c\n", i, seq_pos[i], AminoA[i].SS);
    }
    fclose(fp);

    printf("length from S.S. is %5d\n",Ind_pos);

    xA = (double *) ckalloc((Ind_pos+1)*sizeof(double));
    yA = (double *) ckalloc((Ind_pos+1)*sizeof(double));
    zA = (double *) ckalloc((Ind_pos+1)*sizeof(double));
    xC = (double *) ckalloc((Ind_pos+1)*sizeof(double));
    yC = (double *) ckalloc((Ind_pos+1)*sizeof(double));
    zC = (double *) ckalloc((Ind_pos+1)*sizeof(double));
    xN = (double *) ckalloc((Ind_pos+1)*sizeof(double));
    yN = (double *) ckalloc((Ind_pos+1)*sizeof(double));
    zN = (double *) ckalloc((Ind_pos+1)*sizeof(double));

    /*initialize the coord of the main chain, CB and build atoms*/
    for(i=1;i<=Ind_pos;i++){
        xA[i] = 0.0;
        yA[i] = 0.0;
        zA[i] = 0.0;
        xN[i] = 0.0;
        yN[i] = 0.0;
        zN[i] = 0.0;
        xC[i] = 0.0;
        yC[i] = 0.0;
        zC[i] = 0.0;
    }

    /*identify the coord. of the main chain(CA,N,C) atoms*/
    N=0;
    for(i=1;i<=count;i++){
        if(Type[i] != NULL){
            if(!strcmp(Type[i],"CA")){
                N++;
                identify(D1[i],D2[i],D3[i],&xA[N],&yA[N],&zA[N]);
                strcpy(AminoA[N].Name,Code[i]);
                if(!strcmp(Code[i],"CYS")) AminoA[N].ind = 1;
                else if(!strcmp(Code[i],"PHE")) AminoA[N].ind = 2;
                else if(!strcmp(Code[i],"LEU")) AminoA[N].ind = 3;
                else if(!strcmp(Code[i],"TRP")) AminoA[N].ind = 4;
                else if(!strcmp(Code[i],"VAL")) AminoA[N].ind = 5;
                else if(!strcmp(Code[i],"ILE")) AminoA[N].ind = 6;
                else if(!strcmp(Code[i],"MET")) AminoA[N].ind = 7;
                else if(!strcmp(Code[i],"HIS")) AminoA[N].ind = 8;
                else if(!strcmp(Code[i],"TYR")) AminoA[N].ind = 9;
                else if(!strcmp(Code[i],"ALA")) AminoA[N].ind = 10;
                else if(!strcmp(Code[i],"GLY")) AminoA[N].ind = 11;
                else if(!strcmp(Code[i],"PRO")) AminoA[N].ind = 12;
                else if(!strcmp(Code[i],"ASN")) AminoA[N].ind = 13;
                else if(!strcmp(Code[i],"THR")) AminoA[N].ind = 14;
                else if(!strcmp(Code[i],"SER")) AminoA[N].ind = 15;
                else if(!strcmp(Code[i],"ARG")) AminoA[N].ind = 16;
                else if(!strcmp(Code[i],"GLN")) AminoA[N].ind = 17;
                else if(!strcmp(Code[i],"ASP")) AminoA[N].ind = 18;
                else if(!strcmp(Code[i],"LYS")) AminoA[N].ind = 19;
                else if(!strcmp(Code[i],"GLU")) AminoA[N].ind = 20;
            }
        }
    }

    printf("length from CA is %5d\n", N);

    N=0;
    for(i=1;i<=count;i++){
        if(Type[i] != NULL){
            if(!strcmp(Type[i],"N")){
                N++;
                identify(D1[i],D2[i],D3[i],&xN[N],&yN[N],&zN[N]);
            }
        }
    }

    N=0;
    for(i=1;i<=count;i++){
        if(Type[i] != NULL){
            if(!strcmp(Type[i],"C")){
                N++;
                identify (D1[i],D2[i],D3[i],&xC[N],&yC[N],&zC[N]);
            }
        }
    }


    /*initialize the Phi and the Psi angles*/
    for(i=0;i<L_max;i++){
        Phi[i] = 0.0;
        Psi[i] = 0.0;
    }

    /*determine the Phi and Psi angles*/

    for(i=1;i<=(N-1);i++){
        angle_Phi(xC[i+1],yC[i+1],zC[i+1],xA[i+1],yA[i+1],zA[i+1],
                  xN[i+1],yN[i+1],zN[i+1],xC[i],yC[i],zC[i],&Phi[i]);
        angle_Psi(xN[i],yN[i],zN[i],xA[i],yA[i],zA[i],
                  xC[i],yC[i],zC[i],xN[i+1],yN[i+1],zN[i+1],&Psi[i]);
    }

    /*printf("\nfor protein %6s\n\n", argv[1]);
      for(i=1;i<=(N-1);i++){
      printf("%3d %4s %6.1f %6.1f\n", (i+1), AminoA[i+1].Name, (Phi[i]-180.0),
      (Psi[i+1]-180.0));
      }*/

    /*determine the info about each type of a.a. in each s.s. elem. in
      terms of Phi,Psi angles (good or bad according to Ramachandran plot)*/

    for(i=1;i<=(N-1);i++){

        /*if the amino acid is in a helix*/
        if(AminoA[i+1].SS == 'H'){
            if(strcmp(AminoA[i+1].Name,"GLY") && strcmp(AminoA[i+1].Name,"PRO")){

                /*if wrong angles in helix, but not a GLY or a PRO*/
                if(((Phi[i]-180.0) < -80.0) || ((Phi[i]-180.0) > -48.0) ||
                   ((Psi[i+1]-180.0) < -59.0) || ((Psi[i+1]-180.0) > -27.0)){

                    AminoA[i+1].angle = 50;

                    sprintf(filename,"angle_alpha_hist_%s",AminoA[i+1].Name);
                    fp=fopen(filename,"a");
                    fprintf(fp,"%d\n", AminoA[i+1].angle);
                    fclose(fp);

                }

                else{

                    /*1st interval for Phi*/
                    if(((Phi[i]-180.0) >= -80.0) && ((Phi[i]-180.0) <= -70.0)){

                        if(((Psi[i+1]-180.0) >= -59.0) && ((Psi[i+1]-180.0) <= -49.0)){
                            AminoA[i+1].angle = 1;
                        }

                        else if(((Psi[i+1]-180.0) > -49.0)&&((Psi[i+1]-180.0) <= -38.0)){
                            AminoA[i+1].angle = 2;
                        }

                        else if(((Psi[i+1]-180.0) > -38.0)&&((Psi[i+1]-180.0)<= -27.0)){
                            AminoA[i+1].angle = 3;
                        }

                    }


                    /*2nd interval for Phi*/
                    else if(((Phi[i]-180.0) > -70.0)&&((Phi[i]-180.0) <= -59.0)){

                        if(((Psi[i+1]-180.0) >= -59.0) && ((Psi[i+1]-180.0) <= -49.0)){
                            AminoA[i+1].angle = 4;
                        }

                        else if(((Psi[i+1]-180.0) > -49.0)&&((Psi[i+1]-180.0) <= -38.0)){
                            AminoA[i+1].angle = 5;
                        }

                        else if(((Psi[i+1]-180.0) > -38.0)&&((Psi[i+1]-180.0)<= -27.0)){
                            AminoA[i+1].angle = 6;
                        }

                    }

                    /*3rd (last) interval for Phi*/
                    else if(((Phi[i]-180.0) > -59.0)&&((Phi[i]-180.0) <= -48.0)){

                        if(((Psi[i+1]-180.0) >= -59.0) && ((Psi[i+1]-180.0) <= -49.0)){
                            AminoA[i+1].angle = 7;
                        }

                        else if(((Psi[i+1]-180.0) > -49.0)&&((Psi[i+1]-180.0) <= -38.0)){
                            AminoA[i+1].angle = 8;
                        }

                        else if(((Psi[i+1]-180.0) > -38.0)&&((Psi[i+1]-180.0)<= -27.0)){
                            AminoA[i+1].angle = 9;
                        }

                    }

                    /*AminoA[i+1].angle = 0;*/
                    sprintf(filename,"angle_alpha_hist_%s",AminoA[i+1].Name);
                    fp=fopen(filename,"a");
                    fprintf(fp,"%d\n", AminoA[i+1].angle);
                    fclose(fp);
                }

                /*sprintf(filename,"angle_alpha_hist_%s",AminoA[i+1].Name);
                  fp=fopen(filename,"a");
                  fprintf(fp,"%7.2f %7.2f\n", Phi[i]-180.0, Psi[i+1]-180.0);
                  fclose(fp);*/

            }/*end a.a. is no GLY or PRO*/

            /*if a.a. is GLY or PRO*/
            else{

                AminoA[i+1].angle = 0;
                sprintf(filename,"angle_alpha_hist_%s",AminoA[i+1].Name);
                fp=fopen(filename,"a");
                fprintf(fp,"%d\n", AminoA[i+1].angle);
                fclose(fp);

                /*sprintf(filename,"angle_alpha_hist_%s",AminoA[i+1].Name);
                  fp=fopen(filename,"a");
                  fprintf(fp,"%7.2f %7.2f\n", Phi[i]-180.0, Psi[i+1]-180.0);
                  fclose(fp);*/
            }

        }/*end a.a. in helix*/

        /*if amino acid is in a beta sheet*/
        else if(AminoA[i+1].SS == 'E'){

            if(strcmp(AminoA[i+1].Name,"GLY") && strcmp(AminoA[i+1].Name,"PRO")){

                if(((Phi[i]-180.0) < -150.0) || ((Phi[i]-180.0) > -90.0) ||
                   ((Psi[i+1]-180.0) < 90.0) || ((Psi[i+1]-180.0) > 150.0)){

                    AminoA[i+1].angle = 50;
                    sprintf(filename,"angle_beta_hist_%s",AminoA[i+1].Name);
                    fp=fopen(filename,"a");
                    fprintf(fp,"%d\n", AminoA[i+1].angle);
                    fclose(fp);
                }
                else{

                    /*1st interval for Phi*/
                    if(((Phi[i]-180.0) >= -150.0) && ((Phi[i]-180.0) <= -140.0)){

                        if(((Psi[i+1]-180.0) >= 90.0) && ((Psi[i+1]-180.0) <= 100.0)){
                            AminoA[i+1].angle = 1;
                        }

                        else if(((Psi[i+1]-180.0) > 100.0)&&((Psi[i+1]-180.0) <= 111.0)){
                            AminoA[i+1].angle = 2;
                        }

                        else if(((Psi[i+1]-180.0) > 111.0)&&((Psi[i+1]-180.0) <= 122.0)){
                            AminoA[i+1].angle = 3;
                        }

                        else if(((Psi[i+1]-180.0) > 122.0)&&((Psi[i+1]-180.0) <= 133.0)){
                            AminoA[i+1].angle = 4;
                        }

                        else if(((Psi[i+1]-180.0) > 133.0)&&((Psi[i+1]-180.0) <= 150.0)){
                            AminoA[i+1].angle = 5;
                        }

                    }

                    /*2nd interval for Phi*/
                    else if(((Phi[i]-180.0) > -140.0)&&((Phi[i]-180.0) <= -129.0)){

                        if(((Psi[i+1]-180.0) >= 90.0) && ((Psi[i+1]-180.0) <= 100.0)){
                            AminoA[i+1].angle = 6;
                        }

                        else if(((Psi[i+1]-180.0) > 100.0)&&((Psi[i+1]-180.0) <= 111.0)){
                            AminoA[i+1].angle = 7;
                        }

                        else if(((Psi[i+1]-180.0) > 111.0)&&((Psi[i+1]-180.0) <= 122.0)){
                            AminoA[i+1].angle = 8;
                        }

                        else if(((Psi[i+1]-180.0) > 122.0)&&((Psi[i+1]-180.0) <= 133.0)){
                            AminoA[i+1].angle = 9;
                        }

                        else if(((Psi[i+1]-180.0) > 133.0)&&((Psi[i+1]-180.0) <= 150.0)){
                            AminoA[i+1].angle = 10;
                        }

                    }

                    /*3rd interval for Phi*/
                    else if(((Phi[i]-180.0) > -129.0)&&((Phi[i]-180.0) <= -118.0)){

                        if(((Psi[i+1]-180.0) >= 90.0) && ((Psi[i+1]-180.0) <= 100.0)){
                            AminoA[i+1].angle = 11;
                        }

                        else if(((Psi[i+1]-180.0) > 100.0)&&((Psi[i+1]-180.0) <= 111.0)){
                            AminoA[i+1].angle = 12;
                        }

                        else if(((Psi[i+1]-180.0) > 111.0)&&((Psi[i+1]-180.0) <= 122.0)){
                            AminoA[i+1].angle = 13;
                        }

                        else if(((Psi[i+1]-180.0) > 122.0)&&((Psi[i+1]-180.0) <= 133.0)){
                            AminoA[i+1].angle = 14;
                        }

                        else if(((Psi[i+1]-180.0) > 133.0)&&((Psi[i+1]-180.0) <= 150.0)){
                            AminoA[i+1].angle = 15;
                        }

                    }

                    /*4th interval for Phi*/
                    else if(((Phi[i]-180.0) > -118.0)&&((Phi[i]-180.0) <= -107.0)){

                        if(((Psi[i+1]-180.0) >= 90.0) && ((Psi[i+1]-180.0) <= 100.0)){
                            AminoA[i+1].angle = 16;
                        }

                        else if(((Psi[i+1]-180.0) > 100.0)&&((Psi[i+1]-180.0) <= 111.0)){
                            AminoA[i+1].angle = 17;
                        }

                        else if(((Psi[i+1]-180.0) > 111.0)&&((Psi[i+1]-180.0) <= 122.0)){
                            AminoA[i+1].angle = 18;
                        }

                        else if(((Psi[i+1]-180.0) > 122.0)&&((Psi[i+1]-180.0) <= 133.0)){
                            AminoA[i+1].angle = 19;
                        }

                        else if(((Psi[i+1]-180.0) > 133.0)&&((Psi[i+1]-180.0) <= 150.0)){
                            AminoA[i+1].angle = 20;
                        }

                    }

                    /*5th (last) interval for Phi*/
                    else if(((Phi[i]-180.0) > -107.0)&&((Phi[i]-180.0) <= -90.0)){

                        if(((Psi[i+1]-180.0) >= 90.0) && ((Psi[i+1]-180.0) <= 100.0)){
                            AminoA[i+1].angle = 21;
                        }

                        else if(((Psi[i+1]-180.0) > 100.0)&&((Psi[i+1]-180.0) <= 111.0)){
                            AminoA[i+1].angle = 22;
                        }

                        else if(((Psi[i+1]-180.0) > 111.0)&&((Psi[i+1]-180.0) <= 122.0)){
                            AminoA[i+1].angle = 23;
                        }

                        else if(((Psi[i+1]-180.0) > 122.0)&&((Psi[i+1]-180.0) <= 133.0)){
                            AminoA[i+1].angle = 24;
                        }

                        else if(((Psi[i+1]-180.0) > 133.0)&&((Psi[i+1]-180.0) <= 150.0)){
                            AminoA[i+1].angle = 25;
                        }

                    }

                    /*AminoA[i+1].angle = 0;*/
                    sprintf(filename,"angle_beta_hist_%s",AminoA[i+1].Name);
                    fp=fopen(filename,"a");
                    fprintf(fp,"%d\n", AminoA[i+1].angle);
                    fclose(fp);
                }
            }

            else{
                AminoA[i+1].angle = 0;
                sprintf(filename,"angle_beta_hist_%s",AminoA[i+1].Name);
                fp=fopen(filename,"a");
                fprintf(fp,"%d\n", AminoA[i+1].angle);
                fclose(fp);
            }

            /*sprintf(filename,"angle_beta_hist_%s",AminoA[i+1].Name);
              fp=fopen(filename,"a");
              fprintf(fp,"%7.2f %7.2f\n", Phi[i]-180.0, Psi[i+1]-180.0);
              fclose(fp);*/

        }/*end a.a. is in a sheet*/

    }/*end going over all a.a.*/

    /*----- end getting the info about angles for each type of a.a.--------*/

    /*look for mismatches in sec. struct. elem. (GLY and PRO are special!!)*/
    n_alpha = mism_alpha = 0;
    n_beta = mism_beta = 0;

    /*init. the entries for mismatches*/
    for(i=1;i<=N_class;i++){
        nr_a[i] = 0;
        nr_b[i] = 0;
        mism_a[i] = 0;
        mism_b[i] = 0;
    }

    fp = fopen("Angles","a");
    fprintf(fp,"\n\n\nfor protein %6s\n\n", argv[1]);
    for(i=1;i<=(N-1);i++){

        if(AminoA[i+1].SS == 'H'){

            n_alpha++;
            if(strcmp(AminoA[i+1].Name,"GLY") && strcmp(AminoA[i+1].Name,"PRO")){

                /*if wrong angles in helix, but not a GLY or a PRO*/
                if(((Phi[i]-180.0) < -80.0) || ((Phi[i]-180.0) > -48.0) ||
                   ((Psi[i+1]-180.0) < -59.0) || ((Psi[i+1]-180.0) > -27.0)){
                    fprintf(fp,"%3d %4s %c %6.1f %6.1f M\n",(i+1),AminoA[i+1].Name,
                            AminoA[i+1].SS,(Phi[i]-180.0),(Psi[i+1]-180.0));
                    mism_alpha++;
                    mism_a[AminoA[i+1].ind]++;
                    nr_a[AminoA[i+1].ind]++;
                }

                /*if correct angles in helix, but not a GLY or a PRO*/
                else{
                    fprintf(fp,"%3d %4s %c %6.1f %6.1f\n",(i+1),AminoA[i+1].Name,
                            AminoA[i+1].SS,(Phi[i]-180.0),(Psi[i+1]-180.0));
                    nr_a[AminoA[i+1].ind]++;
                }
            }

            /*if a GLY or a PRO*/
            else{
                fprintf(fp,"%3d %4s %c %6.1f %6.1f\n",(i+1),AminoA[i+1].Name,
                        AminoA[i+1].SS,(Phi[i]-180.0),(Psi[i+1]-180.0));
            }

        }/*end helix*/

        else if(AminoA[i+1].SS == 'E'){

            n_beta++;
            if(strcmp(AminoA[i+1].Name,"GLY") && strcmp(AminoA[i+1].Name,"PRO")){

                /*if wrong angles in sheet, but not a GLY or a PRO*/
                if(((Phi[i]-180.0) < -150.0) || ((Phi[i]-180.0) > -90.0) ||
                   ((Psi[i+1]-180.0) < 90.0) || ((Psi[i+1]-180.0) > 150.0)){
                    fprintf(fp,"%3d %4s %c %6.1f %6.1f M\n",(i+1),AminoA[i+1].Name,
                            AminoA[i+1].SS,(Phi[i]-180.0),(Psi[i+1]-180.0));
                    mism_beta++;
                    mism_b[AminoA[i+1].ind]++;
                    nr_b[AminoA[i+1].ind]++;
                }

                /*if correct angles in helix, but not a GLY or a PRO*/
                else{
                    fprintf(fp,"%3d %4s %c %6.1f %6.1f\n",(i+1),AminoA[i+1].Name,
                            AminoA[i+1].SS,(Phi[i]-180.0),(Psi[i+1]-180.0));
                    nr_b[AminoA[i+1].ind]++;
                }
            }

            /*if a GLY or a PRO*/
            else{
                fprintf(fp,"%3d %4s %c %6.1f %6.1f\n",(i+1),AminoA[i+1].Name,
                        AminoA[i+1].SS,(Phi[i]-180.0),(Psi[i+1]-180.0));
            }

        }/*end beta-sheet*/

    }/*end go over all residues in the protein*/

    fp = fopen("mismatch_in_angles","a");
    fprintf(fp,"\n\n\nfor protein %6s\n\n", argv[1]);
    if(n_alpha != 0){
        fprintf(fp,"in helices = %3.1f (%4d %4d)\n",
                ((double) mism_alpha)/((double) n_alpha)*100.,mism_alpha,n_alpha);
    }
    if(n_beta != 0){
        fprintf(fp,"in strands = %3.1f (%4d %4d)\n",
                ((double) mism_beta)/((double) n_beta)*100.,mism_beta,n_beta);
    }
    fclose(fp);

    /*fp=fopen("mismatch_alpha_aa","a");
      for(i=1;i<=N_class;i++){
      if(nr_a[i] !=0){
      fprintf(fp,"%3d %4s %5.1f\n", i, symbol[i],
      (((double) mism_a[i])/(nr_a[i]))*100.);
      }
      }
      fclose(fp);

      fp=fopen("mismatch_beta_aa","a");
      for(i=1;i<=N_class;i++){
      if(nr_b[i] !=0){
      fprintf(fp,"%3d %4s %5.1f\n", i, symbol[i],
      (((double) mism_b[i])/(nr_b[i]))*100.);
      }
      }
      fclose(fp);

      fp=fopen("mismatch_alpha_angle_aa","a");
      for(i=1;i<=N_class;i++){
      if(nr_a[i] !=0){
      fprintf(fp,"%3d %3d %3d\n", i, mism_a[i], nr_a[i]);
      }
      }
      fclose(fp);

      fp=fopen("mismatch_beta_angle_aa","a");
      for(i=1;i<=N_class;i++){
      if(nr_b[i] !=0){
      fprintf(fp,"%3d %3d %3d\n", i, mism_b[i], nr_b[i]);
      }
      }
      fclose(fp);*/

    /* free all */
    free(D1); free(D2); free(D3);
    for (i=1; i<=count; i++) {
        free(Type[i]);
        free(Code[i]);
        free(Index[i]);
    }

    for (i=0; i<count_h; i++) {
        free(Number1_h[i]);
        free(Number2_h[i]);
    }

    for (i=0; i<count_s; i++) {
        free(Number1_s[i]);
        free(Number2_s[i]);
    }

    free(Code); free(Type);
    free(Number1_h); free(Number2_h); free(Number1_s); free(Number2_s);
    free(Index);

    /*deallocate the space of the positions for the test prot.*/
    free_vector(xA,L_max);
    free_vector(yA,L_max);
    free_vector(zA,L_max);
    free_vector(xN,L_max);
    free_vector(yN,L_max);
    free_vector(zN,L_max);
    free_vector(xC,L_max);
    free_vector(yC,L_max);
    free_vector(zC,L_max);

    exit(0);
}

/* ------------------ library of utility routines -------------------- */

/* fatal - print message and die */
void fatal(char *msg)
{
    fprintf(stderr, "%s\n", msg);
    exit(1);
}


/* fatalf - format message, print it, and die */
void fatalf(char *msg, char *val)
{
    fprintf(stderr, msg, val);
    putc('\n', stderr);
    exit(1);
}

/* ckopen - open file; check for success */
FILE *ckopen(char *name, char *mode)
{
    FILE *fopen(), *fp;

    if ((fp = fopen(name, mode)) == NULL)
        fatalf("Cannot open %s.", name);
    return fp;
}
/* ckalloc - allocate space; check for success */
char *ckalloc(int amount)
{
    char *malloc(), *p;

    if ((p = malloc( (unsigned) amount)) == NULL)
        fatal("Ran out of memory.");
    return p;
}

void nrerror(error_text)
     char error_text[];
     /* Numerical Recipes standard error handler */
{
    void _exit();

    fprintf(stderr,"Numerical Recipes run-time error...\n");
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"...now exiting to system...\n");
    _exit(1);
}


void free_vector(v,nh)
     double *v;
     int nh;
     /* free a double vector allocated with vector() */
{
    free((char*) (v));
}

void identify(d1,d2,d3,x,y,z)
    double d1,d2,d3,*x,*y,*z;
{
    *x = d1;
    *y = d2;
    *z = d3;
}


/*convert a string into an integer*/
int string_to_integer(char string[])
{
    int i, integer_value, result = 0;

    i = 0;
    while(string[i] != '\0'){

        if(string[i] >= '0' && string[i] <= '9'){
            integer_value = string[i] - '0';
            result = result*10 + integer_value;
        }

        i++;

    }

    return(result);

}
