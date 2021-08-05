#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Matrice.h"
#include "EOS_etain.h"



// *******   Newton pour calculer le point triple  ********
/* Valeur d'entree :
    Va,Ea, Vb,Eb, Vc,Ec

    corespodances :  1-a   2-b   3-c

  Valeur sortie
   VA,EA, VB,EB, VC,EC
*/
void Newton_triple(double* A, double* InvA, double Va, double Ea, double Vb, double Eb, double Vc, double Ec, double* VA, double* EA, double* VB, double* EB, double* VC, double* EC){
  int dim=6;
  // Grandeurs thermodynamiques
  double Pa, Ta, Sa, Ga, dPva, dPea, dTva, dTea;
  double Pb, Tb, Sb, Gb, dPvb, dPeb, dTvb, dTeb;
  double Pc, Tc, Sc, Gc, dPvc, dPec, dTvc, dTec;
  // phase 1
  Pa=fP(Va,Ea, 1); Ta=fT(Va,Ea, 1);
  Sa=fS(Va,Ea, 1); Ga=Ea+Pa*Va-Ta*Sa;
  dPva=dPdV(Va,1); dPea=dPdE(1);
  dTva=dTdV(Va,1); dTea=dTdE(1);
  //phase 2
  Pb=fP(Vb,Eb, 2); Tb=fT(Vb,Eb, 2);
  Sb=fS(Vb,Eb, 2); Gb=Eb+Pb*Vb-Tb*Sb;
  dPvb=dPdV(Vb,2); dPeb=dPdE(2);
  dTvb=dTdV(Vb,2); dTeb=dTdE(2);
  // phase 3
  Pc=fP(Vc,Ec, 3); Tc=fT(Vc,Ec, 3);
  Sc=fS(Vc,Ec, 3); Gc=Ec+Pc*Vc-Tc*Sc;
  dPvc=dPdV(Vc,3); dPec=dPdE(3);
  dTvc=dTdV(Vc,3); dTec=dTdE(3);

  double dGva, dGea, dGvb, dGeb, dGvc ,dGec;
  dGva=Va*dPva - Sa*dTva; 
  dGea=Va*dPea - Sa*dTea;

  dGvb=Vb*dPvb - Sb*dTvb; 
  dGeb=Vb*dPeb - Sb*dTeb;

  dGvc=Vc*dPvc - Sc*dTvc; 
  dGec=Vc*dPec - Sc*dTec;


  double* Delta=malloc(6*sizeof(double));
  Delta[0]=Pb-Pa; Delta[1]=Pc-Pa;
  Delta[2]=Tb-Ta; Delta[3]=Tc-Ta;
  Delta[4]=Gb-Ga; Delta[5]=Gc-Ga; 
  
  // on met le second menbre dans la dernière colonne
  //matrice A
  A[dim*0+0]=dPva;  A[dim*0+1]=dPea; A[dim*0+2]=-dPvb;  A[dim*0+3]=-dPeb;  A[dim*0+4]=0.;     A[dim*0+5]=0.;  

  A[dim*1+0]=dPva;  A[dim*1+1]=dPea; A[dim*1+2]=0.;     A[dim*1+3]=0.;     A[dim*1+4]=-dPvc;  A[dim*1+5]=-dPec; 

  A[dim*2+0]=dTva;  A[dim*2+1]=dTea; A[dim*2+2]=-dTvb;  A[dim*2+3]=-dTeb;  A[dim*2+4]=0.;     A[dim*2+5]=0.; 

  A[dim*3+0]=dTva;  A[dim*3+1]=dTea; A[dim*3+2]=0.;     A[dim*3+3]=0.;     A[dim*3+4]=-dTvc;  A[dim*3+5]=-dTec; 

  A[dim*4+0]=dGva;  A[dim*4+1]=dGea; A[dim*4+2]=-dGvb;  A[dim*4+3]=-dGeb;  A[dim*4+4]=0.;     A[dim*4+5]=0.; 

  A[dim*5+0]=dGva;  A[dim*5+1]=dGea; A[dim*5+2]=0.;     A[dim*5+3]=0.;     A[dim*5+4]=-dGvc;  A[dim*5+5]=-dGec;
  

  double det_a=inverse_matrice_pivot(A, dim, InvA);
  printf("detA= %.15lf\n",det_a );
  if ( fabs(det_a)<1e-8){printf("A non inversible : detA= %.15lf\n", det_a);}
  
  /*
  printf("  * InvA = *\n");
  for(int i=0; i<6; i++){
    for(int j=0; j<6; j++){
      printf("InvA[%d][%d]= %.15lf  ",i,j, InvA[dim*i+j] );
    }
    printf("\n");
  }
  printf("\n");
  
  
  // Mutiplcation matricielle de A et InvA
  double res;
  printf("  * InvA*A = *\n");
  for(int i=0; i<6; i++){
    for(int j=0; j<6; j++){
      res=0;
      for(int k=0; k<6; k++){
        res+=InvA[dim*i+k]*A[dim*k+j];
      }
      printf("A*InvA[%d][%d]= %.15lf ",i,j,res );
    }
    printf("\n");
  }
  printf("\n");
 */


  // Mutiplcation matricielle
  double* dX=malloc(6*sizeof(double));

  for(int i=0; i<6; i++){
    dX[i]=0;
    for(int j=0; j<6; j++){
      dX[i]+=InvA[dim*i+j]*Delta[j];
    }
    //printf("  dX[%d]= %.15lf\n",i,dX[i] );
  }
  //printf("\n");
   
  double correc=1.;

  *VA=Va+correc*dX[0];
  *EA=Ea+correc*dX[1];
  *VB=Vb+correc*dX[2];
  *EB=Eb+correc*dX[3];
  *VC=Vc+correc*dX[4];
  *EC=Ec+correc*dX[5];


  free(Delta);
  free(dX);
}





// *** Calcul le point triple en lancant le newton  ****
/*
   Renvoie les 3 couples (V,E) du triangle
*/
void triple(int Nmax, double epsilon, double* Va0, double* Ea0, double* Vb0, double* Eb0, double* Vc0, double* Ec0){
  int n=0;
  int dim=6;

  double *A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));

  double VA, EA, VB, EB, VC, EC;

  double Pa, Ta, Sa, Ga, Ea, Va;
  double Pb, Tb, Sb, Gb, Eb, Vb;
  double Pc, Tc, Sc, Gc, Ec, Vc;

  Pa=fP(*Va0,*Ea0,1);

  // phase 1
  Ea=*Ea0; Va=*Va0;
  Pa=fP(Va,Ea, 1); Ta=fT(Va,Ea, 1);
  Sa=fS(Va,Ea, 1); Ga=Ea+Pa*Va-Ta*Sa;
  //phase 2
  Eb=*Eb0; Vb=*Vb0;
  Pb=fP(Vb,Eb, 2); Tb=fT(Vb,Eb, 2);
  Sb=fS(Vb,Eb, 2); Gb=Eb+Pb*Vb-Tb*Sb;
  // phase 3
  Ec=*Ec0; Vc=*Vc0;
  Pc=fP(Vc,Ec, 3); Tc=fT(Vc,Ec, 3);
  Sc=fS(Vc,Ec, 3); Gc=Ec+Pc*Vc-Tc*Sc;


  double critere=fabs(Ga-Gb)+fabs(Ga-Gc);

  while(critere>=epsilon && n<Nmax){
    printf("n= %d\n",n );
    Newton_triple(A, InvA ,Va, Ea, Vb, Eb, Vc, Ec, &VA, &EA, &VB, &EB, &VC, &EC);
    
    /*
    printf(" dVa= %.15lf, VA= %.15lf (cm^3/kg)\n",(VA-Va)*1e6, VA*1e6 );
    printf(" dEa= %.15lf, EA= %.15lf (kJ/kg)\n",(EA-Ea)*1e-3, EA*1e-3 );
    printf(" dVb= %.15lf, VB= %.15lf (cm^3/kg)\n",(VB-Vb)*1e6, VB*1e6 );
    printf(" dEb= %.15lf, EB= %.15lf (kJ/kg)\n",(EB-Eb)*1e-3, EB*1e-3 );
    printf(" dVc= %.15lf, VC= %.15lf (cm^3/kg)\n",(VC-Vc)*1e6, VC*1e6 );
    printf(" dEc= %.15lf, EC= %.15lf (kJ/kg)\n",(EC-Ec)*1e-3, EC *1e-3);
    */

    Va=VA; Ea=EA;
    Vb=VB; Eb=EB;
    Vc=VC; Ec=EC;

    // phase 1
    Pa=fP(Va,Ea, 1); Ta=fT(Va,Ea, 1);
    Sa=fS(Va,Ea, 1); Ga=Ea+Pa*Va-Ta*Sa;
    //phase 2
    Pb=fP(Vb,Eb, 2); Tb=fT(Vb,Eb, 2);
    Sb=fS(Vb,Eb, 2); Gb=Eb+Pb*Vb-Tb*Sb;
    // phase 3
    Pc=fP(Vc,Ec, 3); Tc=fT(Vc,Ec, 3);
    Sc=fS(Vc,Ec, 3); Gc=Ec+Pc*Vc-Tc*Sc;
 
 /*
 /////////////////////////////////////////////////////////////////////
    printf(" Iteration (triple) :\n");
    printf("   Ta= %.15lf\n",Ta );
    printf("   Tb= %.15lf\n",Tb );
    printf("   Tc= %.15lf\n",Tc );
    printf("\n");
    printf("   Pa= %.15lf (Gpa)\n",Pa*1e-9 );
    printf("   Pb= %.15lf (Gpa)\n",Pb*1e-9 );
    printf("   Pc= %.15lf (Gpa)\n",Pc*1e-9 );
    printf("\n");
    printf("   Ga= %.15lf\n",Ga );
    printf("   Gb= %.15lf\n",Gb );
    printf("   Gc= %.15lf\n",Gc );

    printf("   Sa= %.15lf, Sb= %.15lf, Sc= %.15lf\n",Sa, Sb, Sc );
  */
/////////////////////////////////////////////////////////////////////


    critere=fabs(Ga-Gb)+fabs(Ga-Gc);
    printf(" critere= %.14lf\n",critere );

    n++;
  }


  *Va0=Va; *Ea0=Ea; 
  *Vb0=Vb; *Eb0=Eb; 
  *Vc0=Vc; *Ec0=Ec;

  free(A); free(InvA);

}



// *******  Newton sur les lignes de ghgmt de phase  ********

/* 
  Utile dans la fonction  ligne_double

  Valeur d'entree :
   Va,Ea, Vb,Eb, Vc,Ec

   corespodances :  1-a   2-b   3-c

  Valeur sortie
   VA,EA, VB,EB, VC,EC

*/
void Newton_double(double* A, double* InvA, double Pcst, int phaseA, int phaseB, double Va, double Ea, double Vb, double Eb, double* VA, double* EA, double* VB, double* EB){
  int dim=4;

  // Grandeurs thermodynamiques
  double Pa, Ta, Sa, Ga, dPva, dPea, dTva, dTea;
  double Pb, Tb, Sb, Gb, dPvb, dPeb, dTvb, dTeb;
  // phase A
  Pa=fP(Va,Ea, phaseA); Ta=fT(Va,Ea, phaseA);
  Sa=fS(Va,Ea, phaseA); Ga=Ea+Pa*Va-Ta*Sa;
  dPva=dPdV(Va,phaseA); dPea=dPdE(phaseA);
  dTva=dTdV(Va,phaseA); dTea=dTdE(phaseA);
  // phase B
  Pb=fP(Vb,Eb, phaseB); Tb=fT(Vb,Eb, phaseB);
  Sb=fS(Vb,Eb, phaseB); Gb=Eb+Pb*Vb-Tb*Sb;
  dPvb=dPdV(Vb,phaseB); dPeb=dPdE(phaseB);
  dTvb=dTdV(Vb,phaseB); dTeb=dTdE(phaseB);


  double dGva, dGea, dGvb, dGeb;
  dGva=Va*dPva - Sa*dTva; 
  dGea=Va*dPea - Sa*dTea;

  dGvb=Vb*dPvb - Sb*dTvb; 
  dGeb=Vb*dPeb - Sb*dTeb;


  double* Delta=malloc(dim*sizeof(double));
  Delta[0]=Pa-Pcst; Delta[1]=Pb-Pcst;
  Delta[2]=Tb-Ta;   Delta[3]=Gb-Ga;
   
  /*
  for(int j=0; j<dim; j++){
    printf("Delta[%d]= %.15lf  ",j, Delta[j] );
  }
  printf("\n");
  */


  // on met le second menbre dans la dernière colonne
  //matrice A
  A[dim*0+0]=-dPva; A[dim*0+1]=-dPea; A[dim*0+2]=0.;    A[dim*0+3]=0.;  

  A[dim*1+0]=0.;    A[dim*1+1]=0.;    A[dim*1+2]=-dPvb; A[dim*1+3]=-dPeb;      

  A[dim*2+0]=dTva;  A[dim*2+1]=dTea;  A[dim*2+2]=-dTvb; A[dim*2+3]=-dTeb;  

  A[dim*3+0]=dGva;  A[dim*3+1]=dGea;  A[dim*3+2]=-dGvb; A[dim*3+3]=-dGeb; 
  
  /*
  printf("  * A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      printf("A[%d][%d]= %.15lf  ",i,j, A[dim*i+j] );
    }
    printf("\n");
  }
  printf("\n");  
  */


  double det_a=inverse_matrice_pivot(A, dim, InvA);
  //printf("detA= %.15lf\n",det_a );
  if ( fabs(det_a)<1e-8){printf("A non inversible : detA= %.15lf\n", det_a);}
  
  /*
  printf("  * InvA = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      printf("InvA[%d][%d]= %.15lf  ",i,j, InvA[dim*i+j] );
    }
    printf("\n");
  }
  printf("\n");
  */
  
  /*
  // Mutiplcation matricielle de A et InvA
  double res;
  printf("  * InvA*A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      res=0;
      for(int k=0; k<dim; k++){
        res+=InvA[dim*i+k]*A[dim*k+j];
      }
      printf("A*InvA[%d][%d]= %.15lf ",i,j,res );
    }
    printf("\n");
  }
  printf("\n");
  */


  // Mutiplcation matricielle
  double* dX=malloc(6*sizeof(double));

  for(int i=0; i<6; i++){
    dX[i]=0;
    for(int j=0; j<6; j++){
      dX[i]+=InvA[dim*i+j]*Delta[j];
    }
    //printf("  dX[%d]= %.15lf\n",i,dX[i] );
  }
  //printf("\n");
   
  double correc=1.;

  *VA=Va+correc*dX[0];
  *EA=Ea+correc*dX[1];
  *VB=Vb+correc*dX[2];
  *EB=Eb+correc*dX[3];


  free(Delta);
  free(dX);
}


//  *******  Fonction qui calcule les points (V,E) des lignes de chgmt de phases
/*
Arguments d'entree :
  - Nmax : nombre max d'iteration
  - epsilon : tolérance sur le critère  (ici fabs(G1-G2))
  - [Pdeb, Pfin] : segment de presion sur lequel calcule la zone mixte
  - Np : nomnbre de points de discretisatoin du segment
  - Va0, Ea0, Vb0, Eb0 : Valeur initiale de (Va,Ea) et (Vb,Eb)
   (typiquement les coord du point triple pour chaque phase)

Arguments de sortie :
  - (segP, segT) : courbe de T(P) du duagramme de phase (frontières)
  - (segVa, segEa) : courbe pour la phase A
  - (segVb, segEb) : courbe pour la phase B

*/
void ligne_double(int Nmax, double epsilon, int phaseA, int phaseB, double Pdeb, double Pfin, int Np, double Va0, double Ea0, double Vb0, double Eb0, double* segP, double* segT, double* segVa, double* segEa, double* segVb, double* segEb){
  int n=1;
  int dim=4;
  double Pcst;

  double *A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));

  double hp=(Pfin-Pdeb)/(Np-1);

  double VA, EA, VB, EB;

  double Pa, Ta, Sa, Ga, Ea, Va;
  double Pb, Tb, Sb, Gb, Eb, Vb;

  Ea=Ea0; Va=Va0;
  Eb=Eb0; Vb=Vb0;


  double critere;
  
  for(int i=0; i<Np; i++){
    //printf("\n****  i= %d ******\n",i );
    segP[i]=Pdeb+i*hp;
    Pcst=segP[i];
    //printf("  Pcst= %.15lf, hp= %.15lf, Pdeb= %.15lf\n", Pcst, hp, Pdeb );
    
    //printf("  ok 1\n");
    // phase A
    Pa=fP(Va,Ea, phaseA); Ta=fT(Va,Ea, phaseA);
    Sa=fS(Va,Ea, phaseA); Ga=Ea+Pa*Va-Ta*Sa;
    //phase B
    //printf("  ok 2\n");
    Pb=fP(Vb,Eb, phaseB); Tb=fT(Vb,Eb, phaseB);
    Sb=fS(Vb,Eb, phaseB); Gb=Eb+Pb*Vb-Tb*Sb;
    
    //printf("  ok 3\n");

    critere=fabs(Ga-Gb)+fabs(Pa-Pcst)+fabs(Pb-Pcst);
    //printf("  critere= %.14lf\n",critere );
    
    n=1;
    while(critere>=epsilon && n<Nmax){
      //printf("  --n= %d\n",n );

      Newton_double(A, InvA, Pcst, phaseA, phaseB, Va, Ea, Vb, Eb, &VA, &EA, &VB, &EB);
      
      /*
      printf(" dVa= %.15lf, VA= %.15lf (cm^3/kg)\n",(VA-Va)*1e6, VA*1e6 );
      printf(" dEa= %.15lf, EA= %.15lf (kJ/kg)\n",(EA-Ea)*1e-3, EA*1e-3 );
      printf(" dVb= %.15lf, VB= %.15lf (cm^3/kg)\n",(VB-Vb)*1e6, VB*1e6 );
      printf(" dEb= %.15lf, EB= %.15lf (kJ/kg)\n",(EB-Eb)*1e-3, EB*1e-3 );
      printf(" dVc= %.15lf, VC= %.15lf (cm^3/kg)\n",(VC-Vc)*1e6, VC*1e6 );
      printf(" dEc= %.15lf, EC= %.15lf (kJ/kg)\n",(EC-Ec)*1e-3, EC *1e-3);
      */

      Va=VA; Ea=EA;
      Vb=VB; Eb=EB;

      
      // phase A
      Pa=fP(Va,Ea, phaseA); Ta=fT(Va,Ea, phaseA);
      Sa=fS(Va,Ea, phaseA); Ga=Ea+Pa*Va-Ta*Sa;
      //phase B
      Pb=fP(Vb,Eb, phaseB); Tb=fT(Vb,Eb, phaseB);
      Sb=fS(Vb,Eb, phaseB); Gb=Eb+Pb*Vb-Tb*Sb;

      critere=fabs(Ga-Gb)+fabs(Pa-Pcst)+fabs(Pb-Pcst);
      //printf("    critere= %.14lf\n",critere );

      /*
   /////////////////////////////////////////////////////////////////////
      printf(" Iteration (triple) :\n");
      printf("   Ta= %.15lf\n",Ta );
      printf("   Tb= %.15lf\n",Tb );
      printf("   Tc= %.15lf\n",Tc );
      printf("\n");
      printf("   Pa= %.15lf (Gpa)\n",Pa*1e-9 );
      printf("   Pb= %.15lf (Gpa)\n",Pb*1e-9 );
      printf("   Pc= %.15lf (Gpa)\n",Pc*1e-9 );
      printf("\n");

    */
  /////////////////////////////////////////////////////////////////////

      
      n++;
    }
    //printf("  n fin =%d \n",n);
    //printf("  ok 4\n");

    segT[i]=Tb;
    segVa[i]=Va;
    segEa[i]=Ea;
    segVb[i]=Vb;
    segEb[i]=Eb;

    //printf("  ok 5\n");
    
    /*
    printf("  Iteration (ligne_double) :\n");
    printf("    Ea= %.15lf\n",Ea );
    printf("    Eb= %.15lf\n",Eb );
    printf("\n");
    printf("    Va= %.15lf\n",Va*1e6 );
    printf("    Vb= %.15lf\n",Vb*1e6 );
    printf("\n");
    printf("    Ta= %.15lf\n",Ta );
    printf("    Tb= %.15lf\n",Tb );
    printf("\n");
    printf("    Pcst= %.15lf (GPa)\n",Pcst*1e-9 );
    printf("    Pa= %.15lf (GPa)\n",Pa*1e-9 );
    printf("    Pb= %.15lf (GPa)\n",Pb*1e-9 );
    printf("\n");
    printf("    Ga= %.15lf\n",Ga );
    printf("    Gb= %.15lf\n",Gb );
    */

  }

  //printf("  ok fin boucle i\n");
  free(A); 
  free(InvA);

}

//  ****** Fonction pour imprimer dans des fichiers les points des frontières de domaine P,T  et   V,E
/*
   Utile dans la fonction diag_phase  qui calcul ces points
*/
int print_double(int Np12, double* segP12, double* segT12, double* segVa12, double* segEa12, double* segVb12, double* segEb12,
                 int Np13, double* segP13, double* segT13, double* segVa13, double* segEa13, double* segVb13, double* segEb13,
                 int Np23, double* segP23, double* segT23, double* segVa23, double* segEa23, double* segVb23, double* segEb23,
                 int Nb, int Nh, int Ng, int Nd, double Tb, double Th, double Pg, double Pd,
                 double* segPBAS, double* segPHAUT, double* segTGAUCHE, double* segTDROITE){
  
  FILE* fileP12, *fileT12, *fileVa12, *fileEa12, *fileVb12, *fileEb12;
  FILE* fileP13, *fileT13, *fileVa13, *fileEa13, *fileVb13, *fileEb13;
  FILE* fileP23, *fileT23, *fileVa23, *fileEa23, *fileVb23, *fileEb23;
  FILE* filePb, *filePh, *fileTg, *fileTd;


  if((fileP12 = fopen("P12.txt", "w+")) == NULL){printf("erreur ouverture fichier P\n");return 1;}
  if((fileT12 = fopen("T12.txt", "w+")) == NULL){printf("erreur ouverture fichier T\n");return 1;}
  if((fileVa12 = fopen("Va12.txt", "w+")) == NULL){printf("erreur ouverture fichier Va\n");return 1;}
  if((fileEa12 = fopen("Ea12.txt", "w+")) == NULL){printf("erreur ouverture fichier Ea\n");return 1;}
  if((fileVb12 = fopen("Vb12.txt", "w+")) == NULL){printf("erreur ouverture fichier Vb\n");return 1;}
  if((fileEb12 = fopen("Eb12.txt", "w+")) == NULL){printf("erreur ouverture fichier Eb\n");return 1;}

  if((fileP13 = fopen("P13.txt", "w+")) == NULL){printf("erreur ouverture fichier P\n");return 1;}
  if((fileT13 = fopen("T13.txt", "w+")) == NULL){printf("erreur ouverture fichier T\n");return 1;}
  if((fileVa13 = fopen("Va13.txt", "w+")) == NULL){printf("erreur ouverture fichier Va\n");return 1;}
  if((fileEa13 = fopen("Ea13.txt", "w+")) == NULL){printf("erreur ouverture fichier Ea\n");return 1;}
  if((fileVb13 = fopen("Vb13.txt", "w+")) == NULL){printf("erreur ouverture fichier Vb\n");return 1;}
  if((fileEb13 = fopen("Eb13.txt", "w+")) == NULL){printf("erreur ouverture fichier Eb\n");return 1;}

  if((fileP23 = fopen("P23.txt", "w+")) == NULL){printf("erreur ouverture fichier P\n");return 1;}
  if((fileT23 = fopen("T23.txt", "w+")) == NULL){printf("erreur ouverture fichier T\n");return 1;}
  if((fileVa23 = fopen("Va23.txt", "w+")) == NULL){printf("erreur ouverture fichier Va\n");return 1;}
  if((fileEa23 = fopen("Ea23.txt", "w+")) == NULL){printf("erreur ouverture fichier Ea\n");return 1;}
  if((fileVb23 = fopen("Vb23.txt", "w+")) == NULL){printf("erreur ouverture fichier Vb\n");return 1;}
  if((fileEb23 = fopen("Eb23.txt", "w+")) == NULL){printf("erreur ouverture fichier Eb\n");return 1;}

  if((filePb = fopen("fbas.txt", "w+")) == NULL){printf("erreur ouverture fichier P\n");return 1;}
  if((filePh = fopen("fhaut.txt", "w+")) == NULL){printf("erreur ouverture fichier P\n");return 1;}
  if((fileTg = fopen("fgauche.txt", "w+")) == NULL){printf("erreur ouverture fichier T\n");return 1;}
  if((fileTd = fopen("fdroite.txt", "w+")) == NULL){printf("erreur ouverture fichier T\n");return 1;}


  for(int i=0; i<Np12; i++){
    fprintf(fileP12, "%.15lf ",segP12[i] );
    fprintf(fileT12, "%.15lf ",segT12[i] );
    fprintf(fileVa12, "%.15lf ",segVa12[i] );
    fprintf(fileEa12, "%.15lf ",segEa12[i] );
    fprintf(fileVb12, "%.15lf ",segVb12[i] );
    fprintf(fileEb12, "%.15lf ",segEb12[i] );  
  }

  for(int i=0; i<Np13; i++){
    fprintf(fileP13, "%.15lf ",segP13[i] );
    fprintf(fileT13, "%.15lf ",segT13[i] );
    fprintf(fileVa13, "%.15lf ",segVa13[i] );
    fprintf(fileEa13, "%.15lf ",segEa13[i] );
    fprintf(fileVb13, "%.15lf ",segVb13[i] );
    fprintf(fileEb13, "%.15lf ",segEb13[i] );  
  }

  for(int i=0; i<Np23; i++){
    fprintf(fileP23, "%.15lf ",segP23[i] );
    fprintf(fileT23, "%.15lf ",segT23[i] );
    fprintf(fileVa23, "%.15lf ",segVa23[i] );
    fprintf(fileEa23, "%.15lf ",segEa23[i] );
    fprintf(fileVb23, "%.15lf ",segVb23[i] );
    fprintf(fileEb23, "%.15lf ",segEb23[i] );  
  }

  for(int i=0; i<Nb; i++){
    fprintf(filePb, "%15.lf %.15lf\n",segPBAS[i],Tb );
  }
  for(int i=0; i<Nh; i++){
    fprintf(filePh, "%15.lf %.15lf\n",segPHAUT[i], Th );
  }
  for(int i=0; i<Ng; i++){
    fprintf(fileTg, "%15.lf %.15lf\n",Pg, segTGAUCHE[i] );
  }
  for(int i=0; i<Nd; i++){
    fprintf(fileTd, "%15.lf %.15lf\n",Pd, segTDROITE[i] );
  }



  fclose(fileP12);   fclose(fileT12);   
  fclose(fileVa12);  fclose(fileEa12);   
  fclose(fileVb12);  fclose(fileEb12);   

  fclose(fileP13);   fclose(fileT13);   
  fclose(fileVa13);  fclose(fileEa13);   
  fclose(fileVb13);  fclose(fileEb13);   

  fclose(fileP23);   fclose(fileT23);   
  fclose(fileVa23);  fclose(fileEa23);   
  fclose(fileVb23);  fclose(fileEb23);   

}



// *********    fonction d'inversion    ************
/*
    Obtention de (V,E) à partir de (P,T)
    NEwton classique 
*/
int fVE(int phase, double P, double T, double* pV, double* pE){
  int Nmax=1e2;
  double epsilon=1e-6;
  double ur=1.2;
  int n=0;
  double f, fprime;
  double V,E, dV;
  
  double K0, N0, gamma0, Crv, theta0, T0, P0, rho0, v0, E0, Sr;
  coeff(phase, &K0, &N0, &gamma0, &Crv, &theta0, &T0, &P0, &rho0, &v0, &E0, &Sr);


  // valuer initiale
  if(phase==1){
    V=129.905671875844291e-6;
    E=72.168190127265447e3;
  }
  else if(phase==2){
    V=127.666112185063994e-6;
    E=99.480256119236515e3;
  }
  else if(phase==3){
    V=132.538478805442338e-6;
    E=136.154346228414312e3;
  }
  else{printf("erreur phase= %d (fVE)\n",phase );}

  
  double critere=1.;
  
  // Newton sur E
  while(critere>epsilon && n<Nmax){
    f=fPs(V, phase) + gamma0*rho0*Crv*(T-ur*THETA(V, phase))-P;
    fprime=Ps_prime(V,phase) - gamma0*rho0*Crv*ur*gamma0*exp(gamma0*fx(V,phase));
    
    dV=f/fprime;
    V-=dV;
    critere=fabs(dV);
    //printf("n= %d, V= %.15lf, dV= %.15lf\n",n,V,dV );
    n++;
  }
  
  E = fEs(V,phase) + (1./(rho0*gamma0))*(P-fPs(V,phase));

  *pV=V;
  *pE=E;

  return 0;
}



// *********    Calcul du Point courbe frontières   ************
/*
   Dans le domaine (P,T)    ET    dans le domaine (V,E)
*/
int diag_phase(void){ 

  int Nmax=1e4;
  double epsilon=5e-5;
  double nb_points=20;
  int Np12=2*nb_points; int Np13=nb_points; int Np23=6*nb_points;
  int Nh=3*nb_points;  
  int Ng1=nb_points;    int Ng2=4*nb_points; 
  int Nd1=4*nb_points;  int Nd2=0.5*nb_points; 
  int Nb1=2*nb_points;  int Nb2=4*nb_points; 
  int Ng=Ng1+Ng2;
  int Nd=Nd1+Nd2;
  int Nb=Nb1+Nb2;

  double Ptriple=4.5104217816091e9;
  double Pdeb=Ptriple;


  int phaseA12=1;
  int phaseB12=2;
  double Pfin12=9.4e9;
  
  int phaseA13=1;
  int phaseB13=3;
  double Pfin13=0.;
  
  int phaseA23=2;
  int phaseB23=3;
  double Pfin23=70e9;
  
  double Va012, Vb012, Va013, Vb013, Va023, Vb023;
  double Ea012, Eb012, Ea013, Eb013, Ea023, Eb023;


  // Valeur d'initialisation
  Va012=129.905671875844291e-6;  Ea012=72.168190127265447e3; 
  Vb012=127.666112185063994e-6;  Eb012=99.480256119236515e3;
  
  Va013=129.905671875844291e-6;  Ea013=72.168190127265447e3; 
  Vb013=132.538478805442338e-6;  Eb013=136.154346228414312e3;

  Va023=127.666112185063994e-6;  Ea023=99.480256119236515e3; 
  Vb023=132.538478805442338e-6;  Eb023=136.154346228414312e3;


   
  // Initialisation des tableaux
  double* segP12=malloc(Np12*sizeof(double));
  double* segT12=malloc(Np12*sizeof(double));
  double* segVa12=malloc(Np12*sizeof(double));
  double* segEa12=malloc(Np12*sizeof(double));
  double* segVb12=malloc(Np12*sizeof(double));
  double* segEb12=malloc(Np12*sizeof(double));

  double* segP13=malloc(Np13*sizeof(double));
  double* segT13=malloc(Np13*sizeof(double));
  double* segVa13=malloc(Np13*sizeof(double));
  double* segEa13=malloc(Np13*sizeof(double));
  double* segVc13=malloc(Np13*sizeof(double));
  double* segEc13=malloc(Np13*sizeof(double));

  double* segP23=malloc(Np23*sizeof(double));
  double* segT23=malloc(Np23*sizeof(double));
  double* segVb23=malloc(Np23*sizeof(double));
  double* segEb23=malloc(Np23*sizeof(double));
  double* segVc23=malloc(Np23*sizeof(double));
  double* segEc23=malloc(Np23*sizeof(double));

  double* segPHAUT=malloc(Nh*sizeof(double));
  double* segPBAS=malloc((Nb1+Nb2)*sizeof(double));
  double* segTGAUCHE=malloc((Ng1+Ng2)*sizeof(double));
  double* segTDROITE=malloc((Nd1+Nd2)*sizeof(double));

  // frontiere 12
  //printf("PHASE 12 :\n");
  ligne_double(Nmax, epsilon, phaseA12, phaseB12, Pdeb, Pfin12, Np12, Va012, Ea012, Vb012, Eb012, segP12, segT12, segVa12, segEa12, segVb12, segEb12);

  // frontiere 13
  //printf("PHASE 13 :\n");
  ligne_double(Nmax, epsilon, phaseA13, phaseB13, Pdeb, Pfin13, Np13, Va013, Ea013, Vb013, Eb013, segP13, segT13, segVa13, segEa13, segVc13, segEc13);
  
  // frontiere 23
  //printf("PHASE 23 :\n");
  ligne_double(Nmax, epsilon, phaseA23, phaseB23, Pdeb, Pfin23, Np23, Va023, Ea023, Vb023, Eb023, segP23, segT23, segVb23, segEb23, segVc23, segEc23);
  
  
  //frontières extérieures :
  double Pg=0.;
  double Pd=70e9;
  double Tb=300;
  double Th=2500;

  double Td=2486.138962087622531;
  double Tg=505;
  double Pb=9.4e9;
  // HAUT
  for(int i=0; i<Nh; i++){
    segPHAUT[i]=Pg+i*(Pd-Pg)/(Nh-1);
  }

  // BAS
  for(int i=0; i<Nb1; i++){
    segPBAS[i]=Pg+i*(Pb-Pg)/(Nb1-1);
  }
  for(int i=0; i<Nb2; i++){
    segPBAS[Nb1+i]=Pb+(i+1)*(Pd-Pb)/(Nb2);
  }

  // GAUCHE
  for(int i=0; i<Ng1; i++){
    segTGAUCHE[i]=Tb+i*(Tg-Tb)/(Ng1-1);
  }
  for(int i=0; i<Ng2; i++){
    segTGAUCHE[Ng1+i]=Tg+(i+1)*(Th-Tg)/(Ng2);
  }

  // DROITE
  for(int i=0; i<Nd1; i++){
    segTDROITE[i]=Tb+i*(Td-Tb)/(Nd1-1);
  }
  for(int i=0; i<Nd2; i++){
    segTDROITE[i+Nd1]=Td+(i+1)*(Th-Td)/(Nd2);
  }


  /*    fonction print   */
  int res= print_double(Np12, segP12, segT12, segVa12, segEa12, segVb12, segEb12,  
                        Np13, segP13, segT13, segVa13, segEa13, segVc13, segEc13,    
                        Np23, segP23, segT23, segVb23, segEb23, segVc23, segEc23,
                        Nb, Nh, Ng, Nd, Tb, Th, Pg, Pd,
                        segPBAS, segPHAUT, segTGAUCHE, segTDROITE);    
  if(res){printf("Erreur dans l'impression (diag_phase)"); return 1;}




  int err;
  int phase;
  double P;
  double T;
  double E,V;
  
  FILE* fres1;
  FILE* fres2;
  FILE* fres3;
  if((fres1 = fopen("VEfront1.txt", "w+")) == NULL){printf("erreur ouverture fichier (routine_VE)\n");return 1;}
  if((fres2 = fopen("VEfront2.txt", "w+")) == NULL){printf("erreur ouverture fichier (routine_VE)\n");return 1;}
  if((fres3 = fopen("VEfront3.txt", "w+")) == NULL){printf("erreur ouverture fichier (routine_VE)\n");return 1;}

  

  //  ### Calcul et impression des frontières dans le domaine V,E
  int N;
  double Tdeb, Tfin;
  double Pfin;

  // frontière 1/A
  phase=1;
  N=100;
  P=0;
  Tdeb=505;
  Tfin=300;
  for(int i=0; i<N; i++){
    T=Tdeb + i*(Tfin-Tdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    printf("%.15lf %.15lf\n",V*1e6,E );
    fprintf(fres1, "%.15lf %.15lf\n",V,E );
  }
  T=300;
  Pdeb=0;
  Pfin=9.4e9;
  for(int i=0; i<N; i++){
    P=Pdeb + i*(Pfin-Pdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    printf("%.15lf %.15lf\n",V*1e6,E );
    fprintf(fres1, "%.15lf %.15lf\n",V,E );
  }
  fclose(fres1);



  // frontière 2/B
  phase=2;
  N=100;
  T=300;
  Pdeb=9.4e9;
  Pfin=70e9;
  for(int i=0; i<N; i++){
    P=Pdeb + i*(Pfin-Pdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    printf("%.5lf %.5lf\n",V*1e6,E );
    fprintf(fres2, "%.15lf %.15lf\n",V,E );
  }
  P=70e9;
  Tdeb=300.;
  Tfin=2486.138962087622531;
  for(int i=0; i<N; i++){
    T=Tdeb + i*(Tfin-Tdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    printf("%.5lf %.5lf\n",V*1e6,E );
    fprintf(fres2, "%.15lf %.15lf\n",V,E );
  }
  fclose(fres2);


  // frontière 3/C
  phase=3;
  N=100;
  P=0;
  Tdeb=505;
  Tfin=2500;
  for(int i=0; i<N; i++){
    T=Tdeb + i*(Tfin-Tdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    printf("%.15lf %.15lf\n",V*1e6,E );
    fprintf(fres3, "%.15lf %.15lf\n",V,E );
  }
  T=2500;
  Pdeb=0.;
  Pfin=70e9;
  for(int i=0; i<N; i++){
    P=Pdeb + i*(Pfin-Pdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    printf("%.15lf %.15lf\n",V*1e6,E );
    fprintf(fres3, "%.15lf %.15lf\n",V,E );
  }
  fclose(fres3);
 

  // Liberation
  free(segP12);   free(segT12);
  free(segVa12);  free(segEa12);
  free(segVb12);  free(segEb12);

  free(segP13);   free(segT13);
  free(segVa13);  free(segEa13);
  free(segVc13);  free(segEc13);

  free(segP23);   free(segT23);
  free(segVb23);  free(segEb23);
  free(segVc23);  free(segEc23);

  free(segPHAUT);    free(segPBAS);
  free(segTGAUCHE);  free(segTDROITE);

}



//  ****  Calcul du Point TRIPLE  ****
void point_triple(void){
  int Nmax=1e4;
  double epsilon=1e-8;
  double res;

  double Pa, Ta, Sa, Ga, Ea, Va;
  double Pb, Tb, Sb, Gb, Eb, Vb;
  double Pc, Tc, Sc, Gc, Ec, Vc;

  double Va0=129e-6;
  double Ea0=100e3;

  double Vb0=0.98*Va0;
  double Eb0=1.01*Ea0;
  
  double Vc0=1.02*Va0;
  double Ec0=1.05*Ea0;
  
  
  triple(Nmax, epsilon, &Va0, &Ea0, &Vb0, &Eb0, &Vc0, &Ec0);
  

  // phase 1
  Ea=Ea0; Va=Va0;
  Pa=fP(Va,Ea, 1); Ta=fT(Va,Ea, 1);
  Sa=fS(Va,Ea, 1); Ga=Ea+Pa*Va-Ta*Sa;
  //phase 2
  Eb=Eb0; Vb=Vb0;
  Pb=fP(Vb,Eb, 2); Tb=fT(Vb,Eb, 2);
  Sb=fS(Vb,Eb, 2); Gb=Eb+Pb*Vb-Tb*Sb;
  // phase 3
  Ec=Ec0; Vc=Vc0;
  Pc=fP(Vc,Ec, 3); Tc=fT(Vc,Ec, 3);
  Sc=fS(Vc,Ec, 3); Gc=Ec+Pc*Vc-Tc*Sc;

  printf("Resultats :\n");
  printf(" Ga= %.15lf, Gb= %.15lf,  Gc= %.15lf\n",Ga, Gb, Gc );
  printf(" Pa= %.15lf, Pb= %.15lf,  Pc= %.15lf  (GPa)\n",Pa*1e-9, Pb*1e-9, Pc*1e-9 );
  printf(" Ta= %.15lf, Tb= %.15lf,  Tc= %.15lf\n",Ta, Tb, Tc );
}




// fonction pour determiner si un couple (V,E) est dans le triangle du poit triple
/*  -int *pind : pointeur sur l'indice (valeur retour)  
*/
int ftri(double V, double E){
  double epsilon=1e-6;

  // Points triple
  double VA= 129.905671875844291e-6;
  double EA= 72.168190127265447e3;
  double VB= 127.666112185063994e-6;
  double EB= 99.480256119236515e3;
  double VC= 132.538478805442338e-6;
  double EC= 136.154346228414312e3;

  // TRIANGLE
  double penteBA, penteBC, penteAC;

  penteBA=(EA-EB)/(VA-VB);
  penteBC=(EC-EB)/(VC-VB);
  penteAC=(EC-EA)/(VC-VA);

  double EBA=EB + penteBA*(V-VB);
  double EAC=EA + penteAC*(V-VA);
  double EBC=EB + penteBC*(V-VB);

  // triangle
  if(E>EBA && E<EBC && E>EAC){
    return 1;
  } 
  else{ return 0;}

}


// Calcul de VE
/*
   Discretise les frontière du diagramme de phase (P,T)
   1 fichier pour chaque phase : "VEfront1.txt"  ,  "VEfront2.txt"  et  "VEfront3.txt"
*/
int routine_VE(void){
  int err;
  int phase;
  double P;
  double T;
  double E,V;

  
  FILE* fres1;
  FILE* fres2;
  FILE* fres3;
  if((fres1 = fopen("VEfront1.txt", "w+")) == NULL){printf("erreur ouverture fichier (routine_VE)\n");return 1;}
  if((fres2 = fopen("VEfront2.txt", "w+")) == NULL){printf("erreur ouverture fichier (routine_VE)\n");return 1;}
  if((fres3 = fopen("VEfront3.txt", "w+")) == NULL){printf("erreur ouverture fichier (routine_VE)\n");return 1;}


  int N;
  double Tdeb, Tfin;
  double Pdeb, Pfin;

  // fronitère de la phase 1/A
  phase=1;
  N=100;
  P=0;
  Tdeb=505;
  Tfin=300;
  for(int i=0; i<N; i++){
    T=Tdeb + i*(Tfin-Tdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    printf("%.15lf %.15lf\n",V*1e6,E );
    fprintf(fres1, "%.15lf %.15lf\n",V,E );
  }
  T=300;
  Pdeb=0;
  Pfin=9.4e9;
  for(int i=0; i<N; i++){
    P=Pdeb + i*(Pfin-Pdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    printf("%.15lf %.15lf\n",V*1e6,E );
    fprintf(fres1, "%.15lf %.15lf\n",V,E );
  }
  fclose(fres1);

  //  ******************************************************

  // fronitère de la phase 2/B
  phase=2;
  N=100;
  T=300;
  Pdeb=9.4e9;
  Pfin=70e9;
  for(int i=0; i<N; i++){
    P=Pdeb + i*(Pfin-Pdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    printf("%.5lf %.5lf\n",V*1e6,E );
    fprintf(fres2, "%.15lf %.15lf\n",V,E );
  }
  P=70e9;
  Tdeb=300.;
  Tfin=2486.138962087622531;
  for(int i=0; i<N; i++){
    T=Tdeb + i*(Tfin-Tdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    printf("%.5lf %.5lf\n",V*1e6,E );
    fprintf(fres2, "%.15lf %.15lf\n",V,E );
  }
  fclose(fres2);


  //  ******************************************************

  // fronitère de la phase 3/C
  phase=3;
  N=100;
  P=0;
  Tdeb=505;
  Tfin=2500;
  for(int i=0; i<N; i++){
    T=Tdeb + i*(Tfin-Tdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    printf("%.15lf %.15lf\n",V*1e6,E );
    fprintf(fres3, "%.15lf %.15lf\n",V,E );
  }
  T=2500;
  Pdeb=0.;
  Pfin=70e9;
  for(int i=0; i<N; i++){
    P=Pdeb + i*(Pfin-Pdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    //printf("%.15lf %.15lf\n",V*1e6,E );
    fprintf(fres3, "%.15lf %.15lf\n",V,E );
  }
  fclose(fres3);

}


// fonction qui remplie le tableau S des somments  
/* Arguments :
  - int N13, N12.
  - double* P12, T12 : tableaux de N12 éléments qui contients les points de la courbe C12
                       partant du point triple
  - double* P13, T13 : tableaux de N13 éléments qui contients les points de la courbe C13
                       partant du point triple

  Arguments de sorties: void
  - double** S : les sommets du polygone tendant vers C 
*/
int sommetsA(int N13, int N12, int NTG, int NPB, double PrefG, double PrefD, double TrefH, double TrefB, double* segVa13, double* segEa13, double* segVa12, double* segEa12, double** S){
  
  // le remplisage de S doit se faire à l'étape avant

  int err=0;
  int phase=1;
  int N=N13+(NTG-1)+(NPB-1)+(N12-2);

  double P,T;
  double V,E;

  /*
  Trefbas=T12[N12-1];
  TrefG=T13[N13-1];
  PrefG=P13[N13-1];
  PrefD=P12[N12-1];
  */

  double hTG=(TrefB-TrefH)/(NTG-1);
  double hPB=(PrefD-PrefG)/(NPB-1);

  
  int nref=0;

  // **********   ECRITURE DE LA MATRICE DES SOMMETS   *************
  //  Parcours de C23 duu point triple jusqu'a droite 
  for(int j=0; j<=N13-1; j++){ // N13 points
    /*
    T=T13[j];
    P=P13[j];
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsA)\n");return err;}
    */
    V=segVa13[j];
    E=segEa13[j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }
  
  nref+=N13;

  // Parcours de la frontière droite
  P=PrefG;
  for(int j=1; j<=NTG-1; j++){ // N13 points
    T=TrefH+j*hTG;
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsA)\n");return err;}
    S[nref+j-1][0]=V;
    S[nref+j-1][1]=E;
  }

  nref+=NTG-1;
  
  // Parcours de la frontière haute
  T=TrefB;
  for(int j=1; j<=NPB-2; j++){ // N13 points
    P=PrefG+j*hPB;
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsA)\n");return err;}
    S[nref+j-1][0]=V;
    S[nref+j-1][1]=E;
  }

  nref+=NPB-2;

  // Parcours de C12 de T=Trefbas jusqu'au point triple 
  for(int j=0; j<=N12-2; j++){ // N13 points
    /*
    T=T12[N12-1-j];
    P=P12[N12-1-j];
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsA)\n");return err;}
    */
    V=segVa12[N12-1-j];
    E=segEa12[N12-1-j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }

  S[N][0]=S[0][0];
  S[N][1]=S[0][1];

  return 0;
  
}

// fonction qui remplie le tableau S des somments  
/* Arguments :
  - int N13, N23.
  - double* P12, T12 : tableaux de N12 éléments qui contients les points de la courbe C12
                       partant du point triple
  - double* P13, T13 : tableaux de N13 éléments qui contients les points de la courbe C13
                       partant du point triple

  Arguments de sorties: void
  - double** S : les sommets du polygone tendant vers C 
*/
int sommetsB(int N12, int N23, int NPB, int NTD, double PrefG, double PrefD, double TrefB, double TrefH, double* segVb12, double* segEb12, double* segVb23, double* segEb23, double** S){
  
  // le remplisage de S doit se faire à l'étape avant


  int N=N12+(NPB-1)+(NTD-2)+(N23-1);
  int err=0;
  int phase=2;

  double P,T;
  double V,E;
  
  /*
  TrefB=T12[N12-1];
  TrefH=T23[N23-1];
  PrefG=P12[N12-1];
  PrefD=P23[N23-1];
  */

  double hTD=(TrefH-TrefB)/(NTD-1);
  double hPB=(PrefD-PrefG)/(NPB-1);

  
  int nref=0;

  // **********   ECRITURE DE LA MATRICE DES SOMMETS   *************
  //  Parcours de C12 duu point triple jusqu'a droite 
  for(int j=0; j<=N12-1; j++){ // N12 points
    /*
    T=T12[j];
    P=P12[j];
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsB)\n");return err;}
    */
    V=segVb12[j];
    E=segEb12[j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }
  
  nref+=N12;

  // Parcours de la frontière basse
  T=TrefB;
  for(int j=1; j<=NPB-1; j++){ 
    P=PrefG+j*hPB;
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsB)\n");return err;}
    S[nref+j-1][0]=V;
    S[nref+j-1][1]=E;
  }

  nref+=NPB-1;
  
  // Parcours de la frontière droite
  P=PrefD;
  for(int j=1; j<=NTD-2; j++){ 
    T=TrefB+j*hTD;
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsB)\n");return err;}
    S[nref+j-1][0]=V;
    S[nref+j-1][1]=E;
  }

  nref+=NTD-2;


  //  Parcours de C23 de la droite jusqu'au point triple 
  for(int j=0; j<=N23-2; j++){ 
    /*
    T=T23[N23-1-j];
    P=P23[N23-1-j];
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsB)\n");return err;}
    */
    V=segVb23[N23-1-j];
    E=segEb23[N23-1-j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }

  S[N][0]=S[0][0];
  S[N][1]=S[0][1];

  return 0;
}


// fonction qui remplie le tableau S des somments  
/* Arguments :
  - int N13, N23.
  - double* P12, T12 : tableaux de N12 éléments qui contients les points de la courbe C12
                       partant du point triple
  - double* P13, T13 : tableaux de N13 éléments qui contients les points de la courbe C13
                       partant du point triple

  Arguments de sorties: void
  - double** S : les sommets du polygone tendant vers C 
*/
int sommetsC(int N13, int N23, int NTD, int NPH, int NTG, double PrefG, double PrefD, double TrefG, double TrefD, double Trefhaut, double* segVc13, double* segEc13, double* segVc23, double* segEc23 ,double** S){
  

  int N=N23+(NTD-1)+(NPH-1)+(NTG-1)+(N13-2) ;

  int err=0;
  int phase=3;

  double P,T;
  double V,E;
  
  /*
  Trefhaut=2500;
  TrefD=T23[N23-1];
  TrefG=T13[N13-1];
  PrefG=P13[N13-1];
  PrefD=P23[N23-1];
  */

  double hTD=(Trefhaut-TrefD)/(NTD-1);
  double hPH=(PrefG-PrefD)/(NPH-1);
  double hTG=(TrefG-Trefhaut)/(NTG-1);


  int nref=0;

  // **********   ECRITURE DE LA MATRICE DES SOMMETS   *************
  //  Parcours de C23 duu point triple jusqu'a droite 
  for(int j=0; j<=N23-1; j++){ // N13 points
    /*
    T=T23[j];
    P=P23[j];
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsC)\n");return err;}
    */
    V=segVc23[j];
    E=segEc23[j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }
  
  nref+=N23;

  // Parcours de la frontière droite
  P=PrefD;
  for(int j=1; j<=NTD-1; j++){ 
    T=TrefD+j*hTD;
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsC)\n");return err;}
    S[nref+j-1][0]=V;
    S[nref+j-1][1]=E;
  }

  nref+=NTD-1;
  
  // Parcours de la frontière haute
  T=Trefhaut;
  for(int j=1; j<=NPH-1; j++){ 
    P=PrefD+j*hPH;
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsC)\n");return err;}
    S[nref+j-1][0]=V;
    S[nref+j-1][1]=E;
  }

  nref+=NPH-1;

  // Parcours de la frontière gauche
  P=TrefG;
  for(int j=1; j<=NTG-2; j++){ // N13 points
    T=Trefhaut+j*hTG;
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsC)\n");return err;}
    S[nref+j-1][0]=V;
    S[nref+j-1][1]=E;
  }

  nref+=NTG-2;

  //  Parcours de C13 de la gauche jusqu'au point triple 
  for(int j=0; j<=N13-2; j++){ // N13 points
    /*
    T=T13[N13-1-j];
    P=P13[N13-1-j];
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsC)\n");return err;}
    */
    V=segVc13[N13-1-j];
    E=segEc13[N13-1-j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }

  S[N][0]=S[0][0];
  S[N][1]=S[0][1];
  
  return 0;
}


// fonction qui remplie le tableau S des somments  
/* Arguments :
  - int N13.
  - int ND, NG : nombre de points sur les fornitères droites et gauches
  - double* P12, T12 : tableaux de N12 éléments qui contients les points de la courbe C12
                       partant du point triple


  Arguments de sorties: void
  - double** S : les sommets du polygone tendant vers la zone de mélange 12
*/
int sommetsAB(int N12, double* segVa12, double* segEa12, double* segVb12, double* segEb12, double** S){
  
  int N=2*N12;

  int err=0;

  double V,E;
  

  int nref=0;

  // **********   ECRITURE DE LA MATRICE DES SOMMETS   *************
  //  Parcours de C23 duu point triple jusqu'a gauche (dans le plans VE)
  for(int j=0; j<=N12-1; j++){ 
    V=segVb12[j];
    E=segEb12[j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }
  
  nref+=N12;
  

  //  Parcours de C12 de la gauche jusqu'au point triple 
  for(int j=0; j<=N12-1; j++){ // N13 points
    V=segVa12[N12-1-j];
    E=segEa12[N12-1-j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }


  S[N][0]=S[0][0];
  S[N][1]=S[0][1];
  
  return 0;
}


// fonction qui remplie le tableau S des somments  
/* Arguments :
  - int N13.
  - int ND, NG : nombre de points sur les fornitères droites et gauches
  - double* P12, T12 : tableaux de N12 éléments qui contients les points de la courbe C12
                       partant du point triple


  Arguments de sorties: void
  - double** S : les sommets du polygone tendant vers la zone de mélange 12
*/
int sommetsAC(int N13, double* segVa13, double* segEa13, double* segVc13, double* segEc13, double** S){

  int N=2*N13 ;

  int err=0;

  double V,E;

  int nref=0;

  // **********   ECRITURE DE LA MATRICE DES SOMMETS   *************
  //  Parcours de C13 duu point triple jusqu'a gauche (dans le plans VE)
  for(int j=0; j<=N13-1; j++){ 
    V=segVc13[j];
    E=segEc13[j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }
  
  nref+=N13;
  

  //  Parcours de C13 de la gauche jusqu'au point triple 
  for(int j=0; j<=N13-1; j++){ // N13 points
    V=segVa13[N13-1-j];
    E=segEa13[N13-1-j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }


  S[N][0]=S[0][0];
  S[N][1]=S[0][1];
  
  return 0;
}


// fonction qui remplie le tableau S des somments  
/* Arguments :
  - int N13.
  - int ND, NG : nombre de points sur les fornitères droites et gauches
  - double* P12, T12 : tableaux de N12 éléments qui contients les points de la courbe C12
                       partant du point triple


  Arguments de sorties: void
  - double** S : les sommets du polygone tendant vers la zone de mélange 12
*/
int sommetsBC(int N23, double* segVb23, double* segEb23, double* segVc23, double* segEc23, double** S){
  
  int N=2*N23;

  int err=0;

  double V,E;

  int nref=0;

  // **********   ECRITURE DE LA MATRICE DES SOMMETS   *************
  //  Parcours de C23 duu point triple jusqu'a gauche (dans le plans VE)
  for(int j=0; j<=N23-1; j++){ 
    V=segVc23[j];
    E=segEc23[j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }
  
  nref+=N23;
  

  //  Parcours de C13 de la gauche jusqu'au point triple 
  for(int j=0; j<=N23-1; j++){ // N13 points
    V=segVb23[N23-1-j];
    E=segEb23[N23-1-j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }


  S[N][0]=S[0][0];
  S[N][1]=S[0][1];
  
  return 0;
}



/*
   Renvoie 1 si (x,y) appartient au poly def par S
   et 0 sinon
*/
int intPOLY(double x, double y, int N, double** S){

  double m, pente;
  int cpt=0;
  
  double xi, yi;
  double xi1, yi1;
  double xI, yI;
  double xext, yext;
  double Yext, Xext;

  Yext=S[0][1];
  for(int i=1; i<N; i++){
    if(Yext<S[i][1]){
      Yext=S[i][1];
    }
  }
  Yext+=abs(floor(Yext/2)); // on s'assure que (., Yext) est hors du poly S
  
  xext=x;
  yext=Yext;
  
  for(int i=0; i<N; i++){
    xi=S[i][0]; xi1=S[i+1][0];
    yi=S[i][1]; yi1=S[i+1][1];
    if(xi!=xi1){
      
      pente=(yi1-yi)/(xi1-xi);
      xI=x;
      yI=pente*(x-xi)+yi;
      if( (xI-xi)*(xI-xi1)<0 &&  (yI-y)*(yI-yext)<0){  // (yI-yi)*(yI-yi1)<0
        cpt++;
        //printf("  #######   cpt +1   ###########\n");
      }
      /*
      printf("   --i= %d\n",i );
      printf("        V= %lf  ,   E= %lf\n",x*1e6,y );
      printf("        xI= %lf  ,   yI= %lf\n",xI*1e6,yI );
      printf("        x_ext= %lf  ,   y_ext= %lf\n",xext*1e6,yext );
      printf("       (xI-xi)= %lf , (xI-xi1)= %lf\n",(xI-xi)*1e6, (xI-xi1)*1e6 );
      printf("       (yI-yi)= %lf , (yI-yi1)= %lf\n",(yI-yi)*1e6, (yI-yi1) );
      printf("       (yI-y)= %lf , (yI-yext)= %lf\n",(yI-y)*1e6, (yI-yext) );
      printf("       Vi= %lf Ei= %lf   ||   Vi1= %lf Ei1= %lf\n",xi*1e6,yi ,xi1*1e6,yi1);
      */
    }
  }
  //printf("  * cpt= %d\n",cpt );
  if(cpt%2==1){return 1;} // (x,y) est à l'intérieur  entier impair
  else{return 0;}         // (x,y) extérieur

}



/*  
    Creation des fichiers de maillage
*/
//int sommets_polygone(int NA, int NB, int NC, int NTG1, int NPB1, int NPB2, int NTD2, int NTD3, int NPH3, int NTG3, int N12, int N13, int N23, double ** SA, double** SB, double** SC, int NGAB, int NDAB, int NGAC, int NDAC, int NGBC, int NDBC, double** SAB, double** SAC, double** SBC){
int sommets_polygone(int NA, int NB, int NC, int NTG1, int NPB1, int NPB2, int NTD2, int NTD3, int NPH3, int NTG3, int N12, int N13, int N23, 
                     double ** SA, double** SB, double** SC, 
                     double** SAB, double** SAC, double** SBC,
                     double* segP12, double* segT12, double* segP13, double* segT13, double* segP23, double* segT23,
                     double* segVa12, double* segEa12, double* segVb12, double* segEb12,
                     double* segVa13, double* segEa13, double* segVc13, double* segEc13,
                     double* segVb23, double* segEb23, double* segVc23, double* segEc23){

  int err;
  int Nmax=1e4;
  double epsilon=5e-5;
  
  
  //printf("Sommets A\n");
  //int NA=N13+(NTG1-1)+(NPB1-1)+(N12-2);
  //printf(" -- NA= %d\n",NA );
  double TrefBA=segT12[N12-1];
  double TrefHA=segT13[N13-1];
  double PrefGA=segP13[N13-1];
  double PrefDA=segP12[N12-1];
  err=sommetsA(N13, N12, NTG1, NPB1, PrefGA, PrefDA, TrefHA, TrefBA, segVa13, segEa13, segVa12, segEa12, SA);
  if(err){return err;}
  
  //printf("Sommets B\n");
  //int NB=N12+(NPB2-1)+(NTD2-1)+(N23-2);
  //printf(" -- NB= %d\n",NB );
  double TrefBB=segT12[N12-1];
  double TrefHB=segT23[N23-1];
  double PrefGB=segP12[N12-1];
  double PrefDB=segP23[N23-1];
  err=sommetsB(N12, N23, NPB2, NTD2, PrefGB, PrefDB, TrefBB, TrefHB, segVb12, segEb12, segVb23, segEb23, SB);
  if(err){return err;}
  
  //printf("Sommets C\n");
  //int NC=N23+(NTD3-1)+(NPH3-1)+(NTG3-1)+(N13-2);
  //printf(" -- NC= %d\n",NC );
  double TrefhautC=2500;
  double TrefDC=segT23[N23-1];
  double TrefGC=segT13[N13-1];
  double PrefGC=segP13[N13-1];
  double PrefDC=segP23[N23-1];

  err=sommetsC(N13, N23, NTD3, NPH3, NTG3, PrefGC, PrefDC, TrefGC, TrefDC, TrefhautC, segVc13, segEc13, segVc23, segEc23, SC);
  if(err){return err;}
  
  //printf("Sommets AB\n");
  err=sommetsAB(N12, segVa12, segEa12, segVb12, segEb12, SAB);
  if(err){return err;}
  //printf("Sommets AC\n");
  err=sommetsAC(N13, segVa13, segEa13, segVc13, segEc13, SAC);
  if(err){return err;}
  //printf("Sommets BC\n");
  err=sommetsBC(N23, segVb23, segEb23, segVc23, segEc23, SBC);
  if(err){return err;}
  

  
  return 0;
}





// ZONE 
// repere les zones (triangle pour le moment)
int diag_zone(){
  int err,res;
  int Nv=3e2;
  int Ne=3e2; 
  int Nmax=1e4;
  double epsilon=1e-6;

  double E,V;


  int nb_points=10;
  int NTG1=nb_points, NPB1=2*nb_points;
  int NPB2=5*nb_points, NTD2=4*nb_points;
  int NTD3=ceil(nb_points/3), NPH3=6*nb_points, NTG3=4*nb_points; 
  int N12=nb_points; int N13=nb_points; int N23=5*nb_points;
  
  int NA=N13+(NTG1-1)+(NPB1-1)+(N12-2);
  double** SA=alloctabd(NA+1,2);
  
  int NB=N12+(NPB2-1)+(NTD2-1)+(N23-2);
  double** SB=alloctabd(NB+1,2);
  
  int NC=N23+(NTD3-1)+(NPH3-1)+(NTG3-1)+(N13-2);
  double** SC=alloctabd(NC+1,2); // on ajoute S[N]=S[0]  --> utile pou l'algo
  
  int NAB=2*N12;
  double** SAB=alloctabd(NAB+1,2);

  int NAC=2*N13;
  double** SAC=alloctabd(NAC+1,2);

  int NBC=2*N23;
  double** SBC=alloctabd(NBC+1,2);



  double Ptriple=4.5104217816091e9;
  double Pdeb=Ptriple;


  int phaseA12=1;
  int phaseB12=2;
  double Pfin12=9.4e9;
  
  int phaseA13=1;
  int phaseB13=3;
  double Pfin13=0.;
  
  int phaseA23=2;
  int phaseB23=3;
  double Pfin23=70e9;
  
  double Va012, Vb012, Va013, Vb013, Va023, Vb023;
  double Ea012, Eb012, Ea013, Eb013, Ea023, Eb023;

  /*
  VA= 129.905671875844291 (cm^3/kg)
  EA= 72.168190127265447 (kJ/kg)
  VB= 127.666112185063994 (cm^3/kg)
  EB= 99.480256119236515 (kJ/kg)
  VC= 132.538478805442338 (cm^3/kg)
  EC= 136.154346228414312 (kJ/kg)
  */
  
  // Valeur d'initialisation
  Va012=129.905671875844291e-6;  Ea012=72.168190127265447e3;
  Vb012=127.666112185063994e-6;  Eb012=99.480256119236515e3;
  
  Va013=129.905671875844291e-6;  Ea013=72.168190127265447e3;
  Vb013=132.538478805442338e-6;  Eb013=136.154346228414312e3;
  
  Va023=127.666112185063994e-6;  Ea023=99.480256119236515e3;
  Vb023=132.538478805442338e-6;  Eb023=136.154346228414312e3;
  
   
  // Initialisation des tableaux
  double* segP12=malloc(N12*sizeof(double));
  double* segT12=malloc(N12*sizeof(double));
  double* segVa12=malloc(N12*sizeof(double));
  double* segEa12=malloc(N12*sizeof(double));
  double* segVb12=malloc(N12*sizeof(double));
  double* segEb12=malloc(N12*sizeof(double));
  
  double* segP13=malloc(N13*sizeof(double));
  double* segT13=malloc(N13*sizeof(double));
  double* segVa13=malloc(N13*sizeof(double));
  double* segEa13=malloc(N13*sizeof(double));
  double* segVc13=malloc(N13*sizeof(double));
  double* segEc13=malloc(N13*sizeof(double));
  
  double* segP23=malloc(N23*sizeof(double));
  double* segT23=malloc(N23*sizeof(double));
  double* segVb23=malloc(N23*sizeof(double));
  double* segEb23=malloc(N23*sizeof(double));
  double* segVc23=malloc(N23*sizeof(double));
  double* segEc23=malloc(N23*sizeof(double));
  

  // frontiere 12
  //printf("PHASE 12 :\n");
  ligne_double(Nmax, epsilon, phaseA12, phaseB12, Pdeb, Pfin12, N12, Va012, Ea012, Vb012, Eb012, segP12, segT12, segVa12, segEa12, segVb12, segEb12);
  
  // frontiere 13
  //printf("PHASE 13 :\n");
  ligne_double(Nmax, epsilon, phaseA13, phaseB13, Pdeb, Pfin13, N13, Va013, Ea013, Vb013, Eb013, segP13, segT13, segVa13, segEa13, segVc13, segEc13);
  
  // frontiere 23
  //printf("PHASE 23 :\n");
  ligne_double(Nmax, epsilon, phaseA23, phaseB23, Pdeb, Pfin23, N23, Va023, Ea023, Vb023, Eb023, segP23, segT23, segVb23, segEb23, segVc23, segEc23);
  


  err=sommets_polygone(NA, NB, NC, NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3, N12, N13, N23, SA, SB, SC, SAB, SAC, SBC, segP12, segT12, segP13, segT13, segP23, segT23, segVa12, segEa12, segVb12, segEb12, segVa13, segEa13, segVc13, segEc13, segVb23, segEb23, segVc23, segEc23 );


  if(err){printf("err= %d (sommets_polygone)\n",err ); return err;}
  

  /*
  for(int i=0; i<NA+1; i++){
    printf("   SA[%d][0]= %lf  SA[%d][1]= %lf \n",i,SA[i][0],i,SA[i][1] );
  }
  printf("\n");
  for(int i=0; i<NB+1; i++){
    printf("   SB[%d][0]= %lf  SB[%d][1]= %lf \n",i,SB[i][0],i,SB[i][1] );
  }
  printf("\n");
  for(int i=0; i<NC+1; i++){
    printf("   SC[%d][0]= %lf  SC[%d][1]= %lf \n",i,SC[i][0],i,SC[i][1] );
  }
  */


  // Ecriture des sommets des poly pour vérif
  FILE *fSA, *fSB, *fSC, *fSAB, *fSAC, *fSBC;
  if((fSA = fopen("pSA.txt", "w+")) == NULL){printf("erreur ouverture fichier pSA\n");return 1;}
  if((fSB = fopen("pSB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSB\n");return 1;}
  if((fSC = fopen("pSC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSC\n");return 1;}
  if((fSAB = fopen("pSAB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAB\n");return 1;}
  if((fSAC = fopen("pSAC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAC\n");return 1;}
  if((fSBC = fopen("pSBC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSBC\n");return 1;}
  

  for(int i=0; i<NA+1; i++){
    fprintf(fSA,"%.10lf %.10lf\n",SA[i][0],SA[i][1] );
  }

  for(int i=0; i<NB+1; i++){
    fprintf(fSB," %.10lf %.10lf\n",SB[i][0],SB[i][1] );
  }

  for(int i=0; i<NC+1; i++){
    fprintf(fSC,"%.10lf %.10lf\n",SC[i][0],SC[i][1] );
  }

  for(int i=0; i<NAB+1; i++){
    fprintf(fSAB,"%.10lf %.10lf\n",SAB[i][0],SAB[i][1] );
  }

  for(int i=0; i<NAC+1; i++){
    fprintf(fSAC," %.10lf %.10lf\n",SAC[i][0],SAC[i][1] );
  }

  for(int i=0; i<NBC+1; i++){
    fprintf(fSBC,"%.10lf %.10lf\n",SBC[i][0],SBC[i][1] );
  }
  
  
  fclose(fSA);   fclose(fSB);   fclose(fSC);
  fclose(fSAB);  fclose(fSAC);  fclose(fSBC);
  
  
  double Vmin, Vmax;
  double Emin, Emax;
  /*
  Vmin=100e-6, Vmax=150e-6;
  Emin=-2e3, Emax=1e6;
  */
  Vmin=100e-6, Vmax=190e-6;
  Emin=-5e3, Emax=16e5;
  
  double he=(Emax-Emin)/(Ne-1);
  double hv=(Vmax-Vmin)/(Nv-1);
  
  
  FILE* fZONE; FILE* fE; FILE* fV;
  if((fZONE = fopen("ZONE.txt", "w+")) == NULL){printf("erreur ouverture fichier ZONE\n");return 1;}
  if((fV = fopen("Vzone.txt", "w+")) == NULL){printf("erreur ouverture fichier Vzone \n");return 1;}
  if((fE = fopen("Ezone.txt", "w+")) == NULL){printf("erreur ouverture fichier Ezone \n");return 1;}
  
  for(int j=0; j<Ne; j++){
      E=Emin+j*he;
      fprintf(fE, "%.15lf ",E);
  }
  
  
  for(int i=0; i<Nv; i++){
    V=Vmin+i*hv;
    fprintf(fV, "%.15lf ", V);
    for(int j=0; j<Ne; j++){
      //printf("  -j= %d\n",j );
      E=Emin+j*he;
      //printf("V= %.15lf, E=%.15lf\n",V*1e6,E );
      res=ftri(V, E);
      if(res){fprintf(fZONE, "0 ");}
      else{
        res=intPOLY(V,E,NA,SA);
        if(res){fprintf(fZONE, "1 ");}
        else{
          res=intPOLY(V,E,NB,SB);
          if(res){fprintf(fZONE, "2 ");}
          else{
            res=intPOLY(V,E,NC,SC);
            if(res){fprintf(fZONE, "3 ");}
            else{
              res=intPOLY(V,E,NAB,SAB);
              if(res){fprintf(fZONE, "12 ");}
              else{
                res=intPOLY(V,E,NAC,SAC);
                if(res){fprintf(fZONE, "13 ");}
                else{
                  res=intPOLY(V,E,NBC,SBC);
                  if(res){fprintf(fZONE, "23 ");}
                  else{fprintf(fZONE, "-1 ");}
                }
              }
            }
          }
        }
      }
    }
    fprintf(fZONE, "\n");
  }
  

  fclose(fZONE);
  fclose(fV);
  fclose(fE);

    // Liberation
  free(segP12);  free(segT12);
  free(segVa12);  free(segEa12);
  free(segVb12);  free(segEb12);
  
  free(segP13);  free(segT13);
  free(segVa13);  free(segEa13);
  free(segVc13);  free(segEc13);
  
  free(segP23);  free(segT23);
  free(segVb23);  free(segEb23);
  free(segVc23);  free(segEc23);
  

}


/*  Valeur d'entree :
   Va,Ea, Vb,Eb, x

   corespodances :  1-a   2-b   3-c

   Valeur sortie
  Vx,Ex, Vy,Ey, 

   double* A et InvA permet une seule allocation 

*/
void Newton_mixte(double* A, double* InvA, int phaseX, int phaseY, double Vmel, double Emel, double Vx, double Ex, double Vy, double Ey, double x, double* pVX, double* pEX, double* pVY, double* pEY, double* pX, double* dX){
  int dim=5;

  // Grandeurs thermodynamiques
  double Px, Tx, Sx, Gx, dPvx, dPex, dTvx, dTex;
  double Py, Ty, Sy, Gy, dPvy, dPey, dTvy, dTey;
  // phase X
  Px=fP(Vx,Ex, phaseX); Tx=fT(Vx,Ex, phaseX);
  Sx=fS(Vx,Ex, phaseX); Gx=Ex+Px*Vx-Tx*Sx;
  dPvx=dPdV(Vx,phaseX); dPex=dPdE(phaseX);
  dTvx=dTdV(Vx,phaseX); dTex=dTdE(phaseX);
  // phase Y
  Py=fP(Vy,Ey, phaseY); Ty=fT(Vy,Ey, phaseY);
  Sy=fS(Vy,Ey, phaseY); Gy=Ey+Py*Vy-Ty*Sy;
  dPvy=dPdV(Vy,phaseY); dPey=dPdE(phaseY);
  dTvy=dTdV(Vy,phaseY); dTey=dTdE(phaseY);
  
  /*
  printf("Px= %.10lf, Tx= %.10lf, Sx= %.10lf, Gx= %.10lf\n",Px,Tx,Sx,Gx );
  printf("Py= %.10lf, Ty= %.10lf, Sy= %.10lf, Gy= %.10lf\n",Py,Ty,Sy,Gy );
  */
  
  double dGvx, dGex, dGvy, dGey;
  dGvx=Vx*dPvx - Sx*dTvx; 
  dGex=Vx*dPex - Sx*dTex;
  
  dGvy=Vy*dPvy - Sy*dTvy; 
  dGey=Vy*dPey - Sy*dTey;
  
  
  double* Delta=malloc(dim*sizeof(double));
  Delta[0]=Py-Px;   Delta[1]=Ty-Tx;   Delta[2]=Gy-Gx;
  Delta[3]=Vmel - (1.-x)*Vx - x*Vy;
  Delta[4]=Emel - (1.-x)*Ex - x*Ey;
   
  /*
  for(int j=0; j<dim; j++){
    printf("Delta[%d]= %.15lf  ",j, Delta[j] );
  }
  printf("\n");
  */
  
  
  // on met le second menbre dans la dernière colonne
  //matrice A
  A[dim*0+0]=dPvx; A[dim*0+1]=dPex; A[dim*0+2]=-dPvy; A[dim*0+3]=-dPey; A[dim*0+4]=0.;  
  
  A[dim*1+0]=dTvx; A[dim*1+1]=dTex; A[dim*1+2]=-dTvy; A[dim*1+3]=-dTey; A[dim*1+4]=0.; 
  
  A[dim*2+0]=dGvx; A[dim*2+1]=dGex; A[dim*2+2]=-dGvy; A[dim*2+3]=-dGey; A[dim*2+4]=0.; 
  
  A[dim*3+0]=1.-x;  A[dim*3+1]=0.;   A[dim*3+2]=x;     A[dim*3+3]=0;     A[dim*3+4]=Vy-Vx; 
  
  A[dim*4+0]=0.;   A[dim*4+1]=1.-x;  A[dim*4+2]=0;     A[dim*4+3]=x;     A[dim*4+4]=Ey-Ex; 
  
  /*
  printf("  * A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      printf("A[%d][%d]= %.15lf  ",i,j, A[dim*i+j] );
    }
    printf("\n");
  }
  printf("\n");  
  */


  double det_a=inverse_matrice_pivot(A, dim, InvA);
  //printf("detA= %.15lf\n",det_a );
  if ( fabs(det_a)<1e-8){printf("A non inversible : detA= %.15lf\n", det_a);}
  
  /*
  printf("  * InvA = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      printf("InvA[%d][%d]= %.15lf  ",i,j, InvA[dim*i+j] );
    }
    printf("\n");
  }
  printf("\n");
  */
  
  
  // Mutiplcation matricielle de A et InvA
  /*
  double res;
  printf("  * InvA*A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      res=0;
      for(int k=0; k<dim; k++){
        res+=InvA[dim*i+k]*A[dim*k+j];
      }
      printf("A*InvA[%d][%d]= %.15lf ",i,j,res );
    }
    printf("\n");
  }
  printf("\n");
  */


  // Mutiplcation matricielle
  for(int i=0; i<dim; i++){
    dX[i]=0;
    for(int j=0; j<dim; j++){
      dX[i]+=InvA[dim*i+j]*Delta[j];
    }
    //printf("  dX[%d]= %.15lf\n",i,dX[i] );
  }
  //printf("\n");
   
  double correc=1.;

  *pVX=Vx+correc*dX[0];
  *pEX=Ex+correc*dX[1];
  *pVY=Vy+correc*dX[2];
  *pEY=Ey+correc*dX[3];
  *pX=x+correc*dX[4];

  /*
  printf("  Vx= %.10lf, Ex= %.10lf   Vy= %.10lf, Ey= %.10lf    x= %.10lf\n",Vx,Ex,Vy,Ey,x );
  printf("  (Vmel-Vx)/(Vy-Vx)= %.10lf  ,    (Emel-Ex)/(Ey-Ex)= %.10lf\n", (Vmel-Vx)/(Vy-Vx), (Emel-Ex)/(Ey-Ex));
  printf("  (1-x)*Vx + x*Vy - Vmel= %.10lf  ,    (1-x)*Ex + x*Ey - Emel= %.10lf\n", (1-x)*Vx + x*Vy - Vmel, (1-x)*Ex + x*Ey - Emel);
  */

  free(Delta);
}


/*
Arguments d'entree :
  - Nmax : nombre max d'iteration
  - epsilon : tolérance sur le critère  (ici fabs(G1-G2))
  - Vmel, Emel : Point dont on sait qu'il est dans la zone de mélange (alpha/beta) (X/Y)
  - phaseX, phaseY : les deux phases dans la zone de mélange

Arguments de sortie :
  - Vx, Ex : valeur de la phase X
  - Vy, Ey : valeur de la phase Y
  -   x    : frcation massique de Y
*/
void VE_mixte(int Nmax, double epsilon, int phaseX, int phaseY, double V0x, double E0x, double V0y, double E0y, double Vmel, double Emel, double* pVx, double* pEx, double* pVy, double* pEy, double* px){
  int n=1;
  int dim=5;

  double* A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));
  double* dX=malloc(dim*sizeof(double));
  
  double Px, Tx, Sx, Gx;
  double Py, Ty, Sy, Gy;

  double Vx, Ex, Vy, Ey, x;
  double VX, EX, VY, EY, X;
  

  // POINT TRIPLE
  /*
  VA= 129.905671875844291 (cm^3/kg)
  EA= 72.168190127265447 (kJ/kg)
  VB= 127.666112185063994 (cm^3/kg)
  EB= 99.480256119236515 (kJ/kg)
  VC= 132.538478805442338 (cm^3/kg)
  EC= 136.154346228414312 (kJ/kg)
  */

  double* V0=malloc(3*sizeof(double));
  double* E0=malloc(3*sizeof(double));

  V0[0]=129.905671875844291*1e-6;
  V0[1]=127.666112185063994*1e-6;
  V0[2]=132.538478805442338*1e-6;

  E0[0]=72.168190127265447*1e3;
  E0[1]=99.480256119236515*1e3;
  E0[2]=136.154346228414312*1e3;
  
  //Initialisation
  /*
  Vx=V0[phaseX-1];  Ex=E0[phaseX-1];
  Vy=V0[phaseY-1];  Ey=E0[phaseY-1];
  */
  /*
  Vx=0.99995*Vmel;  Ex=0.99995*Emel;
  Vy=1.000001*Vmel; Ey=1.0000001*Emel;
  */

  Vx=V0x; Ex=E0x; 
  Vy=V0y; Ey=E0y; 


  x=1./2;


  // phase X
  Px=fP(Vx,Ex, phaseX); Tx=fT(Vx,Ex, phaseX);
  Sx=fS(Vx,Ex, phaseX); Gx=Ex+Px*Vx-Tx*Sx;
  // phase Y
  Py=fP(Vy,Ey, phaseY); Ty=fT(Vy,Ey, phaseY);
  Sy=fS(Vy,Ey, phaseY); Gy=Ey+Py*Vy-Ty*Sy;


  double critere;
  
  critere=1.;
  //printf("  critere= %.14lf\n",critere );
    
  n=1;
  while(critere>=epsilon && n<Nmax){
    //printf("  --n= %d\n",n );
    
    Newton_mixte(A, InvA, phaseX, phaseY, Vmel, Emel, Vx, Ex, Vy, Ey, x, &VX, &EX, &VY, &EY, &X, dX);
    
    /*
    printf(" dVa= %.15lf, VA= %.15lf (cm^3/kg)\n",(VA-Va)*1e6, VA*1e6 );
    printf(" dEa= %.15lf, EA= %.15lf (kJ/kg)\n",(EA-Ea)*1e-3, EA*1e-3 );
    printf(" dVb= %.15lf, VB= %.15lf (cm^3/kg)\n",(VB-Vb)*1e6, VB*1e6 );
    printf(" dEb= %.15lf, EB= %.15lf (kJ/kg)\n",(EB-Eb)*1e-3, EB*1e-3 );
    printf(" dVc= %.15lf, VC= %.15lf (cm^3/kg)\n",(VC-Vc)*1e6, VC*1e6 );
    printf(" dEc= %.15lf, EC= %.15lf (kJ/kg)\n",(EC-Ec)*1e-3, EC *1e-3);
    */

    Vx=VX; Ex=EX;
    Vy=VY; Ey=EY;
    x=X;
    
    // phase X
    Px=fP(Vx,Ex, phaseX); Tx=fT(Vx,Ex, phaseX);
    Sx=fS(Vx,Ex, phaseX); Gx=Ex+Px*Vx-Tx*Sx;
    // phase Y
    Py=fP(Vy,Ey, phaseY); Ty=fT(Vy,Ey, phaseY);
    Sy=fS(Vy,Ey, phaseY); Gy=Ey+Py*Vy-Ty*Sy;

    //critere=fabs(Gx-Gy)+fabs(Px-Py)+fabs(Ty-Tx);
    critere=fabs(dX[0])+fabs(dX[1])+fabs(dX[2])+fabs(dX[3])+fabs(dX[4]);
    //printf("    critere= %.14lf\n",critere );

    
    n++;
  }
  printf(" --n= %d, critere= %lf *1e-5\n",n,critere*1e5 );
  
  *pVx=Vx; *pEx=Ex;
  *pVy=Vy; *pEy=Ey;
  *px=x;

  //printf("  ok fin boucle i\n");
  free(A); 
  free(InvA);
  free(V0); free(E0);
  free(dX);

}



// Aire du parallélépipède porté par v1 et v2 =======================================
double aire(double* v1, double* v2){
  return v1[0]*v2[1] - v1[1]*v2[0];
}



//  coorbary  =======================================================================
/*
Renvoie les 3 coordonées barycentrique de (x,y) pour un triangle du maillage
Arguments d'entree :
 - double *xy        : coordonnées d'un point de la discretisation
 - double **sommets  : coordonnées des trois points consistuant un triangle du maillage

Arguments de sortie :
 - double **pcoorbar : pointeur sur coorbar, tableau contenant les trois coordonnées barycentrique du pts
*/
int coorbary(double *xy, double **sommets, double **pcoorbar){
  double epsilon=1e-7;
  double ajax[2], akai[2], ajai[2], akaj[2];
  double *coorbar=*pcoorbar;
  double a1a2a3;
  // (i,j,k)=(1,2,3) à une permutation circulaire pres
  // ai, aj, ak les trois sommets du triangle
  // ajx= vecteur(ak,xy) aji=vecteur(aj,ai)
  int i,j,k;
  for(i=0; i<2; i++){
    j=(i+1)%3; // permutation circulaire de (i,j,k)
    k=(i+2)%3;
    for(int d=0; d<2; d++){
      ajax[d]=xy[d]-sommets[j][d];
      akaj[d]=sommets[j][d]-sommets[k][d];
      akai[d]=sommets[i][d]-sommets[k][d];
      ajai[d]=sommets[i][d]-sommets[j][d];
    }
    a1a2a3=aire(ajai,akai);
    if(fabs(a1a2a3)<epsilon){ return 1;}
    coorbar[i]=aire(ajax,akaj)/a1a2a3;
  }
  coorbar[2]=1-coorbar[0]-coorbar[1];
  return 0;
}


/*
   Fonction qui renvoie les deux indices(corespondant à des droites) entre lesquelles se trouve le couple (V,E)
    --> le but étant de trouver la pression et la température en ce point
*/
int interp_zone_mixte(double V, double E, int phaseALPHA, int phaseBETA, int N, double* Valpha, double* Vbeta, double* Ealpha, double* Ebeta, double* pP, double* pT, double* pVmx, double* pEmx, double* pVmy, double* pEmy){
    double P,T;
    int res;
    double epsilon=1e-5;

    // reparage du quadrilatere
    double** S=alloctabd(5,2);
    int ind=0;
    int testNO=1;
    for(int i=0; i<=N-2; i++){
      S[0][0]=Valpha[i];   S[0][1]=Ealpha[i];
      S[1][0]=Valpha[i+1]; S[1][1]=Ealpha[i+1];
      S[2][0]=Vbeta[i+1];  S[2][1]=Ebeta[i+1];
      S[3][0]=Vbeta[i];    S[3][1]=Ebeta[i];
      S[4][0]=Valpha[i];   S[4][1]=Ealpha[i];

      if(intPOLY(V, E, 4, S)){ind=i; testNO=0;break;} // ou i=N
    }

    if(testNO){
      printf("\n");
      printf("le point est dans aucun des quadrilatere\n");
      printf(" ind= %d\n",ind );
      printf("V= %.10lf , E= %.10lf\n",V*1e6,E); 
      return 1;
    }

    *pVmx=Valpha[ind];  *pEmx=Ealpha[ind];
    *pVmy=Vbeta[ind];   *pEmy=Ebeta[ind];

    
    double pente1, pente2;

    // Calcul du poit milieu
    // diag de [Malpha_i, Mbeta_i+1]
    pente1=(Ebeta[ind+1]-Ealpha[ind])/(Vbeta[ind+1]-Valpha[ind]);
    //Ediag1=Ealpha[i] + pente1*(V-Valpha[i]);

    // diag de [Malpha_i+1, Mbeta_i]
    pente2=(Ebeta[ind]-Ealpha[ind+1])/(Vbeta[ind]-Valpha[ind+1]);
    //Ediag2=Ealpha[i+1] + pente2*(V-Valpha[i+1]);
    
    // On cherche V* tq     Ealpha[i+1] + pente2*(V*-Valpha[i+1]) = Ealpha[i] + pente1*(V*-Valpha[i])
    double VI=(Ealpha[ind]-Ealpha[ind+1] + pente2*Valpha[ind+1] - pente1*Valpha[ind])/(pente2-pente1);
    double EI=Ealpha[ind+1] + pente2*(VI-Valpha[ind+1]);

    //printf("  VI= %.10lf,  EI= %f\n",VI,EI );

    // Coordonées d'un triangle
    double** coorTri=alloctabd(3,2);

    // Calcul des 3 coordonnées barycentriques pour les 4 triangles
    double* coorBar=malloc(3*sizeof(double));  // une ligne = un triangle

    double* VE=malloc(2*sizeof(double));  // une ligne = un triangle
    VE[0]=V; VE[1]=E;


    // Calcul des pressions et températures aux sommets des triangles
    double* Ptri=malloc(3*sizeof(double));  // une ligne = un triangle
    double* Ttri=malloc(3*sizeof(double));  // une ligne = un triangle
    
    int test0, test1, test2;    

     
    //printf("  $$ 1er triangle  $$\n");
    // 1er triangle
    coorTri[0][0]=Valpha[ind];    coorTri[0][1]=Ealpha[ind];
    coorTri[1][0]=Valpha[ind+1];  coorTri[1][1]=Ealpha[ind+1];
    coorTri[2][0]=Vbeta[ind+1];   coorTri[2][1]=Ebeta[ind+1];

    Ptri[0]=fP(coorTri[0][0], coorTri[0][1], phaseALPHA);  
    Ptri[1]=fP(coorTri[1][0], coorTri[1][1], phaseALPHA);  
    Ptri[2]=fP(coorTri[2][0], coorTri[2][1], phaseBETA);  

    Ttri[0]=fT(coorTri[0][0], coorTri[0][1], phaseALPHA);  
    Ttri[1]=fT(coorTri[1][0], coorTri[1][1], phaseALPHA);  
    Ttri[2]=fT(coorTri[2][0], coorTri[2][1], phaseBETA);  


    res=coorbary(VE, coorTri, &coorBar);
    
    /*
    printf(" Triangle 1 :\n");
    for(int i=0; i<3; i++){
      printf(" coorBar[%d]= %lf",i,coorBar[i] );
    }
    printf("\n");
    */

    test0=( 0-epsilon<=coorBar[0] && coorBar[0]<=1+epsilon );
    test1=( 0-epsilon<=coorBar[1] && coorBar[1]<=1+epsilon );
    test2=( 0-epsilon<=coorBar[2] && coorBar[2]<=1+epsilon );
    if(test0 && test1 && test2){  // on est dans ce triangle
      P=0; T=0;
      for(int i=0; i<3; i++){
        P+=coorBar[i]*Ptri[i];  
        T+=coorBar[i]*Ttri[i];
      }
      *pP=P;
      *pT=T;
      return 0;
    }
    else{

      //printf("  $$ 2nd triangle  $$\n");
      // 2nd triangle
      coorTri[0][0]=Valpha[ind];    coorTri[0][1]=Ealpha[ind];
      coorTri[1][0]=Vbeta[ind];     coorTri[1][1]=Ebeta[ind];
      coorTri[2][0]=Vbeta[ind+1];   coorTri[2][1]=Ebeta[ind+1];

      Ptri[0]=fP(coorTri[0][0], coorTri[0][1], phaseALPHA);  
      Ptri[1]=fP(coorTri[1][0], coorTri[1][1], phaseBETA);  
      Ptri[2]=fP(coorTri[2][0], coorTri[2][1], phaseBETA);  

      Ttri[0]=fT(coorTri[0][0], coorTri[0][1], phaseALPHA);  
      Ttri[1]=fT(coorTri[1][0], coorTri[1][1], phaseBETA);  
      Ttri[2]=fT(coorTri[2][0], coorTri[2][1], phaseBETA);  

     
      res=coorbary(VE, coorTri, &coorBar);
      
      /*
      printf(" Triangle 2 :\n");
      for(int i=0; i<3; i++){
        printf(" coorBar[%d]= %lf",i,coorBar[i] );
      }
      printf("\n");
      */

      test0=( 0-epsilon<=coorBar[0] && coorBar[0]<=1+epsilon );
      test1=( 0-epsilon<=coorBar[1] && coorBar[1]<=1+epsilon );
      test2=( 0-epsilon<=coorBar[2] && coorBar[2]<=1+epsilon );
      if(test0 && test1 && test2){  // on est dans ce triangle
        P=0; T=0;
        for(int i=0; i<3; i++){
          P+=coorBar[i]*Ptri[i];  
          T+=coorBar[i]*Ttri[i];  
        }
        *pP=P;
        *pT=T;
        return 0;
      }
      else{
        printf("\nLe triangle n'a pas été trouvé (interp_zone_mixte)\n");
        
        printf("  ind= %d       N-1= %d\n",ind,N-1 );
        printf("V= %.10lf, E= %lf\n",V,E );
        printf("  4 sommets du quadrilatere :\n");
        printf("Valpha[%d]= %.10lf , Ealpha[%d]= %lf\n", ind, Valpha[ind], ind, Ealpha[ind]);
        printf("Valpha[%d]= %.10lf , Ealpha[%d]= %lf\n", ind+1, Valpha[ind+1], ind+1, Ealpha[ind+1]);
        printf("Vbeta[%d]= %.10lf , Ebeta[%d]= %lf\n", ind+1, Vbeta[ind+1], ind+1, Ebeta[ind+1]);
        printf("Vbeta[%d]= %.10lf , Ebeta[%d]= %lf\n", ind, Vbeta[ind], ind, Ebeta[ind]);
    
        for(int i=0; i<3; i++){
          printf("    coorBar[%d]= %lf\n",i,coorBar[i] );
        }
        return 1;  
      }
    }

  free(S);
  free(Ptri); free(Ttri);
  freetab(coorTri); 
  free(coorBar);
  free(VE);

}



/*
   Calcule dans quel quadrilatere d'une zone mixte donné se trouve le couple (V,E)
*/
int graine_zone_mixte(double V, double E, int N, double* Valpha, double* Vbeta, double* Ealpha, double* Ebeta, double* pVmx, double* pEmx, double* pVmy, double* pEmy){
    int res;

    // reparage du quadrilatere
    double** S=alloctabd(5,2);
    int ind=0;
    int testNO=1;
    for(int i=0; i<=N-2; i++){
      S[0][0]=Valpha[i];   S[0][1]=Ealpha[i];
      S[1][0]=Valpha[i+1]; S[1][1]=Ealpha[i+1];
      S[2][0]=Vbeta[i+1];  S[2][1]=Ebeta[i+1];
      S[3][0]=Vbeta[i];    S[3][1]=Ebeta[i];
      S[4][0]=Valpha[i];   S[4][1]=Ealpha[i];

      if(intPOLY(V, E, 4, S)){ind=i; testNO=0;break;} // ou i=N
    }

    if(testNO){
      printf("\n");
      printf("le point est dans aucun des quadrilatere (graine_zone_mixte)\n");
      printf(" ind= %d\n",ind );
      printf("V= %.10lf , E= %.10lf\n",V*1e6,E); 
      return 1;
    }

    *pVmx=Valpha[ind];  *pEmx=Ealpha[ind];
    *pVmy=Vbeta[ind];   *pEmy=Ebeta[ind];

    


  free(S);

}



/* Fonction pour mapper la (P,T) sur (V,E)

*/
int  mapsPT(void){
  
  int Nmax=1e4;
  double epsilon=1e-6;
  double Vx, Ex, Vy, Ey, x;
  double P,T;
  double PT,TT;

  double VAt= 129.905671875844291e-6;
  double EAt= 72.168190127265447e3;
  PT=fP(VAt, EAt,1);
  TT=fT(VAt, EAt,1);
  

  int err,res;
  int Nv=1e2;
  int Ne=1e2;

  double Vmin, Vmax;
  double Emin, Emax;
  
  Vmin=80e-6, Vmax=210e-6;
  Emin=-2e3, Emax=16e5;
  
  /*
  Vmin=110e-6, Vmax=150e-6;
  Emin=-5e3, Emax=5e5; 
  */

  double E,V;


  int nb_points=20;
  int NTG1=nb_points, NPB1=2*nb_points;
  int NPB2=5*nb_points, NTD2=4*nb_points;
  int NTD3=ceil(nb_points/3), NPH3=6*nb_points, NTG3=4*nb_points; 
  int N12=nb_points; int N13=nb_points; int N23=5*nb_points;
  
  int NA=N13+(NTG1-1)+(NPB1-1)+(N12-2);
  double** SA=alloctabd(NA+1,2);
  
  int NB=N12+(NPB2-1)+(NTD2-1)+(N23-2);
  double** SB=alloctabd(NB+1,2);
  
  int NC=N23+(NTD3-1)+(NPH3-1)+(NTG3-1)+(N13-2);
  double** SC=alloctabd(NC+1,2); // on ajoute S[N]=S[0]  --> utile pour l'algo
  
  int NAB=2*N12;
  double** SAB=alloctabd(NAB+1,2);

  int NAC=2*N13;
  double** SAC=alloctabd(NAC+1,2);

  int NBC=2*N23;
  double** SBC=alloctabd(NBC+1,2);


  double Ptriple=4.5104217816091e9;
  double Pdeb=Ptriple;


  int phaseA12=1;
  int phaseB12=2;
  double Pfin12=9.4e9;
  
  int phaseA13=1;
  int phaseB13=3;
  double Pfin13=0.;
  
  int phaseA23=2;
  int phaseB23=3;
  double Pfin23=70e9;
  
  double Va012, Vb012, Va013, Vb013, Va023, Vb023;
  double Ea012, Eb012, Ea013, Eb013, Ea023, Eb023;

  
  // Valeur d'initialisation
  Va012=129.905671875844291e-6;  Ea012=72.168190127265447e3;
  Vb012=127.666112185063994e-6;  Eb012=99.480256119236515e3;
  
  Va013=129.905671875844291e-6;  Ea013=72.168190127265447e3;
  Vb013=132.538478805442338e-6;  Eb013=136.154346228414312e3;
  
  Va023=127.666112185063994e-6;  Ea023=99.480256119236515e3;
  Vb023=132.538478805442338e-6;  Eb023=136.154346228414312e3;
  
   
  // Initialisation des tableaux
  double* segP12=malloc(N12*sizeof(double));
  double* segT12=malloc(N12*sizeof(double));
  double* segVa12=malloc(N12*sizeof(double));
  double* segEa12=malloc(N12*sizeof(double));
  double* segVb12=malloc(N12*sizeof(double));
  double* segEb12=malloc(N12*sizeof(double));
  
  double* segP13=malloc(N13*sizeof(double));
  double* segT13=malloc(N13*sizeof(double));
  double* segVa13=malloc(N13*sizeof(double));
  double* segEa13=malloc(N13*sizeof(double));
  double* segVc13=malloc(N13*sizeof(double));
  double* segEc13=malloc(N13*sizeof(double));
  
  double* segP23=malloc(N23*sizeof(double));
  double* segT23=malloc(N23*sizeof(double));
  double* segVb23=malloc(N23*sizeof(double));
  double* segEb23=malloc(N23*sizeof(double));
  double* segVc23=malloc(N23*sizeof(double));
  double* segEc23=malloc(N23*sizeof(double));
  
  // frontiere 12
  ligne_double(Nmax, epsilon, phaseA12, phaseB12, Pdeb, Pfin12, N12, Va012, Ea012, Vb012, Eb012, segP12, segT12, segVa12, segEa12, segVb12, segEb12);
  
  // frontiere 13
  ligne_double(Nmax, epsilon, phaseA13, phaseB13, Pdeb, Pfin13, N13, Va013, Ea013, Vb013, Eb013, segP13, segT13, segVa13, segEa13, segVc13, segEc13);
  
  // frontiere 23
  ligne_double(Nmax, epsilon, phaseA23, phaseB23, Pdeb, Pfin23, N23, Va023, Ea023, Vb023, Eb023, segP23, segT23, segVb23, segEb23, segVc23, segEc23);
  

  // CALCUL DES SOMMETS  
  err=sommets_polygone(NA, NB, NC, NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3, N12, N13, N23, SA, SB, SC, SAB, SAC, SBC, segP12, segT12, segP13, segT13, segP23, segT23, 
                       segVa12, segEa12, segVb12, segEb12, segVa13, segEa13, segVc13, segEc13, segVb23, segEb23, segVc23, segEc23 );

  if(err){printf("err= %d (sommets_polygone)\n",err ); return err;}



  ///////////
  // Ecriture des sommets des poly pour vérif
  FILE *fSA, *fSB, *fSC, *fSAB, *fSAC, *fSBC;
  if((fSA = fopen("pSA.txt", "w+")) == NULL){printf("erreur ouverture fichier pSA\n");return 1;}
  if((fSB = fopen("pSB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSB\n");return 1;}
  if((fSC = fopen("pSC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSC\n");return 1;}
  if((fSAB = fopen("pSAB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAB\n");return 1;}
  if((fSAC = fopen("pSAC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAC\n");return 1;}
  if((fSBC = fopen("pSBC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSBC\n");return 1;}
  
  for(int i=0; i<NA+1; i++){
    fprintf(fSA,"%.10lf %.10lf\n",SA[i][0],SA[i][1] );
  }

  for(int i=0; i<NB+1; i++){
    fprintf(fSB," %.10lf %.10lf\n",SB[i][0],SB[i][1] );
  }

  for(int i=0; i<NC+1; i++){
    fprintf(fSC,"%.10lf %.10lf\n",SC[i][0],SC[i][1] );
  }

  for(int i=0; i<NAB+1; i++){
    fprintf(fSAB,"%.10lf %.10lf\n",SAB[i][0],SAB[i][1] );
  }

  for(int i=0; i<NAC+1; i++){
    fprintf(fSAC," %.10lf %.10lf\n",SAC[i][0],SAC[i][1] );
  }

  for(int i=0; i<NBC+1; i++){
    fprintf(fSBC,"%.10lf %.10lf\n",SBC[i][0],SBC[i][1] );
  }
  fclose(fSA);   fclose(fSB);   fclose(fSC);
  fclose(fSAB);  fclose(fSAC);  fclose(fSBC);
 
  /////////////////

  double he=(Emax-Emin)/(Ne-1);
  double hv=(Vmax-Vmin)/(Nv-1);


  FILE *fmaps, *fmapsGNU, *fE, *fV;
  if((fmaps = fopen("MAPS.txt", "w+")) == NULL){printf("erreur ouverture fichier MAPS\n");return 1;}
  if((fmapsGNU = fopen("P_GNU.txt", "w+")) == NULL){printf("erreur ouverture fichier MAPS\n");return 1;}
  if((fV = fopen("Vmaps.txt", "w+")) == NULL){printf("erreur ouverture fichier Vmaps \n");return 1;}
  if((fE = fopen("Emaps.txt", "w+")) == NULL){printf("erreur ouverture fichier Emaps \n");return 1;}


  for(int j=0; j<Ne; j++){
      E=Emin+j*he;
      fprintf(fE, "%.15lf ",E);
  }
  
  double V_unuse_x, E_unuse_x;
  double V_unuse_y, E_unuse_y;

  for(int i=0; i<Nv; i++){
    V=Vmin+i*hv;
    fprintf(fV, "%.15lf ", V);
    for(int j=0; j<Ne; j++){
      //printf("  -j= %d\n",j );
      E=Emin+j*he;
      //printf("V= %.15lf, E=%.15lf\n",V*1e6,E );
      res=intPOLY(V,E,NC,SC);
      if(res){ 
        P=fP(V,E,3);
        T=fT(V,E,3);        
        fprintf(fmaps, "%.15lf %.15lf\n", P, T);
        fprintf(fmapsGNU, "%.15lf %.15lf %.15lf\n", V, E, P);
      }
      else{
        res=intPOLY(V,E,NA,SA);
        if(res){
          P=fP(V,E,1);
          T=fT(V,E,1);
          fprintf(fmaps, "%.15lf %.15lf\n", P, T);
          fprintf(fmapsGNU, "%.15lf %.15lf %.15lf\n", V, E, P);
        }
        else{
          res=intPOLY(V,E,NB,SB);
          if(res){
            P=fP(V,E,2);
            T=fT(V,E,2);
            fprintf(fmaps, "%.15lf %.15lf\n", P, T);
            fprintf(fmapsGNU, "%.15lf %.15lf %.15lf\n", V, E, P);
          }
          else{
            res=ftri(V, E); // triangle triple
            if(res){
              P=PT;
              T=TT;
              fprintf(fmaps, "%.15lf %.15lf\n", P, T);
              fprintf(fmapsGNU, "%.15lf %.15lf %.15lf\n", V, E, P);
            }
            else{
              res=intPOLY(V,E,NAB,SAB);
              if(res){
                /*
                VE_mixte(Nmax, epsilon, 1, 2, V, E, &Vx, &Ex, &Vy, &Ey, &x);
                P=fP(Vx,Ex,1);
                T=fT(Vx,Ex,1);
                */
                err=interp_zone_mixte(V, E, 1, 2, N12, segVa12, segVb12, segEa12, segEb12, &P, &T, &V_unuse_x, &E_unuse_x, &V_unuse_y, &E_unuse_y);
                //printf("P= %lf, T= %lf   mixte 1/2\n",P*1e-9,T );
                fprintf(fmaps, "%.15lf %.15lf\n", P,T);
                fprintf(fmapsGNU, "%.15lf %.15lf %.15lf\n", V, E, P);
              }
              else{
                res=intPOLY(V,E,NAC,SAC);
                if(res){
                  /*
                  VE_mixte(Nmax, epsilon, 1, 3, V, E, &Vx, &Ex, &Vy, &Ey, &x);
                  P=fP(Vx,Ex,1);
                  T=fT(Vx,Ex,1);
                  */
                  err=interp_zone_mixte(V, E, 1, 3, N13, segVa13, segVc13, segEa13, segEc13, &P, &T, &V_unuse_x, &E_unuse_x, &V_unuse_y, &E_unuse_y);
                  //printf("P= %lf, T= %lf   mixte 1/3\n",P*1e-9,T );
                  fprintf(fmaps, "%.15lf %.15lf\n", P,T);
                  fprintf(fmapsGNU, "%.15lf %.15lf %.15lf\n", V, E, P);
                }
                else{
                  res=intPOLY(V,E,NBC,SBC);
                  if(res){
                    /*
                    VE_mixte(Nmax, epsilon, 2, 3, V, E, &Vx, &Ex, &Vy, &Ey, &x);
                    P=fP(Vx,Ex,2);
                    T=fT(Vx,Ex,2);
                    */
                    err=interp_zone_mixte(V, E, 2, 3, N23, segVb23, segVc23, segEb23, segEc23, &P, &T, &V_unuse_x, &E_unuse_x, &V_unuse_y, &E_unuse_y);
                    //printf("P= %lf, T= %lf   mixte 2/3\n",P*1e-9,T );
                    fprintf(fmaps, "%.15lf %.15lf\n", P, T);
                    fprintf(fmapsGNU, "%.15lf %.15lf %.15lf\n", V, E, P);
                  }
                  else{
                    P=-10e9;
                    T=0;
                    fprintf(fmaps, "%.15lf %.15lf\n", P, T);
                    fprintf(fmapsGNU, "%.15lf %.15lf %.15lf\n", V, E, P);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  

  fclose(fmaps);
  fclose(fV);
  fclose(fE);
  fclose(fmapsGNU);

  // Liberation
  free(segP12);   free(segT12);
  free(segVa12);  free(segEa12);
  free(segVb12);  free(segEb12);
  
  free(segP13);   free(segT13);
  free(segVa13);  free(segEa13);
  free(segVc13);  free(segEc13);
  
  free(segP23);   free(segT23);
  free(segVb23);  free(segEb23);
  free(segVc23);  free(segEc23);


  return 0;

}









/* Fonction pour avoir une benchmark de (P,T) sur (V,E)

*/
int ecr_benchmark_PT(int Nv ,int Ne, double Vmin, double Vmax, double Emin, double Emax){

  double he=(Emax-Emin)/(Ne-1);
  double hv=(Vmax-Vmin)/(Nv-1);


  int Nmax=1e4;
  double epsilon=1e-6;
  double Vx, Ex, Vy, Ey, x;
  

  int err,res;
  

  double E,V;

  double Vmx, Emx;
  double Vmy, Emy;


  int nb_points=20;
  int NTG1=nb_points, NPB1=2*nb_points;
  int NPB2=5*nb_points, NTD2=4*nb_points;
  int NTD3=ceil(nb_points/3), NPH3=6*nb_points, NTG3=4*nb_points; 
  int N12=nb_points; int N13=nb_points; int N23=5*nb_points;
  
  int NA=N13+(NTG1-1)+(NPB1-1)+(N12-2);
  double** SA=alloctabd(NA+1,2);
  
  int NB=N12+(NPB2-1)+(NTD2-1)+(N23-2);
  double** SB=alloctabd(NB+1,2);
  
  int NC=N23+(NTD3-1)+(NPH3-1)+(NTG3-1)+(N13-2);
  double** SC=alloctabd(NC+1,2); // on ajoute S[N]=S[0]  --> utile pour l'algo
  
  int NAB=2*N12;
  double** SAB=alloctabd(NAB+1,2);

  int NAC=2*N13;
  double** SAC=alloctabd(NAC+1,2);

  int NBC=2*N23;
  double** SBC=alloctabd(NBC+1,2);


  double Ptriple=4.5104217816091e9;
  double Pdeb=Ptriple;


  int phaseA12=1;
  int phaseB12=2;
  double Pfin12=9.4e9;
  
  int phaseA13=1;
  int phaseB13=3;
  double Pfin13=0.;
  
  int phaseA23=2;
  int phaseB23=3;
  double Pfin23=70e9;
  
  double Va012, Vb012, Va013, Vb013, Va023, Vb023;
  double Ea012, Eb012, Ea013, Eb013, Ea023, Eb023;

  
  // Valeur d'initialisation
  Va012=129.905671875844291e-6;  Ea012=72.168190127265447e3;
  Vb012=127.666112185063994e-6;  Eb012=99.480256119236515e3;
  
  Va013=129.905671875844291e-6;  Ea013=72.168190127265447e3;
  Vb013=132.538478805442338e-6;  Eb013=136.154346228414312e3;
  
  Va023=127.666112185063994e-6;  Ea023=99.480256119236515e3;
  Vb023=132.538478805442338e-6;  Eb023=136.154346228414312e3;
  
   
  // Initialisation des tableaux
  double* segP12=malloc(N12*sizeof(double));
  double* segT12=malloc(N12*sizeof(double));
  double* segVa12=malloc(N12*sizeof(double));
  double* segEa12=malloc(N12*sizeof(double));
  double* segVb12=malloc(N12*sizeof(double));
  double* segEb12=malloc(N12*sizeof(double));
  
  double* segP13=malloc(N13*sizeof(double));
  double* segT13=malloc(N13*sizeof(double));
  double* segVa13=malloc(N13*sizeof(double));
  double* segEa13=malloc(N13*sizeof(double));
  double* segVc13=malloc(N13*sizeof(double));
  double* segEc13=malloc(N13*sizeof(double));
  
  double* segP23=malloc(N23*sizeof(double));
  double* segT23=malloc(N23*sizeof(double));
  double* segVb23=malloc(N23*sizeof(double));
  double* segEb23=malloc(N23*sizeof(double));
  double* segVc23=malloc(N23*sizeof(double));
  double* segEc23=malloc(N23*sizeof(double));
  
  // frontiere 12
  ligne_double(Nmax, epsilon, phaseA12, phaseB12, Pdeb, Pfin12, N12, Va012, Ea012, Vb012, Eb012, segP12, segT12, segVa12, segEa12, segVb12, segEb12);
  
  // frontiere 13
  ligne_double(Nmax, epsilon, phaseA13, phaseB13, Pdeb, Pfin13, N13, Va013, Ea013, Vb013, Eb013, segP13, segT13, segVa13, segEa13, segVc13, segEc13);
  
  // frontiere 23
  ligne_double(Nmax, epsilon, phaseA23, phaseB23, Pdeb, Pfin23, N23, Va023, Ea023, Vb023, Eb023, segP23, segT23, segVb23, segEb23, segVc23, segEc23);
  

  // CALCUL DES SOMMETS  
  err=sommets_polygone(NA, NB, NC, NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3, N12, N13, N23, SA, SB, SC, SAB, SAC, SBC, segP12, segT12, segP13, segT13, segP23, segT23, 
                       segVa12, segEa12, segVb12, segEb12, segVa13, segEa13, segVc13, segEc13, segVb23, segEb23, segVc23, segEc23 );
  if(err){printf("err= %d (sommets_polygone)\n",err ); return err;}

  
  // Pour tester les valeurs
  FILE *fileV, *fileE;
  if((fileV = fopen("Vbench.txt", "w+")) == NULL){printf("erreur ouverture fichier V bench \n"); return 1;}
  if((fileE = fopen("Ebench.txt", "w+")) == NULL){printf("erreur ouverture fichier E bench \n"); return 1;}
  for(int j=0; j<Ne; j++){
      E=Emin+j*he;
      fprintf(fileE, "%.15lf ",E);
  }
  for(int j=0; j<Nv; j++){
      V=Vmin+j*hv;
      fprintf(fileV, "%.15lf ",V);
  }
  fclose(fileV);  fclose(fileE);
  /////////////
  

  double* tabVx=malloc(Nv*Ne*sizeof(double)); double* tabEx=malloc(Nv*Ne*sizeof(double));
  double* tabVy=malloc(Nv*Ne*sizeof(double)); double* tabEy=malloc(Nv*Ne*sizeof(double));
  
  int* tabZONE=malloc(Nv*Ne*sizeof(int));
  

  int k=0;
  for(int i=0; i<Nv; i++){
    V=Vmin+i*hv;
    for(int j=0; j<Ne; j++){
      E=Emin+j*he;
      res=intPOLY(V,E,NC,SC);
      if(res){ 
        tabVx[k]=0; tabEx[k]=-1e-7;
        tabVy[k]=0; tabEy[k]=-1e-7;
        tabZONE[k]=3;
        //printf(" C : P= %lf , T= %lf\n",P*1e-9,T);
        //printf("     tabP= %lf , tabT= %lf\n",tabP[k]*1e-9,tabT[k]);
        /*
        fwrite(&P,sizeof(double),1,fileP);
        fwrite(&T,sizeof(double),1,fileT);
        */
      }
      else{
        res=intPOLY(V,E,NA,SA);
        if(res){
          tabVx[k]=0; tabEx[k]=-1e-7;
          tabVy[k]=0; tabEy[k]=-1e-7;
          tabZONE[k]=1;
          //printf(" A : P= %lf , T= %lf\n",P*1e-9,T);
          //printf("     tabP= %lf , tabT= %lf\n",tabP[k]*1e-9,tabT[k]);
          /*
          fwrite(&P,sizeof(double),1,fileP);
          fwrite(&T,sizeof(double),1,fileT);
          */
        }
        else{
          res=intPOLY(V,E,NB,SB);
          if(res){
            tabVx[k]=0; tabEx[k]=-1e-7;
            tabVy[k]=0; tabEy[k]=-1e-7;
            tabZONE[k]=2;
            //printf(" B : P= %lf , T= %lf\n",P*1e-9,T);
            //printf("     tabP= %lf , tabT= %lf\n",tabP[k]*1e-9,tabT[k]);            
            /*
            fwrite(&P,sizeof(double),1,fileP);
            fwrite(&T,sizeof(double),1,fileT);
            */
          }
          else{
            res=ftri(V, E); // triangle triple
            if(res){
              tabVx[k]=0; tabEx[k]=-1e-7;
              tabVy[k]=0; tabEy[k]=-1e-7;
              tabZONE[k]=0;
              //printf(" Triangle : P= %lf , T= %lf\n",P*1e-9,T);
              //printf("            tabP= %lf , tabT= %lf\n",tabP[k]*1e-9,tabT[k]);
              /*
              fwrite(&P,sizeof(double),1,fileP);
              fwrite(&T,sizeof(double),1,fileT);
              */
            }
            else{
              res=intPOLY(V,E,NAB,SAB);
              if(res){
                err=graine_zone_mixte(V, E, N12, segVa12, segVb12, segEa12, segEb12, &Vmx, &Emx, &Vmy, &Emy);
                if(err){printf("err= %d (benchmark_PT)\n",err ); return err;}
                tabVx[k]=Vmx; tabEx[k]=Emx;
                tabVy[k]=Vmy; tabEy[k]=Emy;
                tabZONE[k]=12;
              }
              else{
                res=intPOLY(V,E,NAC,SAC);
                if(res){
                  err=graine_zone_mixte(V, E, N13, segVa13, segVc13, segEa13, segEc13, &Vmx, &Emx, &Vmy, &Emy);
                  if(err){printf("err= %d (benchmark_PT)\n",err ); return err;}
                  tabVx[k]=Vmx; tabEx[k]=Emx;
                  tabVy[k]=Vmy; tabEy[k]=Emy;
                  tabZONE[k]=13;
                }
                else{
                  res=intPOLY(V,E,NBC,SBC);
                  if(res){
                    err=graine_zone_mixte(V, E, N23, segVb23, segVc23, segEb23, segEc23, &Vmx, &Emx, &Vmy, &Emy);
                    if(err){printf("err= %d (benchmark_PT)\n",err ); return err;}
                    tabVx[k]=Vmx; tabEx[k]=Emx;
                    tabVy[k]=Vmy; tabEy[k]=Emy;
                    tabZONE[k]=23;
                  }
                  else{
                    tabVx[k]=0; tabEx[k]=-1e-7;
                    tabVy[k]=0; tabEy[k]=-1e-7;
                    tabZONE[k]=-1;
                  }
                }
              }
            }
          }
        }
      }
      k++;
    }
  }
  
  /*
  for(int i=0; i<Ne*Nv; i++){
    printf(" i= %d | P= %lf , T= %f\n",i,tabP[i]*1e-9,tabT[i] );
  }
  */

  // fichier binaire    couple (V,E) pour les zones mixte
  FILE *fileVmx, *fileEmx;
  if((fileVmx = fopen("Vbench_mixte_x.bin", "wb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmx = fopen("Ebench_mixte_x.bin", "wb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  FILE *fileVmy, *fileEmy;
  if((fileVmy = fopen("Vbench_mixte_y.bin", "wb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmy = fopen("Ebench_mixte_y.bin", "wb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  
  // Ecriture des tableaux
  fwrite(tabVx,sizeof(double),Ne*Nv,fileVmx);
  fwrite(tabEx,sizeof(double),Ne*Nv,fileEmx);
  fwrite(tabVy,sizeof(double),Ne*Nv,fileVmy);
  fwrite(tabEy,sizeof(double),Ne*Nv,fileEmy);

  fclose(fileVmx);  fclose(fileEmx);
  fclose(fileVmy);  fclose(fileEmy);
  
  // fichier binaire ZONE
  FILE *fileZONE;
  if((fileZONE = fopen("ZONEbench.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE bench\n"); return 1;}
  
  // Ecriture des tableaux
  fwrite(tabZONE,sizeof(int),Ne*Nv,fileZONE);

  fclose(fileZONE);


  // Liberation
  free(tabVx);   free(tabEx);
  free(tabVy);   free(tabEy);
  free(tabZONE);


  free(segP12);   free(segT12);
  free(segVa12);  free(segEa12);
  free(segVb12);  free(segEb12);
  
  free(segP13);   free(segT13);
  free(segVa13);  free(segEa13);
  free(segVc13);  free(segEc13);
  
  free(segP23);   free(segT23);
  free(segVb23);  free(segEb23);
  free(segVc23);  free(segEc23);


  return 0;

}


/* Lire le benchmark

   le fichier est en binaire
*/
int lec_benchmark_PT(double V, double E, double* tabVx, double* tabEx, double* tabVy, double* tabEy, int* tabZONE, int Nv, int Ne, double Vmin, double Vmax, double Emin, double Emax, double* pP, double* pT, double* px){
  int Nmax=1e4;
  double epsilon=1e-5;

  int err; 
  double P,T;
  double Vmx, Emx, Vmy, Emy;
  double Vx, Ex, Vy, Ey;

  double x;
  
  // Point triple
  double PT= 4510421781.609182357788085;
  double TT= 624.63982538492518;

  if(V<Vmin || Vmax<V){printf("V hors domaine (lec_benchmark_PT)\n");   return 1;}  
  if(E<Emin || Emax<E){printf("E hors domaine (lec_benchmark_PT)\n");   return 1;}
 
  double hv=(Vmax-Vmin)/(Nv-1);
  double he=(Emax-Emin)/(Ne-1);
  int i=(V-Vmin)/hv;
  int j=(E-Emin)/he;
  //printf("  ** i= %d , j= %d  (lec_benchmark_PT)\n",i,j );
  
  int pos_fic=i*Ne + j;
  /* INUTILE
  P=tabP[pos_fic]; 
  T=tabT[pos_fic]; 
  */
  Vmx=tabVx[pos_fic];  Emx=tabEx[pos_fic];
  Vmy=tabVy[pos_fic];  Emy=tabEy[pos_fic];
  
  int ind_ZONE=tabZONE[pos_fic];


  if(ind_ZONE==-1){
    P=-10e9; T=0.;
    x=-1;
  }
  else if(ind_ZONE==3){
    P=fP(V,E,3);  T=fT(V,E,3);
    x=-1;
  }
  else  if(ind_ZONE==2){
    P=fP(V,E,2);  T=fT(V,E,2);
    x=-1;
  }
  else if(ind_ZONE==1){
    P=fP(V,E,1);  T=fT(V,E,1);
    x=-1;
  }
  else if(ind_ZONE==23){
    VE_mixte(Nmax, epsilon, 2, 3, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, &x);
    P=fP(Vx,Ex,2);  T=fT(Vx,Ex,2);
  }
  else if(ind_ZONE==13){
    VE_mixte(Nmax, epsilon, 1, 3, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, &x);
    P=fP(Vx,Ex,1);  T=fT(Vx,Ex,1);
  }
  else if(ind_ZONE==12){
    VE_mixte(Nmax, epsilon, 1, 2, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, &x);
    P=fP(Vx,Ex,1);  T=fT(Vx,Ex,1);
  }
  else if(ftri(V, E)){
    P=PT;
    T=TT;
  }
  else{printf("Erreur dans le choix de ind_err (\n"); }
  
  *pP=P;
  *pT=T;
  *px=x;

  return 0;
}




// Script pour tester le Benchmark
int use_benchmark_PT(){
  int err;
  double E,V;
  double P,T;
  double x;
  
  // Discrétisation du benchmark/tabulation
  int Nv_dis=4e2;
  int Ne_dis=4e2;

  // Discrétisation du test
  int Nv_test=4e2;
  int Ne_test=4e2;

  double Vmin_dis, Vmax_dis;
  double Emin_dis, Emax_dis;

  double Vmin_test, Vmax_test;
  double Emin_test, Emax_test;
  
  Vmin_dis=80e-6, Vmax_dis=210e-6;
  Emin_dis=-2e3,  Emax_dis=16e5;

  Vmin_test=110e-6, Vmax_test=150e-6;
  Emin_test=-2e3,   Emax_test=4e5;


  double hv_test=(Vmax_test-Vmin_test)/(Nv_test-1);
  double he_test=(Emax_test-Emin_test)/(Ne_test-1);
  
  printf(" deb ecr_benchmark_PT\n");
  err=ecr_benchmark_PT(Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis);
  if(err){return err;}
  printf(" fin ecr_benchmark_PT\n");


  //  ****************************************
  // Pour tester les valeurs
  FILE *fPtest, *fTtest, *fXtest;
  if((fPtest = fopen("Ptest.txt", "w+")) == NULL){printf("erreur ouverture fichier P test\n");return 1;}
  if((fTtest = fopen("Ttest.txt", "w+")) == NULL){printf("erreur ouverture fichier T test\n");return 1;}
  if((fXtest = fopen("Xtest.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  FILE *fVtest, *fEtest;
  if((fVtest = fopen("Vtest.txt", "w+")) == NULL){printf("erreur ouverture fichier V test \n");return 1;}
  if((fEtest = fopen("Etest.txt", "w+")) == NULL){printf("erreur ouverture fichier E test \n");return 1;}
  for(int j=0; j<Ne_test; j++){
      E=Emin+j*he_test;
      fprintf(fEtest, "%.15lf ",E);
  }
  for(int j=0; j<Nv_test; j++){
      V=Vmin+j*hv_test;
      fprintf(fVtest, "%.15lf ",V);
  }
  fclose(fVtest);  fclose(fEtest);
  //  ****************************************


  printf("ouverture V,E\n");
  // fichier binaire    couple (V,E) pour les zones mixte
  FILE *fileVmx, *fileEmx, *fileVmy, *fileEmy;
  if((fileVmx = fopen("Vbench_mixte_x.bin", "rb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmx = fopen("Ebench_mixte_x.bin", "rb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  if((fileVmy = fopen("Vbench_mixte_y.bin", "rb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmy = fopen("Ebench_mixte_y.bin", "rb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  
  double* tabVx=malloc(Nv_dis*Ne_dis*sizeof(double)); double* tabEx=malloc(Nv_dis*Ne_dis*sizeof(double));
  double* tabVy=malloc(Nv_dis*Ne_dis*sizeof(double)); double* tabEy=malloc(Nv_dis*Ne_dis*sizeof(double));

  // lecture des tableaux
  fread(tabVx,sizeof(double),Ne_dis*Nv_dis,fileVmx);
  fread(tabEx,sizeof(double),Ne_dis*Nv_dis,fileEmx);
  fread(tabVy,sizeof(double),Ne_dis*Nv_dis,fileVmy);
  fread(tabEy,sizeof(double),Ne_dis*Nv_dis,fileEmy);

  fclose(fileVmx);  fclose(fileEmx);
  fclose(fileVmy);  fclose(fileEmy);
  
  printf("ouverture ZONE\n");
  // Ouverture du fichier ZONE
  FILE *fileZONE;
  if((fileZONE = fopen("ZONEbench.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench\n");return 2;}
  rewind(fileZONE); 
  int* tabZONE=malloc(Nv_dis*Ne_dis*sizeof(int));
  fread(tabZONE,sizeof(int),Ne_dis*Nv_dis,fileZONE);
  fclose(fileZONE);




  // ***** Boucle *****
  for(int i=0; i<Nv_test; i++){
    printf("  - i= %d\n",i );
    V=Vmin+i*hv_test;
    for(int j=0; j<Ne_test; j++){
      //printf(" - i= %d, j= %d \n",i,j );
      E=Emin+j*he_test;
      //err=lec_benchmark_PT(V, E, fileP, fileT, Nv, Ne, Vmin, Vmax, Emin, Emax, &P, &T);
      err=lec_benchmark_PT(V, E, tabVx, tabEx, tabVy, tabEy, tabZONE, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, &P, &T, &x);

      if(err){return err;}
      //printf("     P= %lf , T= %f\n",P*1e-9,T );
      fprintf(fPtest,"%.10lf\n",P);
      fprintf(fTtest,"%.10lf\n",T);
      fprintf(fXtest,"%.10lf\n",x);
    }
  }
  

  fclose(fPtest);  fclose(fTtest);  fclose(fXtest);

  return 0;
}



/*   CHerche la zone grace aux enthalpies
  Processus  :
  Pour un couple (V,E) donné, 
  De la on peut calcuer les 3 enthalpie corespondandes :
    - Si les 3 sont égales on est dans le triangle
    - Si 2 seulements sont égales, on est dans une zones de mélanges coresp. au phase qui ont le meme G
    - Si les 3 enthalpies sont diffénrentes aloron est dans une zone de phase pur coresp. ala phase au G minimum é

   Risque de ne pas fonctionner : pourquoi les enthalipes seraient égales dans le triangle ???
         A TESTER
*/
/*
int minG(double V, double E, double* pP, double* pT, double* pX){
  double epsilon=1e-7;

  double VA, EA, VB, EB, VC, EC;

  double Pa, Ta, Sa, Ga;
  double Pb, Tb, Sb, Gb;
  double Pc, Tc, Sc, Gc;

  // phase 1/A
  Pa=fP(V,E, 1); Ta=fT(V,E, 1);
  Sa=fS(V,E, 1); Ga=E+Pa*V-Ta*Sa;
  //phase 2/B
  Pb=fP(V,E, 2); Tb=fT(V,E, 2);
  Sb=fS(V,E, 2); Gb=E+Pb*V-Tb*Sb;
  // phase 3/C
  Pc=fP(V,E, 3); Tc=fT(V,E, 3);
  Sc=fS(V,E, 3); Gc=E+Pc*V-Tc*Sc;

  if(isnan(Ga) && isnan(Gb) && isnan(Gc)){printf("Hors domaine\n"); return 1;}


  if(isnan(Ga)){Ga=1e70;}
  if(isnan(Gb)){Gb=1e70;}
  if(isnan(Gc)){Gc=1e70;}

  // Point triple
  if (fabs(Ga-Gb)<epsilon && fabs(Ga-Gc)<epsilon){
    P=fP(V,E,1)
  }




  return 0;
}
*/


/*  

*/
/*
void Newton_all(double* A, double* InvA, double Vmel, double Emel, double Va, double Ea, double Vb, double Eb, double Vc, double Ec, double xa, double xb, double* pVA, double* pEA, double* pVB, double* pEB, double* pVC, double* pEC, double* pXa, double* pXb, double* dX){
  int dim=8;

  // Grandeurs thermodynamiques
  double Pa, Ta, Sa, Ga, dPva, dPea, dTva, dTea;
  double Pb, Tb, Sb, Gb, dPvb, dPeb, dTvb, dTeb;
  double Pc, Tc, Sc, Gc, dPvc, dPec, dTvc, dTec;
  // phase 1/A
  Pa=fP(Va,Ea, 1); Ta=fT(Va,Ea, 1);
  Sa=fS(Va,Ea, 1); Ga=Ea+Pa*Va-Ta*Sa;
  dPva=dPdV(Va, 1); dPea=dPdE(1);
  dTva=dTdV(Va, 1); dTea=dTdE(1);
  // phase 2/B
  Pb=fP(Vb,Eb, 2); Tb=fT(Vb,Eb, 2);
  Sb=fS(Vb,Eb, 2); Gb=Eb+Pb*Vb-Tb*Sb;
  dPvb=dPdV(Vb,2); dPeb=dPdE(2);
  dTvb=dTdV(Vb,2); dTeb=dTdE(2);
  // phase 3/C
  Pc=fP(Vc,Ec, 3); Tc=fT(Vc,Ec, 3);
  Sc=fS(Vc,Ec, 3); Gc=Ec+Pc*Vc-Tc*Sc;
  dPvc=dPdV(Vc,3); dPec=dPdE(3);
  dTvc=dTdV(Vc,3); dTec=dTdE(3);
  
  
  double dGvx, dGex, dGvb, dGeb, dGvc, dGec;
  dGvx=Vx*dPvx - Sx*dTvx; 
  dGex=Vx*dPex - Sx*dTex;
  
  dGvb=Vb*dPvb - Sb*dTvb; 
  dGeb=Vb*dPeb - Sb*dTeb;

  dGvc=Vc*dPvc - Sc*dTvc; 
  dGec=Vc*dPec - Sc*dTec;
  
  
  double* Delta=malloc(dim*sizeof(double));
  Delta[0]=Pb-Pa;   Delta[1]=Pc-Pa;   
  Delta[2]=Tb-Ta;   Delta[3]=Tc-Ta;   
  Delta[4]=Gb-Ga;   Delta[5]=Gc-Ga;  
  Delta[6]=Vmel - xa*Va - xb*Vb - (1.-xa-xb)*Vc;
  Delta[7]=Emel - xa*Ea - xb*Eb - (1.-xa-xb)*Ec;
   
  /*
  for(int j=0; j<dim; j++){
    printf("Delta[%d]= %.15lf  ",j, Delta[j] );
  }
  printf("\n");
  */
  
  /*
  // on met le second menbre dans la dernière colonne
  //matrice A
  A[dim*0+0]=dPvx; A[dim*0+1]=dPex; A[dim*0+2]=-dPvy; A[dim*0+3]=-dPey; A[dim*0+4]=0.; A[dim*0+5]=dPvx; A[dim*0+6]=0; A[dim*0+7]=0; 
  
  A[dim*1+0]=dTvx; A[dim*1+1]=dTex; A[dim*1+2]=-dTvy; A[dim*1+3]=-dTey; A[dim*1+4]=0.; A[dim*1+5]=dPvx; A[dim*1+6]=0; A[dim*1+7]=0; 
  
  A[dim*2+0]=dGvx; A[dim*2+1]=dGex; A[dim*2+2]=-dGvy; A[dim*2+3]=-dGey; A[dim*2+4]=0.; A[dim*2+5]=dPvx; A[dim*2+6]=dPex; A[dim*2+7]=0; 
  
  A[dim*3+0]=1.-x;  A[dim*3+1]=0.;   A[dim*3+2]=x;     A[dim*3+3]=0;     A[dim*3+4]=Vy-Vx; A[dim*3+5]=dPvx; A[dim*3+6]=dPex; A[dim*3+7]=-dPvy; 
  
  A[dim*4+0]=0.;   A[dim*4+1]=1.-x;  A[dim*4+2]=0;     A[dim*4+3]=x;     A[dim*4+4]=Ey-Ex; A[dim*4+5]=dPvx; A[dim*4+6]=dPex; A[dim*4+7]=-dPvy;

  A[dim*5+0]=0.;   A[dim*5+1]=1.-x;  A[dim*5+2]=0;     A[dim*5+3]=x;     A[dim*5+4]=Ey-Ex; A[dim*5+5]=dPvx; A[dim*5+6]=dPex; A[dim*5+7]=-dPvy; 

  A[dim*6+0]=0.;   A[dim*6+1]=1.-x;  A[dim*6+2]=0;     A[dim*6+3]=x;     A[dim*6+4]=Ey-Ex; A[dim*6+5]=dPvx; A[dim*6+6]=dPex; A[dim*6+7]=-dPvy; 
  
  A[dim*7+0]=0.;   A[dim*7+1]=1.-x;  A[dim*7+2]=0;     A[dim*7+3]=x;     A[dim*7+4]=Ey-Ex; A[dim*7+5]=dPvx; A[dim*7+6]=dPex; A[dim*7+7]=-dPvy; 
  
  /*
  printf("  * A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      printf("A[%d][%d]= %.15lf  ",i,j, A[dim*i+j] );
    }
    printf("\n");
  }
  printf("\n");  
  */

/*
  double det_a=inverse_matrice_pivot(A, dim, InvA);
  //printf("detA= %.15lf\n",det_a );
  if ( fabs(det_a)<1e-8){printf("A non inversible : detA= %.15lf\n", det_a);}
  
  /*
  printf("  * InvA = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      printf("InvA[%d][%d]= %.15lf  ",i,j, InvA[dim*i+j] );
    }
    printf("\n");
  }
  printf("\n");
  */
  
  
  // Mutiplcation matricielle de A et InvA
  /*
  double res;
  printf("  * InvA*A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      res=0;
      for(int k=0; k<dim; k++){
        res+=InvA[dim*i+k]*A[dim*k+j];
      }
      printf("A*InvA[%d][%d]= %.15lf ",i,j,res );
    }
    printf("\n");
  }
  printf("\n");
  */

/*
  // Mutiplcation matricielle
  for(int i=0; i<dim; i++){
    dX[i]=0;
    for(int j=0; j<dim; j++){
      dX[i]+=InvA[dim*i+j]*Delta[j];
    }
    //printf("  dX[%d]= %.15lf\n",i,dX[i] );
  }
  //printf("\n");
   
  double correc=1.;

  *pVX=Vx+correc*dX[0];
  *pEX=Ex+correc*dX[1];
  *pVY=Vy+correc*dX[2];
  *pEY=Ey+correc*dX[3];
  *pX=x+correc*dX[4];

  /*
  printf("  Vx= %.10lf, Ex= %.10lf   Vy= %.10lf, Ey= %.10lf    x= %.10lf\n",Vx,Ex,Vy,Ey,x );
  printf("  (Vmel-Vx)/(Vy-Vx)= %.10lf  ,    (Emel-Ex)/(Ey-Ex)= %.10lf\n", (Vmel-Vx)/(Vy-Vx), (Emel-Ex)/(Ey-Ex));
  printf("  (1-x)*Vx + x*Vy - Vmel= %.10lf  ,    (1-x)*Ex + x*Ey - Emel= %.10lf\n", (1-x)*Vx + x*Vy - Vmel, (1-x)*Ex + x*Ey - Emel);
  */
/*
  free(Delta);
}
*/

int main(){
  int err;
  
  //     Calcul des fonction thermo pour chaque phase
  //void courb_thermo(void);
  
  //     Calcul du Point TRIPLE
  //void point_triple(void);
  
  //     Calcul du Point courbe frontières
  //err= diag_phase();
  
  // ZONE 
  // repere les zones (triangle pour le moment)
  //err=diag_zone();
  
  
  // Frontière VE
  //err=routine_VE();
  //if(err){return err;}
  
  
  //**  MAILLAGE  **
  //err=maillage();
  
  
  // Zone mixte
  /*
  mapsPT();
  
  int Nmax=1e4;
  double epsilon=1e-7;
  int phaseX=2;
  int phaseY=3;
  double Vmel=110e-6;
  double Emel=4.768e5;
  double Vx, Ex, Vy, Ey, x;
  VE_mixte(Nmax, epsilon, phaseX, phaseY, Vmel, Emel, &Vx, &Ex, &Vy, &Ey, &x);
  
  printf("phaseX= %d\n",phaseX );
  printf("Vx= %.10lf, Ex= %.10lf\n",Vx,Ex);
  printf("phaseY= %d\n",phaseY );
  printf("Vy= %.10lf, Ey= %.10lf\n",Vy,Ey);
  printf(" x= %.10lf\n",x );
  */

  // Map de la Pression et de la temperature e, fct de (V,E)
  //err=mapsPT();
  

  // Pression et température du point triple
  /*
  double VAt= 129.905671875844291e-6;
  double EAt= 72.168190127265447e3;
  double PT=fP(VAt, EAt,1);
  double TT=fT(VAt, EAt,1);
  printf("P= %.20lf  ,  T= %.20lf\n",PT,TT);
  */

  err=use_benchmark_PT();
  if(err){return err;}

}
