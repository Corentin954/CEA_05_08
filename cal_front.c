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
void Newton_triple(double* A, double* InvA, double* Delta, double* dX, double Va, double Ea, double Vb, double Eb, double Vc, double Ec, double* VA, double* EA, double* VB, double* EB, double* VC, double* EC){
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
  double* Delta = (double *) malloc(sizeof(double) * (dim));
  double* dX = (double *) malloc(sizeof(double) * (dim));

  double VA, EA, VB, EB, VC, EC;

  double Pa, Ta, Sa, Ga, Ea, Va;
  double Pb, Tb, Sb, Gb, Eb, Vb;
  double Pc, Tc, Sc, Gc, Ec, Vc;


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


  double critere=1.;

  while(critere>=epsilon && n<Nmax){
    //printf("n= %d\n",n );
    Newton_triple(A, InvA, Delta, dX, Va, Ea, Vb, Eb, Vc, Ec, &VA, &EA, &VB, &EB, &VC, &EC);

    Va=VA; Ea=EA;
    Vb=VB; Eb=EB;
    Vc=VC; Ec=EC;
     
    /*
    // phase 1
    Pa=fP(Va,Ea, 1); Ta=fT(Va,Ea, 1);
    Sa=fS(Va,Ea, 1); Ga=Ea+Pa*Va-Ta*Sa;
    //phase 2
    Pb=fP(Vb,Eb, 2); Tb=fT(Vb,Eb, 2);
    Sb=fS(Vb,Eb, 2); Gb=Eb+Pb*Vb-Tb*Sb;
    // phase 3
    Pc=fP(Vc,Ec, 3); Tc=fT(Vc,Ec, 3);
    Sc=fS(Vc,Ec, 3); Gc=Ec+Pc*Vc-Tc*Sc;
    */

    //critere=fabs(Ga-Gb)+fabs(Ga-Gc);
    critere=fabs(dX[0]/Va) + fabs(dX[1]/Ea) + fabs(dX[2]/Vb) + fabs(dX[3]/Eb) + fabs(dX[4]/Vc) + fabs(dX[5]/Ec);
    printf(" critere= %.14lf\n",critere );

    n++;
  }


  *Va0=Va; *Ea0=Ea; 
  *Vb0=Vb; *Eb0=Eb; 
  *Vc0=Vc; *Ec0=Ec;

  free(A); free(InvA);
  free(Delta);
  free(dX);
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
void Newton_double(double* A, double* InvA, double* Delta, double* dX, double Pcst, int phaseA, int phaseB, double Va, double Ea, double Vb, double Eb, double* VA, double* EA, double* VB, double* EB){
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
  for(int i=0; i<dim; i++){
    dX[i]=0;
    for(int j=0; j<dim; j++){
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


}



// Fonction qui renvoie le nombre de points sur les frontières en focntions d'un seul entier
void fnb_pts(int nb_points, int* pNTG1, int* pNPB1, int* pNPB2, int* pNTD2 , int* pNTD3, int* pNPH3, int* pNTG3, 
            int* pN12 ,int* pN13, int* pN23, int* pNA, int* pNB, int* pNC, int* pNAB, int* pNAC, int* pNBC ){
  

  int nb_points_front=ceil(nb_points); // /3    La dicretisation des frontières a moins d'importance NOP
  // A
  *pNTG1=nb_points_front;
  *pNPB1=2*nb_points;
  // B
  *pNPB2=12*nb_points_front;
  *pNTD2=3*nb_points_front;
  // C
  *pNTD3=2*nb_points_front;
  *pNPH3=15*nb_points_front;
  *pNTG3=8*nb_points_front;
  
  *pN12=5*nb_points; *pN13=3*nb_points; *pN23=40*nb_points;

  *pNA=*pN13+(*pNTG1-1)+(*pNPB1-1)+(*pN12-2);
  *pNB=*pN12+(*pNPB2-1)+(*pNTD2-1)+(*pN23-2);
  *pNC=*pN23+(*pNTD3-1)+(*pNPH3-1)+(*pNTG3-1)+(*pN13-2);
  *pNAB=2*(*pN12);
  *pNAC=2*(*pN13);
  *pNBC=2*(*pN23);
}



//  *******  Fonction qui calcule les points (V,E) des lignes de chgmt de phases
/*
Arguments d'entree :
  - Nmax : nombre max d'iteration
  - epsilon : tolérance sur le critère  (ici fabs(G1-G2))
  - [Pdeb, Pfin] : segment de presion sur lequel calcule la zone mixte
  - Np : nomnbre de points de discretisatoin du segment
çgoukjb// *********    fonction d'inversion    ************
/*
    Obtention de (V,E) à partir de (P,T)
    NEwton classique 
*/
int fVE(int phase, double P, double T, double* pV, double* pE){
  int Nmax=1e2;
  double epsilon=1e-10;
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


  int nb_points=20;
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
  if((fres1 = fopen("fichiers/VEfront1.txt", "w+")) == NULL){printf("erreur ouverture fichier (routine_VE)\n");return 1;}
  if((fres2 = fopen("fichiers/VEfront2.txt", "w+")) == NULL){printf("erreur ouverture fichier (routine_VE)\n");return 1;}
  if((fres3 = fopen("fichiers/VEfront3.txt", "w+")) == NULL){printf("erreur ouverture fichier (routine_VE)\n");return 1;}

  

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
    //printf("%.15lf %.15lf\n",V*1e6,E );
    fprintf(fres1, "%.15lf %.15lf\n",V,E );
  }
  T=300;
  Pdeb=0;
  Pfin=9.4e9;
  for(int i=0; i<N; i++){
    P=Pdeb + i*(Pfin-Pdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    //printf("%.15lf %.15lf\n",V*1e6,E );
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
    //printf("%.5lf %.5lf\n",V*1e6,E );
    fprintf(fres2, "%.15lf %.15lf\n",V,E );
  }
  P=70e9;
  Tdeb=300.;
  Tfin=2486.138962087622531;
  for(int i=0; i<N; i++){
    T=Tdeb + i*(Tfin-Tdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    //printf("%.5lf %.5lf\n",V*1e6,E );
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
    //printf("%.15lf %.15lf\n",V*1e6,E );
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
  if((fres1 = fopen("fichiers/VEfront1.txt", "w+")) == NULL){printf("erreur ouverture fichier (routine_VE)\n");return 1;}
  if((fres2 = fopen("fichiers/VEfront2.txt", "w+")) == NULL){printf("erreur ouverture fichier (routine_VE)\n");return 1;}
  if((fres3 = fopen("fichiers/VEfront3.txt", "w+")) == NULL){printf("erreur ouverture fichier (routine_VE)\n");return 1;}


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
int sommetsC(int N13, int N23, int NTD, int NPH, int NTG, double PrefGB, double PrefGH, double PrefDH, double PrefD, double TrefG, double TrefD, double Trefhaut, double* segVc13, double* segEc13, double* segVc23, double* segEc23 ,double** S){
  

  int N=N23+(NTD-1)+(NPH-1)+(NTG-2)+(N13-1) ;

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
  double hPD=(PrefDH-PrefD)/(NTD-1);
  
  double hPH=(PrefGH-PrefDH)/(NPH-1);
  // la frontière gauche n'est pas horizontale
  double hTG=(TrefG-Trefhaut)/(NTG-1);
  double hPG=(PrefGB-PrefGH)/(NTG-1);


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
  
  for(int j=1; j<=NTD-1; j++){ 
    T=TrefD+j*hTD;
    P=PrefD+j*hPD;
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsC)\n");return err;}
    S[nref+j-1][0]=V;
    S[nref+j-1][1]=E;
  }

  nref+=NTD-1;
  
  // Parcours de la frontière haute
  T=Trefhaut;
  for(int j=1; j<=NPH-1; j++){ 
    P=PrefDH+j*hPH;
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsC)\n");return err;}
    S[nref+j-1][0]=V;
    S[nref+j-1][1]=E;
  }
  nref+=NPH-1;

  /*
      Calculer les points et les mettre sur la frontières
  */



  // Parcours de la frontière gauche
  for(int j=1; j<=NTG-2; j++){ // N13 points
    T=Trefhaut+j*hTG;
    P=PrefGH+j*hPG;
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



// Test si le segment [(xa,ya),(xb,yb)] s'intersecte avec un des segment du polygone S
/*
   Renvoie  1 si intersection
   et       0 si pas d'intersection 
           -1 en cas d'erreur 

*/


int inter_seg(double** REC, int N, double** S, int* pInd){
  
  double penteP, pente;
  double xi,yi,xi1,yi1;
  double xI,yI; // coordonées de l'intersection
  
  double xa, ya, xb, yb;

  for(int ind_rec=0; ind_rec<4; ind_rec++){
    //printf(" ind_rec= %d\n",ind_rec);

    xa=REC[ind_rec][0];    ya=REC[ind_rec][1];
    xb=REC[ind_rec+1][0];  yb=REC[ind_rec+1][1];

    if(xa==xb && ya==yb){printf("deux points du rectangle sont egaux (inter_seg)\n"); return -1;}

    //printf("    xa= %lf  ,   ya= %lf\n",xa,ya );
    //printf("    xb= %lf  ,   yb= %lf\n",xb,yb );

    if(xa==xb){ // pente du segment P non-défini
      xI=xa; // segement vertical
      for(int i=0; i<N; i++){
        xi=S[i][0]; xi1=S[i+1][0];
        yi=S[i][1]; yi1=S[i+1][1];
        
        //printf("   --i= %d\n",i );
        //printf("       xi= %lf  ,   yi= %lf\n",xi,yi );
        //printf("       xi1= %lf  ,   yi1= %lf\n",xi1,yi1 );
        
        if(xi1==xi){printf(" segment tout deux verticaux : xi1= %g, xi= %g, yi1= %g, yi= %g\n",xi1,xi,yi1,yi); } // mettre une commande pour passer au i suivant 
        else{
          pente=(yi1-yi)/(xi1-xi);
          yI=pente*(xI-xi)+yi;  
          //printf("       xI= %lf  ,   yI= %lf\n",xI,yI );
          //printf("       (xI-xi)= %lf , (xI-xi1)= %lf\n",(xI-xi), (xI-xi1) );
          //printf("       (yI-yi)= %lf , (yI-yi1)= %lf\n",(yI-yi), (yI-yi1) );
          //printf("       (xI-xa)= %lf , (xI-xb)= %lf\n",(xI-xa), (xI-xb) );
          //printf("       (yI-ya)= %lf , (yI-yb)= %lf\n",(yI-ya), (yI-yb) );
          if( (yI-ya)*(yI-yb)<0  &&  (xI-xi)*(xI-xi1)<0 ){ // on sait que ya!=yb
            if(yi==yi1){*pInd=i; return 1;}
            else{
              if((yI-yi)*(yI-yi1)<0){*pInd=i; return 1;} 
            }
          }
        }  
      }
    }
    else{ // pente défini car xa=!xb
      penteP=(yb-ya)/(xb-xa);
      for(int i=0; i<N; i++){
        xi=S[i][0]; xi1=S[i+1][0];
        yi=S[i][1]; yi1=S[i+1][1];
        //printf("   --i= %d\n",i );
        //printf("        xi= %lf  ,  yi= %lf\n",xi,yi );
        //printf("        xi1= %lf ,  yi1= %lf\n",xi1,yi1 );
        if(xi1==xi){ // pas de pente pour le segment S
          xI=xi;
          yI=penteP*(xI-xa)+ya;
          //printf("       xI= %lf  ,   yI= %lf\n",xI,yI );
          //printf("       (xI-xi)= %lf , (xI-xi1)= %lf\n",(xI-xi), (xI-xi1) );
          //printf("       (yI-yi)= %lf , (yI-yi1)= %lf\n",(yI-yi), (yI-yi1) );
          //printf("       (xI-xa)= %lf , (xI-xb)= %lf\n",(xI-xa), (xI-xb) );
          //printf("       (yI-ya)= %lf , (yI-yb)= %lf\n",(yI-ya), (yI-yb) );
          if( (xI-xa)*(xI-xb)<0  &&  (yI-yi)*(yI-yi1)<0 ){  // on sait que yi!=yi1
            if(ya==yb){*pInd=i; return 1;}
            else{
              if((yI-ya)*(yI-yb)<0){*pInd=i; return 1;}
            }
          }
        }
        else{ // les deux segment ont une pente  ouf !!
          pente=(yi1-yi)/(xi1-xi);
          xI=(yi-ya+penteP*xa-pente*xi)/(penteP-pente);
          yI=pente*(xI-xi)+yi; 
          //printf("       xI= %lf  ,   yI= %lf\n",xI,yI );
          //printf("       (xI-xi)= %lf , (xI-xi1)= %lf\n",(xI-xi), (xI-xi1) );
          //printf("       (yI-yi)= %lf , (yI-yi1)= %lf\n",(yI-yi), (yI-yi1) );
          //printf("       (xI-xa)= %lf , (xI-xb)= %lf\n",(xI-xa), (xI-xb) );
          //printf("       (yI-ya)= %lf , (yI-yb)= %lf\n",(yI-ya), (yI-yb) );
          //if( (xI-xa)*(xI-xb)<0 &&  (yI-ya)*(yI-yb)<0  &&  (xI-xi)*(xI-xi1)<0 &&  (yI-yi)*(yI-yi1)<0 ){return 1;}
          
          if(penteP==pente){} // segment parralelle //
          else if(yi1==yi){
            if( (xI-xa)*(xI-xb)<0 &&  (yI-ya)*(yI-yb)<0  &&  (xI-xi)*(xI-xi1)<0 ){*pInd=i; return 1;}
          }
          else if(ya==yb){
            if( (xI-xa)*(xI-xb)<0 &&  (xI-xi)*(xI-xi1)<0 &&  (yI-yi)*(yI-yi1)<0 ){*pInd=i; return 1;}
          }
          else{
            if( (xI-xa)*(xI-xb)<0 &&  (yI-ya)*(yI-yb)<0  &&  (xI-xi)*(xI-xi1)<0 &&  (yI-yi)*(yI-yi1)<0 ){*pInd=i; return 1;}
          }
        }
      }
    }
  }


  return 0;  // pas d'intsection

}




/*  
    Creation des fichiers de maillage
*/
int sommets_polygone(int NA, int NB, int NC, int NTG1, int NPB1, int NPB2, int NTD2, int NTD3, int NPH3, int NTG3, int N12, int N13, int N23, 
                     double ** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC,
                     double* segP12, double* segT12, double* segP13, double* segT13, double* segP23, double* segT23,
                     double* segVa12, double* segEa12, double* segVb12, double* segEb12,
                     double* segVa13, double* segEa13, double* segVc13, double* segEc13,
                     double* segVb23, double* segEb23, double* segVc23, double* segEc23){

  int err;
  int Nmax=1e4;
  double epsilon=5e-5;

  double Ptriple=4.5104217816091e9;
  double Pdeb=Ptriple;

  int phaseA12=1;
  int phaseB12=2;
  double Pfin12=14.62153e9; //9.4e9;
  
  int phaseA13=1;
  int phaseB13=3;
  double Pfin13=-3e9;   // 0
  
  int phaseA23=2;
  int phaseB23=3;
  double Pfin23=80e9;
  
  double Va012, Vb012, Va013, Vb013, Va023, Vb023;
  double Ea012, Eb012, Ea013, Eb013, Ea023, Eb023;

  
  // Valeur d'initialisation
  Va012=129.905671875844291e-6;  Ea012=72.168190127265447e3;
  Vb012=127.666112185063994e-6;  Eb012=99.480256119236515e3;
  
  Va013=129.905671875844291e-6;  Ea013=72.168190127265447e3;
  Vb013=132.538478805442338e-6;  Eb013=136.154346228414312e3;
  
  Va023=127.666112185063994e-6;  Ea023=99.480256119236515e3;
  Vb023=132.538478805442338e-6;  Eb023=136.154346228414312e3;


  // frontiere 12
  //printf(" ligne 12  N12= %d\n",N12);
  ligne_double(Nmax, epsilon, phaseA12, phaseB12, Pdeb, Pfin12, N12, Va012, Ea012, Vb012, Eb012, segP12, segT12, segVa12, segEa12, segVb12, segEb12);
  
  // frontiere 13
  //printf(" ligne 13   N13= %d\n",N13);
  ligne_double(Nmax, epsilon, phaseA13, phaseB13, Pdeb, Pfin13, N13, Va013, Ea013, Vb013, Eb013, segP13, segT13, segVa13, segEa13, segVc13, segEc13);
  
  // frontiere 23
  //printf(" ligne 23   N23= %d\n",N23);
  ligne_double(Nmax, epsilon, phaseA23, phaseB23, Pdeb, Pfin23, N23, Va023, Ea023, Vb023, Eb023, segP23, segT23, segVb23, segEb23, segVc23, segEc23);
  
  
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
  double TrefhautC=25e3;
  double TrefDC=segT23[N23-1];
  double TrefGC=segT13[N13-1];
  double PrefGBC=segP13[N13-1];
  double PrefGHC=71e9;
  double PrefDHC=170e9;
  double PrefDC=segP23[N23-1];
  err=sommetsC(N13, N23, NTD3, NPH3, NTG3, PrefGBC,PrefGHC,PrefDHC, PrefDC, TrefGC, TrefDC, TrefhautC, segVc13, segEc13, segVc23, segEc23, SC);
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
  int Nmax=1e2;
  double epsilon=1e-8;

  double E,V;

  // nb_points est un argument de la fonction
  int nb_points=10;
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3; 
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;
  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );
  
  
  double** SA=alloctabd(NA+1,2);
  double** SB=alloctabd(NB+1,2);
  double** SC=alloctabd(NC+1,2); // on ajoute S[N]=S[0]  --> utile pou l'algo
  double** SAB=alloctabd(NAB+1,2);
  double** SAC=alloctabd(NAC+1,2);
  double** SBC=alloctabd(NBC+1,2);

  // Initialisation des tableaux
  double* segP12=malloc(N12*sizeof(double));   double* segT12=malloc(N12*sizeof(double));
  double* segVa12=malloc(N12*sizeof(double));  double* segEa12=malloc(N12*sizeof(double));
  double* segVb12=malloc(N12*sizeof(double));  double* segEb12=malloc(N12*sizeof(double));
  
  double* segP13=malloc(N13*sizeof(double));   double* segT13=malloc(N13*sizeof(double));
  double* segVa13=malloc(N13*sizeof(double));  double* segEa13=malloc(N13*sizeof(double));
  double* segVc13=malloc(N13*sizeof(double));  double* segEc13=malloc(N13*sizeof(double));
  
  double* segP23=malloc(N23*sizeof(double));   double* segT23=malloc(N23*sizeof(double));
  double* segVb23=malloc(N23*sizeof(double));  double* segEb23=malloc(N23*sizeof(double));
  double* segVc23=malloc(N23*sizeof(double));  double* segEc23=malloc(N23*sizeof(double));
  

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
   


  // DIagramme (V,E) 
  // Ecriture des sommets des poly pour vérif
  FILE *fSA, *fSB, *fSC, *fSAB, *fSAC, *fSBC;
  if((fSA = fopen("fichiers/pSA.txt", "w+")) == NULL){printf("erreur ouverture fichier pSA\n");return 1;}
  if((fSB = fopen("fichiers/pSB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSB\n");return 1;}
  if((fSC = fopen("fichiers/pSC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSC\n");return 1;}
  if((fSAB = fopen("fichiers/pSAB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAB\n");return 1;}
  if((fSAC = fopen("fichiers/pSAC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAC\n");return 1;}
  if((fSBC = fopen("fichiers/pSBC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSBC\n");return 1;}
  

  for(int i=0; i<NA+1; i++){  fprintf(fSA,"%.10lf %.10lf\n",SA[i][0],SA[i][1] );  }
  for(int i=0; i<NB+1; i++){  fprintf(fSB,"%.10lf %.10lf\n",SB[i][0],SB[i][1] );  }
  for(int i=0; i<NC+1; i++){  fprintf(fSC,"%.10lf %.10lf\n",SC[i][0],SC[i][1] );  }
  for(int i=0; i<NAB+1; i++){  fprintf(fSAB,"%.10lf %.10lf\n",SAB[i][0],SAB[i][1] );  }
  for(int i=0; i<NAC+1; i++){  fprintf(fSAC,"%.10lf %.10lf\n",SAC[i][0],SAC[i][1] );  }
  for(int i=0; i<NBC+1; i++){  fprintf(fSBC,"%.10lf %.10lf\n",SBC[i][0],SBC[i][1] );  }
  
  
  fclose(fSA);   fclose(fSB);   fclose(fSC);
  fclose(fSAB);  fclose(fSAC);  fclose(fSBC);


  // DIagramme (P,T)
  // Ecriture des sommets des poly pour vérif
  //FILE *fSA, *fSB, *fSC, *fSAB, *fSAC, *fSBC;
  if((fSA = fopen("fichiers/pSA_PT.txt", "w+")) == NULL){printf("erreur ouverture fichier pSA\n");return 1;}
  if((fSB = fopen("fichiers/pSB_PT.txt", "w+")) == NULL){printf("erreur ouverture fichier pSB\n");return 1;}
  if((fSC = fopen("fichiers/pSC_PT.txt", "w+")) == NULL){printf("erreur ouverture fichier pSC\n");return 1;}
  

  for(int i=0; i<NA+1; i++){  fprintf(fSA,"%.10lf %.10lf\n",fP(SA[i][0],SA[i][1],1 ),fT(SA[i][0],SA[i][1],1 )) ;  }
  for(int i=0; i<NB+1; i++){  fprintf(fSB,"%.10lf %.10lf\n",fP(SB[i][0],SB[i][1],2 ),fT(SB[i][0],SB[i][1],2 )) ;  }
  for(int i=0; i<NC+1; i++){  fprintf(fSC,"%.10lf %.10lf\n",fP(SC[i][0],SC[i][1],3 ),fT(SC[i][0],SC[i][1],3 )) ;  }
  
  
  fclose(fSA);   fclose(fSB);   fclose(fSC);
  
  


  //  ****  Permet uen répresentation des zones (calcul la zone de chaque points)  **********
  /*
  double Vmin, Vmax;
  double Emin, Emax;

  int Nv=1e3;
  int Ne=1e3;

  Vmin=115e-6, Vmax=1450e-6;
  Emin=-5e3,   Emax=4e5;
  
  double he=(Emax-Emin)/(Ne-1);
  double hv=(Vmax-Vmin)/(Nv-1);
  
  
  FILE* fZONE; FILE* fE; FILE* fV;
  if((fZONE = fopen("ZONE.txt", "w+")) == NULL){printf("erreur ouverture fichier ZONE\n"); return 1;}
  if((fV = fopen("Vzone.txt", "w+")) == NULL){printf("erreur ouverture fichier Vzone \n"); return 1;}
  if((fE = fopen("Ezone.txt", "w+")) == NULL){printf("erreur ouverture fichier Ezone \n"); return 1;}
  
  for(int j=0; j<Ne; j++){
      E=Emin+j*he;
      fprintf(fE, "%.15lf ",E);
  }
  
  
  for(int i=0; i<Nv; i++){
    printf(" -- i= %d\n",i );
    V=Vmin+i*hv;
    fprintf(fV, "%.15lf ", V);
    for(int j=0; j<Ne; j++){
      //printf("  -j= %d\n",j );
      E=Emin+j*he;
      //printf("V= %.15lf, E=%.15lf\n",V*1e6,E );
      if(ftri(V, E)){fprintf(fZONE, "0 ");}
      else if(intPOLY(V,E,NA,SA)){fprintf(fZONE, "1 ");}
      else if(intPOLY(V,E,NB,SB)){fprintf(fZONE, "2 ");}
      else if(intPOLY(V,E,NC,SC)){fprintf(fZONE, "3 ");}
      else if(intPOLY(V,E,NAB,SAB)){fprintf(fZONE, "12 ");}
      else if(intPOLY(V,E,NAC,SAC)){fprintf(fZONE, "13 ");}
      else if(intPOLY(V,E,NBC,SBC)){fprintf(fZONE, "23 ");}
      else{fprintf(fZONE, "-1 ");}
    }
    fprintf(fZONE, "\n");
  }  

  fclose(fZONE);
  fclose(fV);
  fclose(fE);
  */


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
  

}


/* Valeur d'entree :
    - Va,Ea, Vb,Eb
    - Vmel, Emel

    - A,InvA,dX : pour limiter les allocations 


  Valeur sortie
    - Vx,Ex, Vy,Ey, 
    - x : fraction massique de la phase Y

   double* A et InvA permet une seule allocation 

*/
void Newton_mixte(double* A, double* InvA, double* Delta, int phaseX, int phaseY, double Vmel, double Emel, double Vx, double Ex, double Vy, double Ey, double x, double* pVX, double* pEX, double* pVY, double* pEY, double* pX, double* dX){
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
  
  
  Delta[0]=Py-Px;   Delta[1]=Ty-Tx;   Delta[2]=Gy-Gx;
  Delta[3]=Vmel - (1.-x)*Vx - x*Vy;
  Delta[4]=Emel - (1.-x)*Ex - x*Ey;
   
  /*
  for(int j=0; j<dim; j++){
    printf("Delta[%d]= %.g  ",j, Delta[j] );
  }
  printf("\n");
  */
  
  
  // on met le second menbre dans la dernière colonne
  //matrice A
  A[dim*0+0]=dPvx; A[dim*0+1]=dPex; A[dim*0+2]=-dPvy; A[dim*0+3]=-dPey; A[dim*0+4]=0.;  
  
  A[dim*1+0]=dTvx; A[dim*1+1]=dTex; A[dim*1+2]=-dTvy; A[dim*1+3]=-dTey; A[dim*1+4]=0.; 
  
  A[dim*2+0]=dGvx; A[dim*2+1]=dGex; A[dim*2+2]=-dGvy; A[dim*2+3]=-dGey; A[dim*2+4]=0.; 
  
  A[dim*3+0]=1.-x; A[dim*3+1]=0.;   A[dim*3+2]=x;     A[dim*3+3]=0;     A[dim*3+4]=Vy-Vx; 
  
  A[dim*4+0]=0.;   A[dim*4+1]=1.-x; A[dim*4+2]=0;     A[dim*4+3]=x;     A[dim*4+4]=Ey-Ex; 
  
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
  //printf("  * InvA*A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      res=0;
      for(int k=0; k<dim; k++){
        res+=InvA[dim*i+k]*A[dim*k+j];
      }
      //printf("A*InvA[%d][%d]= %.15lf ",i,j,res );
    }
    //printf("\n");
  }
  //printf("\n");
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
int VE_mixte(int Nmax, double epsilon, int phaseX, int phaseY, double V0x, double E0x, double V0y, double E0y, double Vmel, double Emel, double* pVx, double* pEx, double* pVy, double* pEy, double* px){
  int n=1;
  int dim=5;

  double* A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));
  double* dX=malloc(dim*sizeof(double));
  double* Delta=malloc(dim*sizeof(double));
  
  double Px, Tx, Sx, Gx;
  double Py, Ty, Sy, Gy;

  double Vx, Ex, Vy, Ey, x;
  double VX, EX, VY, EY, X;


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
    
    Newton_mixte(A, InvA, Delta, phaseX, phaseY, Vmel, Emel, Vx, Ex, Vy, Ey, x, &VX, &EX, &VY, &EY, &X, dX);
    
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
    
    /*
    // phase X
    Px=fP(Vx,Ex, phaseX); Tx=fT(Vx,Ex, phaseX);
    Sx=fS(Vx,Ex, phaseX); Gx=Ex+Px*Vx-Tx*Sx;
    // phase Y
    Py=fP(Vy,Ey, phaseY); Ty=fT(Vy,Ey, phaseY);
    Sy=fS(Vy,Ey, phaseY); Gy=Ey+Py*Vy-Ty*Sy;
    */

    //critere=fabs(Gx-Gy)+fabs(Px-Py)+fabs(Ty-Tx);
    critere=fabs(dX[0]) +fabs(dX[1]) +fabs(dX[2]) +fabs(dX[3])+fabs(dX[4]);
    //critere=fabs(dX[0]) +fabs(dX[1]) +fabs(dX[2]) +fabs(dX[3])+fabs(dX[4]);
    //printf("    critere= %.14lf\n",critere );

    
    n++;
  }
  //printf(" --n= %d, critere= %g\n",n,critere );
  //for(int i=0; i<dim; i++){  printf(" dX[%d]= %g  ",i,dX[i] );  } printf("\n");
  
  if(n==Nmax){printf("Newton_mixte non convergent : -n= %d , critere= %g\n",n,critere );return 1;}

  *pVx=Vx; *pEx=Ex;
  *pVy=Vy; *pEy=Ey;
  *px=x;

  /*
  if(x<0 || 1<x){
    printf("Erreur fraction massique= %g (VE_mixte)\n",x );
    printf("Newton_mixte : -n= %d , critere= %lf \n",n,critere );
    return 1;
  }
  */

  //printf("  ok fin boucle i\n");
  free(A); 
  free(InvA);
  free(dX);
  free(Delta);
  
  return 0;
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


  // nb_points est l'argument d'entree de la fonction
  int nb_points=20;
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3; 
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;
  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );
  
  
  double** SA=alloctabd(NA+1,2);
  double** SB=alloctabd(NB+1,2);
  double** SC=alloctabd(NC+1,2); // on ajoute S[N]=S[0]  --> utile pour l'algo
  double** SAB=alloctabd(NAB+1,2);
  double** SAC=alloctabd(NAC+1,2);
  double** SBC=alloctabd(NBC+1,2);

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
  

  // CALCUL DES SOMMETS  
  err=sommets_polygone(NA, NB, NC, NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3, N12, N13, N23, SA, SB, SC, SAB, SAC, SBC, segP12, segT12, segP13, segT13, segP23, segT23, 
                       segVa12, segEa12, segVb12, segEb12, segVa13, segEa13, segVc13, segEc13, segVb23, segEb23, segVc23, segEc23 );

  if(err){printf("err= %d (sommets_polygone)\n",err ); return err;}



  ///////////
  // Ecriture des sommets des poly pour vérif
  FILE *fSA, *fSB, *fSC, *fSAB, *fSAC, *fSBC;
  if((fSA = fopen("fichiers/pSA.txt", "w+")) == NULL){printf("erreur ouverture fichier pSA\n");return 1;}
  if((fSB = fopen("fichiers/pSB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSB\n");return 1;}
  if((fSC = fopen("fichiers/pSC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSC\n");return 1;}
  if((fSAB = fopen("fichiers/pSAB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAB\n");return 1;}
  if((fSAC = fopen("fichiers/pSAC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAC\n");return 1;}
  if((fSBC = fopen("fichiers/pSBC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSBC\n");return 1;}
  
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
  if((fmaps = fopen("fichiers/MAPS.txt", "w+")) == NULL){printf("erreur ouverture fichier MAPS\n");return 1;}
  if((fmapsGNU = fopen("fichiers/P_GNU.txt", "w+")) == NULL){printf("erreur ouverture fichier MAPS\n");return 1;}
  if((fV = fopen("fichiers/Vmaps.txt", "w+")) == NULL){printf("erreur ouverture fichier Vmaps \n");return 1;}
  if((fE = fopen("fichiers/Emaps.txt", "w+")) == NULL){printf("erreur ouverture fichier Emaps \n");return 1;}


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
      if(intPOLY(V,E,NC,SC)){ 
        P=fP(V,E,3);
        T=fT(V,E,3);        
        fprintf(fmaps, "%.15lf %.15lf\n", P, T);
        fprintf(fmapsGNU, "%.15lf %.15lf %.15lf\n", V, E, P);
      }
      else if(intPOLY(V,E,NA,SA)){
        P=fP(V,E,1);
        T=fT(V,E,1);
        fprintf(fmaps, "%.15lf %.15lf\n", P, T);
        fprintf(fmapsGNU, "%.15lf %.15lf %.15lf\n", V, E, P);
      }
      else if(intPOLY(V,E,NB,SB)){
        P=fP(V,E,2);
        T=fT(V,E,2);
        fprintf(fmaps, "%.15lf %.15lf\n", P, T);
        fprintf(fmapsGNU, "%.15lf %.15lf %.15lf\n", V, E, P);
      }
      else if(ftri(V, E)){
        P=PT;
        T=TT;
        fprintf(fmaps, "%.15lf %.15lf\n", P, T);
        fprintf(fmapsGNU, "%.15lf %.15lf %.15lf\n", V, E, P);
      }
      else if(intPOLY(V,E,NAB,SAB)){
        /*
        err=VE_mixte(Nmax, epsilon, 1, 2, V, E, &Vx, &Ex, &Vy, &Ey, &x);
        P=fP(Vx,Ex,1);
        T=fT(Vx,Ex,1);
        */
        err=interp_zone_mixte(V, E, 1, 2, N12, segVa12, segVb12, segEa12, segEb12, &P, &T, &V_unuse_x, &E_unuse_x, &V_unuse_y, &E_unuse_y);
        //printf("P= %lf, T= %lf   mixte 1/2\n",P*1e-9,T );
        fprintf(fmaps, "%.15lf %.15lf\n", P,T);
        fprintf(fmapsGNU, "%.15lf %.15lf %.15lf\n", V, E, P);
      }
      else if(intPOLY(V,E,NAC,SAC)){
        /*
        err=VE_mixte(Nmax, epsilon, 1, 3, V, E, &Vx, &Ex, &Vy, &Ey, &x);
        P=fP(Vx,Ex,1);
        T=fT(Vx,Ex,1);
        */
        err=interp_zone_mixte(V, E, 1, 3, N13, segVa13, segVc13, segEa13, segEc13, &P, &T, &V_unuse_x, &E_unuse_x, &V_unuse_y, &E_unuse_y);
        //printf("P= %lf, T= %lf   mixte 1/3\n",P*1e-9,T );
        fprintf(fmaps, "%.15lf %.15lf\n", P,T);
        fprintf(fmapsGNU, "%.15lf %.15lf %.15lf\n", V, E, P);
      }
      else if(intPOLY(V,E,NBC,SBC)){
        /*
        err=VE_mixte(Nmax, epsilon, 2, 3, V, E, &Vx, &Ex, &Vy, &Ey, &x);
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



// Calcul de la célérite du son
/*  Pp=P(V+dV,E+dE(+dV) )  
    Pm=P(V-dV,E+dE(-dV) )  

    Argument de sortie :
     - double* pc : pointeur sur la celerité du son
*/
int c_celerite_son(int Nmax, double epsilon, int zone, double P, double V, double E, double dV, 
                   double V0x, double E0x, double V0y, double E0y, double *pc, double *pG){

  int err;

  double Pp,Pm;
  double dE=-P*dV;

  double Vx,Ex,Vy,Ey,x;

  int phaseA,phaseB;

  if(zone==-1){printf("Erreur out of domain : zone= %d (c_celerite_son)\n",zone); return 1;  }
  else if(zone==1 || zone==2 || zone==3){  
    Pp=fP(V+dV,E+dE,zone); 
    Pm=fP(V-dV,E-dE,zone); 
  }
  else if(zone==12 || zone==13 || zone==23){
    phaseB=zone%10;
    phaseA=(zone-phaseB)/10;
    err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V+dV, E+dE, &Vx, &Ex, &Vy, &Ey, &x);  if(err){return 1;}
    Pp=fP(Vx,Ex,phaseA);

    err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V-dV, E-dE, &Vx, &Ex, &Vy, &Ey, &x);  if(err){return 1;}
    Pm=fP(Vx,Ex,phaseA);  
  }
  else if(zone==0){
    *pc=0;
    return 0;
  }

  double dPpdV=(Pp-P)/dV;
  double dPmdV=(Pm-P)/(-dV);
  double dPdV=(dPpdV + dPmdV)/2;
  double c2=-V*V*dPdV;  

  double d2PdV2=(Pm+Pp-2*P)/(dV*dV);
  *pG=-(V/2)*d2PdV2/dPdV;
  
  if(c2<0){
    printf("Erreur c2= %g negatif\n",c2);
    printf("  V= %g , E= %g\n",V*1e6,E );  
    printf("  dV/V= %g , dE/E= %g \n",dV/V,dE/E );
    printf("  P= %g , Pp= %g , Pm= %g\n",P,Pp,Pm );
    printf("  dPp/dV= %g , dPm/dV= %g , dPdV= %g \n",dPpdV,dPmdV,dPdV );
    //*pc=3.3e2;
    return 1;
  }
  *pc=sqrt(c2);
  return 0;
}


// Fonction pour avoir une benchmark de (P,T) sur (V,E)
/*

*/
int ecr_benchmark_PT(int nb_points, int Nv ,int Ne, double Vmin, double Vmax, double Emin, double Emax){

  // Ecriture dans un fichier texte des paramètres de la discretisation
  // pour pouvoir utiliser la discretisations juste en lisant les fichiers
  FILE *fparams;
  if((fparams = fopen("fichiers/PARAMS_dis.txt", "w+")) == NULL){printf("erreur ouverture fichier PARAMS_dis.txt \n"); return 1;}
  fprintf(fparams, "%d \n",nb_points );
  fprintf(fparams, "%d %d\n",Nv,Ne );
  fprintf(fparams, "%.15lf %.15lf\n",Vmin,Vmax );
  fprintf(fparams, "%.15lf %.15lf\n",Emin,Emax );
  fprintf(fparams, "\n\n Organisation du fichier\n");
  fprintf(fparams, "nb_points\n");
  fprintf(fparams, "Nv_dis  Ne_dis\n");
  fprintf(fparams, "Vmin_dis Vmax_dis\n");
  fprintf(fparams, "Emin_dis Emax_dis\n");


  double he=(Emax-Emin)/(Ne-1);
  double hv=(Vmax-Vmin)/(Nv-1);


  int Nmax=1e4;
  double epsilon=1e-6;
  double Vx, Ex, Vy, Ey, x;
  

  int err,res;
  

  double E,V;

  double Vmx, Emx;
  double Vmy, Emy;


  // Varaibles pour discretisation de la maille
  double V_sub, E_sub;
  int N_rec=10;
  double he_sub=he/(N_rec-1);
  double hv_sub=hv/(N_rec-1);
  double** VE_front_rec=alloctabd(4*(N_rec-1),2);
  int count;
  

  // nb_points est un argument de la fonction
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3; 
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;
  
  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );
  

  double** SA=alloctabd(NA+1,2);
  double** SB=alloctabd(NB+1,2);
  double** SC=alloctabd(NC+1,2); // on ajoute S[N]=S[0]  --> utile pour l'algo
  double** SAB=alloctabd(NAB+1,2);
  double** SAC=alloctabd(NAC+1,2);
  double** SBC=alloctabd(NBC+1,2);
  
   
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
  


  // CALCUL DES SOMMETS  
  printf(" Sommets_polygone \n");
  err=sommets_polygone(NA, NB, NC, NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3, N12, N13, N23, SA, SB, SC, SAB, SAC, SBC, segP12, segT12, segP13, segT13, segP23, segT23, 
                       segVa12, segEa12, segVb12, segEb12, segVa13, segEa13, segVc13, segEc13, segVb23, segEb23, segVc23, segEc23 );
  if(err){printf("err= %d (sommets_polygone)\n",err ); return err;}


  double** STRI_12=alloctabd(2,2); 
  STRI_12[0][0]=129.905671875844291e-6;   STRI_12[0][1]=72.168190127265447e3;
  STRI_12[1][0]=127.666112185063994e-6;   STRI_12[1][1]=99.480256119236515e3;
  
  double** STRI_13=alloctabd(2,2); 
  STRI_13[0][0]=129.905671875844291e-6;   STRI_13[0][1]=72.168190127265447e3;
  STRI_13[1][0]=132.538478805442338e-6;   STRI_13[1][1]=136.154346228414312e3;
  
  double** STRI_23=alloctabd(2,2); 
  STRI_23[0][0]=127.666112185063994e-6;   STRI_23[0][1]=99.480256119236515e3;
  STRI_23[1][0]=132.538478805442338e-6;   STRI_23[1][1]=136.154346228414312e3;
  
  // Coordonnées de la ligne Sa12  (sépartion de A et de la zone mixte 12 ou A/B)
  //  et de Sb12
  double** Sa12=alloctabd(N12,2);  
  double** Sb12=alloctabd(N12,2); 
  for(int i=0; i<N12; i++){
    Sa12[i][0]=segVa12[i]; 
    Sa12[i][1]=segEa12[i];

    Sb12[i][0]=segVb12[i]; 
    Sb12[i][1]=segEb12[i];  
  }

  // Coordonnées de la ligne Sa13  (sépartion de A et de la zone mixte 13 ou A/C)
  //  et de Sb12
  double** Sa13=alloctabd(N13,2);  
  double** Sc13=alloctabd(N13,2); 
  for(int i=0; i<N13; i++){
    Sa13[i][0]=segVa13[i];
    Sa13[i][1]=segEa13[i];

    Sc13[i][0]=segVc13[i]; 
    Sc13[i][1]=segEc13[i];  
  }

  // Coordonnées de la ligne Sa13  (sépartion de A et de la zone mixte 13 ou A/C)
  //  et de Sb12
  double** Sb23=alloctabd(N23,2);  
  double** Sc23=alloctabd(N23,2); 
  for(int i=0; i<N23; i++){
    Sb23[i][0]=segVb23[i]; 
    Sb23[i][1]=segEb23[i];

    Sc23[i][0]=segVc23[i]; 
    Sc23[i][1]=segEc23[i];  
  }

  // frontières extéieures des phases
  int NextA=NA-N12-(N13-1);
  double** SextA=alloctabd(NextA,2);  
  for(int i=0; i<NextA; i++){
    SextA[i][0]=SA[N13+i][0]; 
    SextA[i][1]=SA[N13+i][1]; 
  }

  int NextB=NB-N12-(N23-1);
  double** SextB=alloctabd(NextB,2);  
  for(int i=0; i<NextB; i++){
    SextB[i][0]=SB[N12+i][0]; 
    SextB[i][1]=SB[N12+i][1]; 
  }

  int NextC=NC-N23-(N13-1);
  double** SextC=alloctabd(NextC,2);  
  for(int i=0; i<NextC; i++){
    SextC[i][0]=SC[N23+i][0]; 
    SextC[i][1]=SC[N23+i][1]; 
  }
  
  
  // Pour tester les valeurs (coord des points du maillage)
  FILE *fileV, *fileE;
  if((fileV = fopen("fichiers/Vbench.txt", "w+")) == NULL){printf("erreur ouverture fichier V bench \n"); return 1;}
  if((fileE = fopen("fichiers/Ebench.txt", "w+")) == NULL){printf("erreur ouverture fichier E bench \n"); return 1;}
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
  

  double* tabVx=malloc((Nv-1)*(Ne-1)*sizeof(double)); double* tabEx=malloc((Nv-1)*(Ne-1)*sizeof(double));
  double* tabVy=malloc((Nv-1)*(Ne-1)*sizeof(double)); double* tabEy=malloc((Nv-1)*(Ne-1)*sizeof(double));
  
  int* tabNbZONE=malloc((Nv-1)*(Ne-1)*sizeof(int)); // nombre de zones présentent sur le triangle
  int** tabZONE=alloctab(5,(Nv-1)*(Ne-1));  // mémorise les zones présentent dans le rectangle
  int** tabINTER=alloctab(9,(Nv-1)*(Ne-1)); // mémorise les lignes avec lesquelles il y a une intersection
  int** tabInd_Mixte=alloctab(9,(Nv-1)*(Ne-1)); // indice du segment de la ligne de chgmt de phase qui s'intersecte avec le rectangle
  

  // une ligne un carré puis sur les colonnes, les zones qui le composent on rempli par -99 si rien dans les zones 
  // (ou 0 peut-être plus simple à implementer mais il faudra changer l'indice du tri) 

  int** numZONE=alloctab(234,2);
  numZONE[121][0]=1;  numZONE[121][1]=12;
  numZONE[122][0]=2;  numZONE[122][1]=12;
  numZONE[131][0]=1;  numZONE[131][1]=13;
  numZONE[133][0]=3;  numZONE[133][1]=13;
  numZONE[232][0]=2;  numZONE[232][1]=23;
  numZONE[233][0]=3;  numZONE[233][1]=23;
  numZONE[120][0]=0;  numZONE[120][1]=12;
  numZONE[130][0]=0;  numZONE[130][1]=13;
  numZONE[230][0]=0;  numZONE[230][1]=23;
  numZONE[11][0]=1;   numZONE[11][1]=-1;
  numZONE[12][0]=2;   numZONE[12][1]=-1;
  numZONE[13][0]=3;   numZONE[13][1]=-1;
 


  double** REC=alloctabd(5,2); // rectangle du cadrillage
  
  int ind_mixte;
  int nb_zone;
  int nb_inter=0;
  int ind_inter;
  int test;

  int k=0;
  for(int i=0; i<Nv-1; i++){
    printf("  --i= %d\n",i );
    V=Vmin+i*hv + hv/2; // pour prendre la valeur au centre du rectangle
    for(int j=0; j<Ne-1; j++){
      //printf("  -- i= %d, j= %d\n",i,j );

      E=Emin+j*he + he/2;  
      
      nb_inter=0;

      //   Calcul du chevanchement des zones
      REC[0][0]=Vmin+i*hv;       REC[0][1]=Emin+j*he;
      REC[1][0]=Vmin+(i+1)*hv;   REC[1][1]=Emin+j*he;
      REC[2][0]=Vmin+(i+1)*hv;   REC[2][1]=Emin+(j+1)*he;
      REC[3][0]=Vmin+i*hv;       REC[3][1]=Emin+(j+1)*he;
      REC[4][0]=REC[0][0];       REC[4][1]=REC[0][1];

      
      res=inter_seg(REC, N12-1, Sa12, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=121; tabInd_Mixte[nb_inter][k]=ind_mixte; nb_inter++;}
      //else if(res==0){tabZONE[1][k]=-99}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}
      
      res=inter_seg(REC, N12-1, Sb12, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=122; tabInd_Mixte[nb_inter][k]=ind_mixte; nb_inter++;}
      //else if(res==0){tabZONE[1][k]=-99}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}

      res=inter_seg(REC, N13-1, Sa13, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=131; tabInd_Mixte[nb_inter][k]=ind_mixte; nb_inter++;}
      //else if(res==0){tabZONE[1][k]=-99}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}

      res=inter_seg(REC, N13-1, Sc13, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=133; tabInd_Mixte[nb_inter][k]=ind_mixte; nb_inter++;}
      //else if(res==0){tabZONE[1][k]=-99}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}

      res=inter_seg(REC, N23-1, Sb23, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=232; tabInd_Mixte[nb_inter][k]=ind_mixte; nb_inter++;}
      //else if(res==0){tabZONE[1][k]=-99}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}

      res=inter_seg(REC, N23-1, Sc23, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=233; tabInd_Mixte[nb_inter][k]=ind_mixte; nb_inter++;}
      //else if(res==0){tabZONE[1][k]=-99}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}
      
      // Pour les frontières triangles, on mets ind_mixte=0 pour avoir les pointd triples comme graine
      res=inter_seg(REC, 1, STRI_12, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=120; ind_mixte=0; nb_inter++;}
      //else if(res==0){tabZONE[1][k]=-99}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}

      res=inter_seg(REC, 1, STRI_13, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=130; ind_mixte=0; nb_inter++;}
      //else if(res==0){tabZONE[1][k]=-99}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}

      res=inter_seg(REC, 1, STRI_23, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=230; ind_mixte=0; nb_inter++;}
      //else if(res==0){tabZONE[1][k]=-99}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}

      res=inter_seg(REC, NextA-1, SextA, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=11; ind_mixte=0; nb_inter++;}
      //else if(res==0){tabZONE[1][k]=-99}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}

      res=inter_seg(REC, NextB-1, SextB, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=12; ind_mixte=0; nb_inter++;}
      //else if(res==0){tabZONE[1][k]=-99}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}

      res=inter_seg(REC, NextC-1, SextC, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=13; ind_mixte=0; nb_inter++;}
      //else if(res==0){tabZONE[1][k]=-99}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}
      

  
      //printf("   Nb intersection= %d\n",nb_inter );
      //for(int i=0; i<nb_inter; i++){
      //  printf("  %d : %d\n",i,tabINTER[i][k] );
      //}

      
      // Cas de rectangle sans chevauchement
      if(nb_inter==0){
        tabNbZONE[k]=1;
        // Attribution d'une zone aux points milieux des rectangles
        if(intPOLY(V,E,NC,SC)){ 
          tabZONE[0][k]=3;
        }
        else if(intPOLY(V,E,NA,SA)){
          tabZONE[0][k]=1;
        }
        else if(intPOLY(V,E,NB,SB)){
          tabZONE[0][k]=2;
        }
        else if(ftri(V,E)){
          tabZONE[0][k]=0;
        }
        else if(intPOLY(V,E,NAB,SAB)){
          err=graine_zone_mixte(V, E, N12, segVa12, segVb12, segEa12, segEb12, &Vmx, &Emx, &Vmy, &Emy);
          if(err){printf("err= %d (benchmark_PT)\n",err ); return err;}
          tabVx[k]=Vmx; tabEx[k]=Emx;
          tabVy[k]=Vmy; tabEy[k]=Emy;
          tabZONE[0][k]=12;
        }
        else if(intPOLY(V,E,NAC,SAC)){
          err=graine_zone_mixte(V, E, N13, segVa13, segVc13, segEa13, segEc13, &Vmx, &Emx, &Vmy, &Emy);
          if(err){printf("err= %d (benchmark_PT)\n",err ); return err;}
          tabVx[k]=Vmx; tabEx[k]=Emx;
          tabVy[k]=Vmy; tabEy[k]=Emy;
          tabZONE[0][k]=13;
        }
        else if(intPOLY(V,E,NBC,SBC)){
          err=graine_zone_mixte(V, E, N23, segVb23, segVc23, segEb23, segEc23, &Vmx, &Emx, &Vmy, &Emy);
          if(err){printf("err= %d (benchmark_PT)\n",err ); return err;}
          tabVx[k]=Vmx; tabEx[k]=Emx;
          tabVy[k]=Vmy; tabEy[k]=Emy;
          tabZONE[0][k]=23;
        }
        else{  // hors domaine --> on ne se soucis des pb de chevauchement des rectangles sur les zones
          tabZONE[0][k]=-1;
        }
      }

      // ## 1 intersection donc seulment deux zones ##
      else if(nb_inter==1){
        tabNbZONE[k]=2;
        tabZONE[0][k]=numZONE[ tabINTER[nb_inter-1][k] ][0];
        tabZONE[1][k]=numZONE[ tabINTER[nb_inter-1][k] ][1];
        //printf("tabINTER[%d][%d]= %d,  tabZONE[0][%d]= %d, tabZONE[1][%d]= %d \n",nb_inter-1,k, tabINTER[nb_inter-1][k],k,tabZONE[0][k],k,tabZONE[1][k]);
        // On met les graines pour le Newton mixte
        if(tabZONE[1][k]==12){
          tabVx[k]=segVa12[ind_mixte];
          tabEx[k]=segEa12[ind_mixte];
          tabVy[k]=segVb12[ind_mixte];
          tabEy[k]=segEb12[ind_mixte];
        }
        else if(tabZONE[1][k]==13){
          tabVx[k]=segVa13[ind_mixte];
          tabEx[k]=segEa13[ind_mixte];
          tabVy[k]=segVc13[ind_mixte];
          tabEy[k]=segEc13[ind_mixte];
        }
        else if(tabZONE[1][k]==23){
          tabVx[k]=segVb23[ind_mixte];
          tabEx[k]=segEb23[ind_mixte];
          tabVy[k]=segVc23[ind_mixte];
          tabEy[k]=segEc23[ind_mixte];
        }
        else if(tabZONE[1][k]==-1){}
        else{  printf("Erreur dans tabZONE[1]= %d (ecr_benchmark_PT)\n",tabZONE[1][k]); return 1; }
      }
      // ## Plusieurs intersections ##
      else if(nb_inter==2 || nb_inter==3  || nb_inter==4 || nb_inter==5){
        
        //printf("  !!!  zone multi-multi-phase  nb_inter= %d !!!!\n",nb_inter);
        //for(int i=0; i<nb_inter; i++){
        //  printf("  %d : %d\n",i,tabINTER[i][k] );
        //}
        

        nb_zone=0;
        // On cree un tableau de la discretisation des fronitères du rectangle
        //printf(" Rectangle :\n");
        //printf(" Vmin= %.10lf  Vmax= %.10lf , Emin= %.10lf, Emax= %.10lf\n",(Vmin+i*hv)*1e6,(Vmin+(i+1)*hv)*1e6,Emin+j*he,Emin+(j+1)*he );
        count=0;
        for(int iv=0; iv<2; iv++){
          for(int ie=0; ie<N_rec; ie++){
            VE_front_rec[count][0]=Vmin+(i+iv)*hv;
            VE_front_rec[count][1]=Emin+j*he+ie*he_sub;
            count++;
          }
        }
        for(int ie=0; ie<2; ie++){
          for(int iv=1; iv<N_rec-1; iv++){
            VE_front_rec[count][0]=Vmin+i*hv+iv*hv_sub;
            VE_front_rec[count][1]=Emin+(j+ie)*he;
            count++;
          }
        }
        // On parcours les frontières du rectangle
        for(int ir=0; ir<4*(N_rec-1); ir++){
          //printf("  ir/(4*N_rec-2)= %d / %d\n", ir, 4*N_rec-2 );
          V_sub=VE_front_rec[ir][0];
          E_sub=VE_front_rec[ir][1];
          if(intPOLY(V_sub,E_sub,NC,SC)){ 
            test=1;
            for(int i=0;i<nb_zone;i++){
              if(tabZONE[i][k]==3){test=0;}
            }
            if(test==1){
              tabZONE[nb_zone][k]=3;  
              nb_zone++;
             }
          }
          else if(intPOLY(V_sub,E_sub,NA,SA)){
            test=1;
            for(int i=0;i<nb_zone;i++){
              if(tabZONE[i][k]==1){test=0;}
            }
            if(test==1){tabZONE[nb_zone][k]=1;  nb_zone++; }
          }
          else if(intPOLY(V_sub,E_sub,NB,SB)){
            test=1;
            for(int i=0;i<nb_zone;i++){
              if(tabZONE[i][k]==2){test=0;}
            }
            if(test==1){tabZONE[nb_zone][k]=2;  nb_zone++; }
          }
          else if(ftri(V_sub, E_sub)){
            test=1;
            for(int i=0;i<nb_zone;i++){
              if(tabZONE[i][k]==0){test=0;}
            }
            if(test==1){tabZONE[nb_zone][k]=0;  nb_zone++; }
          }
          else if(intPOLY(V_sub,E_sub,NAB,SAB)){
            test=1;
            for(int i=0;i<nb_zone;i++){
              if(tabZONE[i][k]==12){test=0;}
            }
            if(test==1){tabZONE[nb_zone][k]=12;  nb_zone++; }
          }
          else if(intPOLY(V_sub,E_sub,NAC,SAC)){
            test=1;
            for(int i=0;i<nb_zone;i++){
              if(tabZONE[i][k]==13){test=0;}
            }
            if(test==1){tabZONE[nb_zone][k]=13;  nb_zone++; }
          }
          else if(intPOLY(V_sub,E_sub,NBC,SBC)){
            test=1;
            for(int i=0;i<nb_zone;i++){
              if(tabZONE[i][k]==23){test=0;}
            }
            if(test==1){tabZONE[nb_zone][k]=23;  nb_zone++; }
          }
          else{
            test=1;
            for(int i=0;i<nb_zone;i++){
              if(tabZONE[i][k]==-1){test=0;}
            }
            if(test==1){tabZONE[nb_zone][k]=-1;  nb_zone++; }
          }
        }
        tabNbZONE[k]=nb_zone;
      } 
      else{printf("Erreur dans la valeur de nb_inter= %d (ecr_benchmark_PT)\n",nb_inter); return 1;}
      if (tabNbZONE[k]==0){printf("Erreur aucune zone atribuer tabNbZONE[%d]= %d  V=%.10lf , E=%.10lf\n",k,tabNbZONE[k],V,E); return 1; }
      k++;
    }
  }


  // fichier binaire   couple (V,E) pour les zones mixte
  FILE *fileVmx, *fileEmx;
  if((fileVmx = fopen("fichiers/Vbench_mixte_x.bin", "wb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmx = fopen("fichiers/Ebench_mixte_x.bin", "wb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  FILE *fileVmy, *fileEmy;
  if((fileVmy = fopen("fichiers/Vbench_mixte_y.bin", "wb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmy = fopen("fichiers/Ebench_mixte_y.bin", "wb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  
  // Ecriture des tableaux
  fwrite(tabVx,sizeof(double),(Nv-1)*(Ne-1),fileVmx);
  fwrite(tabEx,sizeof(double),(Nv-1)*(Ne-1),fileEmx);
  fwrite(tabVy,sizeof(double),(Nv-1)*(Ne-1),fileVmy);
  fwrite(tabEy,sizeof(double),(Nv-1)*(Ne-1),fileEmy);

  fclose(fileVmx);  fclose(fileEmx);
  fclose(fileVmy);  fclose(fileEmy);
  
  // fichier binaire ZONE
  FILE *fileZONE_0,*fileZONE_1,*fileZONE_2,*fileZONE_3,*fileZONE_4;
  if((fileZONE_0 = fopen("fichiers/ZONEbench_0.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE bench 0\n"); return 1;}
  if((fileZONE_1 = fopen("fichiers/ZONEbench_1.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE bench 1\n"); return 1;}
  if((fileZONE_2 = fopen("fichiers/ZONEbench_2.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE bench 2\n"); return 1;}
  if((fileZONE_3 = fopen("fichiers/ZONEbench_3.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE bench 3\n"); return 1;}
  if((fileZONE_4 = fopen("fichiers/ZONEbench_4.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE bench 4\n"); return 1;}
  
  // Ecriture des tableaux
  fwrite(tabZONE[0],sizeof(int),(Nv-1)*(Ne-1),fileZONE_0);
  fwrite(tabZONE[1],sizeof(int),(Nv-1)*(Ne-1),fileZONE_1);
  fwrite(tabZONE[2],sizeof(int),(Nv-1)*(Ne-1),fileZONE_2);
  fwrite(tabZONE[3],sizeof(int),(Nv-1)*(Ne-1),fileZONE_3);
  fwrite(tabZONE[4],sizeof(int),(Nv-1)*(Ne-1),fileZONE_4);

  fclose(fileZONE_0);  fclose(fileZONE_1);  fclose(fileZONE_2);
  fclose(fileZONE_3);  fclose(fileZONE_4);


  // fichier binaire NOMBRE DE ZONE
  FILE *file_Nb_ZONE;
  if((file_Nb_ZONE = fopen("fichiers/Nb_ZONEbench.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE bench 0\n"); return 1;}

  // Ecriture des tableaux
  fwrite(tabNbZONE,sizeof(int),(Nv-1)*(Ne-1),file_Nb_ZONE);

  fclose(file_Nb_ZONE);


  // fichier binaire POLYGONES Frontières de zones
  FILE *fileA_V,*fileB_V,*fileC_V,*fileAB_V,*fileAC_V,*fileBC_V;
  FILE *fileA_E,*fileB_E,*fileC_E,*fileAB_E,*fileAC_E,*fileBC_E;
  if((fileA_V = fopen("fichiers/A_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE A (ecr_benchmark_PT)\n"); return 1;}     if((fileA_E = fopen("fichiers/A_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE A (ecr_benchmark_PT)\n"); return 1;}
  if((fileB_V = fopen("fichiers/B_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE B (ecr_benchmark_PT)\n"); return 1;}     if((fileB_E = fopen("fichiers/B_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE B (ecr_benchmark_PT)\n"); return 1;}
  if((fileC_V = fopen("fichiers/C_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE C (ecr_benchmark_PT)\n"); return 1;}     if((fileC_E = fopen("fichiers/C_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE C (ecr_benchmark_PT)\n"); return 1;}
  if((fileAB_V = fopen("fichiers/AB_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE AB (ecr_benchmark_PT)\n"); return 1;}  if((fileAB_E = fopen("fichiers/AB_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE AB (ecr_benchmark_PT)\n"); return 1;}
  if((fileAC_V = fopen("fichiers/AC_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE AC (ecr_benchmark_PT)\n"); return 1;}  if((fileAC_E = fopen("fichiers/AC_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE AC (ecr_benchmark_PT)\n"); return 1;}
  if((fileBC_V = fopen("fichiers/BC_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE BC (ecr_benchmark_PT)\n"); return 1;}  if((fileBC_E = fopen("fichiers/BC_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE BC (ecr_benchmark_PT)\n"); return 1;}
  
  double* SA_V=malloc((NA+1)*sizeof(double));    double* SA_E=malloc((NA+1)*sizeof(double));
  double* SB_V=malloc((NB+1)*sizeof(double));    double* SB_E=malloc((NB+1)*sizeof(double));
  double* SC_V=malloc((NC+1)*sizeof(double));    double* SC_E=malloc((NC+1)*sizeof(double));
  double* SAB_V=malloc((NAB+1)*sizeof(double));  double* SAB_E=malloc((NAB+1)*sizeof(double));
  double* SAC_V=malloc((NAC+1)*sizeof(double));  double* SAC_E=malloc((NAC+1)*sizeof(double));
  double* SBC_V=malloc((NBC+1)*sizeof(double));  double* SBC_E=malloc((NBC+1)*sizeof(double));

  for(int i=0; i<NA+1; i++){  SA_V[i]=SA[i][0];  SA_E[i]=SA[i][1];  }
  for(int i=0; i<NB+1; i++){  SB_V[i]=SB[i][0];  SB_E[i]=SB[i][1];  }
  for(int i=0; i<NC+1; i++){  SC_V[i]=SC[i][0];  SC_E[i]=SC[i][1];  }

  for(int i=0; i<NAB+1; i++){  SAB_V[i]=SAB[i][0];  SAB_E[i]=SAB[i][1];  }
  for(int i=0; i<NAC+1; i++){  SAC_V[i]=SAC[i][0];  SAC_E[i]=SAC[i][1];  }
  for(int i=0; i<NBC+1; i++){  SBC_V[i]=SBC[i][0];  SBC_E[i]=SBC[i][1];  }

  // Ecriture des tableaux
  fwrite(SA_V,sizeof(double),NA+1,fileA_V);     fwrite(SA_E,sizeof(double),NA+1,fileA_E);
  fwrite(SB_V,sizeof(double),NB+1,fileB_V);     fwrite(SB_E,sizeof(double),NB+1,fileB_E);
  fwrite(SC_V,sizeof(double),NC+1,fileC_V);     fwrite(SC_E,sizeof(double),NC+1,fileC_E);
  fwrite(SAB_V,sizeof(double),NAB+1,fileAB_V);  fwrite(SAB_E,sizeof(double),NAB+1,fileAB_E);
  fwrite(SAC_V,sizeof(double),NAC+1,fileAC_V);  fwrite(SAC_E,sizeof(double),NAC+1,fileAC_E);
  fwrite(SBC_V,sizeof(double),NBC+1,fileBC_V);  fwrite(SBC_E,sizeof(double),NBC+1,fileBC_E);


  fclose(fileA_V);   fclose(fileB_V);   fclose(fileC_V);
  fclose(fileA_E);   fclose(fileB_E);   fclose(fileC_E);
  fclose(fileAB_V);  fclose(fileAC_V);  fclose(fileBC_V);
  fclose(fileAB_E);  fclose(fileAC_E);  fclose(fileBC_E);

  free(SA_V); free(SB_V); free(SC_V); free(SAB_V); free(SAC_V); free(SBC_V);
  free(SA_E); free(SB_E); free(SC_E); free(SAB_E); free(SAC_E); free(SBC_E);



  // Liberation
  free(tabVx);   free(tabEx);
  free(tabVy);   free(tabEy);
  
  freetab(tabZONE);
  freetab(REC);

  free(tabNbZONE);
  freetab(tabINTER);
  freetab(tabInd_Mixte);
  freetab(numZONE);
  freetab(VE_front_rec);


  free(segP12);   free(segT12);
  free(segVa12);  free(segEa12);
  free(segVb12);  free(segEb12);
  
  free(segP13);   free(segT13);
  free(segVa13);  free(segEa13);
  free(segVc13);  free(segEc13);
  
  free(segP23);   free(segT23);
  free(segVb23);  free(segEb23);
  free(segVc23);  free(segEc23);

  freetab(SA);   freetab(SB);   freetab(SC);
  freetab(SAB);  freetab(SAC);  freetab(SBC);
  freetab(Sa12); freetab(Sb12);
  freetab(Sa13); freetab(Sc13);
  freetab(Sb23); freetab(Sc23);

  freetab(STRI_12);  freetab(STRI_13);  freetab(STRI_23);

  fprintf(fparams, "\nEcriture bien terminé\n");
  fclose(fparams);


  return 0;

}



// Renvoie la valeur de P et T à partir du benchmark et du newton sur les zones mixtes
/*  Arguments d'entrée :
     - int affichage : affichage==1 -> sort P=-10e9 et T=0 lorsque on est hors domaine ~~ PERMET L'AFFICHAGE DE P ET T ~~
                       afichage==0 --> renvoie un erreur lorsque on est hors domaine ~~ PERMET L'UTILISATION DE LA FONCTION DANS LES SCHEMAS ~~
     - V,E
     - nb_points : nombre de points de discretisation des frontières
     - SA,SB,SC etc.. contiennent les discretisations des frontières
     - tabVx, tabEx, ect.. valeur d'initialisation pour le newton sur les zones mixtes
     - tab_Nb_ZONE : nombre de zones sur chaque case
     - tabZONE : case(s) présente(nt) sur chaque case
     - Nv,Ne, Vmin,Vmax,etc ... paramétrent du cadrillage
    Arguments de sorties :
     - pP,pT,px,pzone : pointeur sur la P,T,fraction massique et la zone de (V,E)
     - pV0x, pE0x, etc...  pointeur sur la valeur initiale choisi pour lancer une zone mixte (On peut y avoir acces à la fonction avant grace à tabVx)
*/
int lec_benchmark_PT(int affichage, double V, double E, int nb_points, double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC, 
                     double* tabVx, double* tabEx, double* tabVy, double* tabEy, int** tabZONE, int* tab_Nb_ZONE, 
                     int Nv, int Ne, double Vmin, double Vmax, double Emin, double Emax, double* pP, double* pT, double* px, int* pzone, 
                     double* pV0x, double* pE0x, double* pV0y, double* pE0y){
  int Nmax=1e2;
  double epsilon=1e-8;

  int err; 
  double P,T;
  double Vmx, Emx, Vmy, Emy;
  double Vx, Ex, Vy, Ey;

  double x;

  // nb_points est un argument de la fonction
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3; 
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;

  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );




  // Point triple
  double PT= 4510421781.609182357788085;
  double TT= 624.63982538492518;

  if(V<Vmin || Vmax<V){printf("V hors domaine : Vmin= %.15lf  Vmax= %.15lf  et V= %.15lf (lec_benchmark_PT)\n",Vmin,Vmax,V);   return 1;}  
  if(E<Emin || Emax<E){printf("E hors domaine : Emin= %.15lf  Emax= %.15lf  et E= %.15lf (lec_benchmark_PT)\n",Emin,Emax,E);   return 1;}
 
  double hv=(Vmax-Vmin)/(Nv-1);
  double he=(Emax-Emin)/(Ne-1);
  int i=(V-Vmin)/hv;
  int j=(E-Emin)/he;


 //  printf("  Vrec(i)= %.15lf  Vrec(i+1)= %.15lf\n",Vmin+i*hv, Vmin+(i+1)*hv );
 //  printf("  Erec(j)= %.15lf  Erec(j+1)= %.15lf\n",Emin+j*he, Emin+(j+1)*he );


  int pos_fic=i*(Ne-1) + j;
  //printf("  --i_dis= %d, j_dis= %d \n",i,j );


  Vmx=tabVx[pos_fic];  Emx=tabEx[pos_fic];
  Vmy=tabVy[pos_fic];  Emy=tabEy[pos_fic];
  //printf("Vmx= %lf , Emx= %lf \n",Vmx,Emx );
  //printf("Vmy= %lf , Emy= %lf \n",Vmy,Emy );


  *pV0x=Vmx; *pE0x=Emx;
  *pV0y=Vmy; *pE0y=Emy;

  // Valeur d'initialisation
  double VaT=129.905671875844291e-6;  double EaT=72.168190127265447e3;
  double VbT=127.666112185063994e-6;  double EbT=99.480256119236515e3;
  double VcT=132.538478805442338e-6;  double EcT=136.154346228414312e3;

  
  int nb_zone=tab_Nb_ZONE[pos_fic];
  
  /*
  printf("  nb zone= %d\n",nb_zone );
  for(int i=0; i<nb_zone; i++){
    printf("  tabZONE[%d][pos_fic]= %d",i,tabZONE[i][pos_fic]);
  }
  printf("\n");
  */
  

  int zone;
  int phaseA, phaseB;

  
  // On regarde si lorsque qu'il ya 3 zones ou plus, on a -1 dans les zones
  // Dans ce cas on mets comme graines segV[N] 
  // sinon c'est qu'on est au abord du point triple et dans ce cas on met comme graines les sommets du point triple 
  int isout=0;
  if(nb_zone>=3){  
    for(int i=0; i<nb_zone; i++){
      if(tabZONE[i][pos_fic]==-1){isout=1;}
    }
  }
  // on se sert du isout dans la boucle ou nb_zone>=3



  //printf(" nb_zone= %d\n",nb_zone );
  
  //  ##  1 seule zone ##
  if(nb_zone==1){
    zone=tabZONE[0][pos_fic];
    //printf("   - zone= %d\n",zone );

    if(zone==-1){ 
     if(affichage==1){*pP=-10e9; *pT=0.;  *px=-1;   *pzone=-1;  return 0;}
     else if(affichage==0){printf("Error out of domain : V= %.15lf , E= %.15lf\n",V,E); return 1;}
     else{printf("Erreur dans la valeur de affichage== %d",affichage); return 1;} 
    }
    else if(zone==3){  *pP=fP(V,E,3);  *pT=fT(V,E,3);  *px=-1;  *pzone=3;  return 0;  }
    else if(zone==2){  *pP=fP(V,E,2);  *pT=fT(V,E,2);  *px=-1;  *pzone=2;  return 0;  }
    else if(zone==1){  *pP=fP(V,E,1);  *pT=fT(V,E,1);  *px=-1;  *pzone=1;  return 0;  }
    else if(zone==23){
      err=VE_mixte(Nmax, epsilon, 2, 3, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px);  if(err){return 1;}
      *pP=fP(Vx,Ex,2);  *pT=fT(Vx,Ex,2);  *pzone=23;
      return 0;
    }
    else if(zone==13){
      err=VE_mixte(Nmax, epsilon, 1, 3, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px);  if(err){return 1;}
      *pP=fP(Vx,Ex,1);  *pT=fT(Vx,Ex,1);    *pzone=13;
      return 0;
    }
    else if(zone==12){
      err=VE_mixte(Nmax, epsilon, 1, 2, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px);  if(err){return 1;}
      *pP=fP(Vx,Ex,1);  *pT=fT(Vx,Ex,1);  *pzone=12;
      return 0;
    }
    else if(zone==0){  *pP=PT;  *pT=TT;  *px=-1;  *pzone=0;  return 0; }
  }
  //  ##  2 zones   ##   la première zone est forcément la zone pur (ou triple)
  else if(nb_zone==2){
    zone=tabZONE[0][pos_fic];
    //printf("   - zone= %d\n",zone );
    if(zone==3){
      if(intPOLY(V,E,NC,SC)){  *pP=fP(V,E,3);  *pT=fT(V,E,3);  *px=-1;  *pzone=3;  return 0;  }
      else if(tabZONE[1][pos_fic]==-1){  
        if(affichage==1){*pP=-10e9; *pT=0.;  *px=-1;   *pzone=-1;  return 0;}
        else if(affichage==0){printf("Error out of domain : V= %.15lf , E= %.15lf\n",V,E); return 1;}
      }
      else{
        zone=tabZONE[1][pos_fic];
        phaseB=zone%10;  phaseA=(zone-phaseB)/10;

        err=VE_mixte(Nmax, epsilon, phaseA, phaseB, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px);
        if(err){return 1;}
        *pP=fP(Vx,Ex,phaseA);  *pT=fT(Vx,Ex,phaseA);  *pzone=zone;

        return 0;
      } 
    }
    else if(zone==2){
      if(intPOLY(V,E,NB,SB)){  *pP=fP(V,E,2);  *pT=fT(V,E,2);  *px=-1;   *pzone=2;  return 0;  }
      else if(tabZONE[1][pos_fic]==-1){  
        if(affichage==1){*pP=-10e9; *pT=0.;  *px=-1;   *pzone=-1;  return 0;}
        else if(affichage==0){printf("Error out of domain : V= %.15lf , E= %.15lf\n",V,E); return 1;}
      }
      else{
        zone=tabZONE[1][pos_fic];
        phaseB=zone%10;  phaseA=(zone-phaseB)/10;

        err=VE_mixte(Nmax, epsilon, phaseA, phaseB, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px);
        if(err){return 1;}
        *pP=fP(Vx,Ex,phaseA);  *pT=fT(Vx,Ex,phaseA); *pzone=zone;
        
        return 0;
      } 
    }
    else if(zone==1){
      if(intPOLY(V,E,NA,SA)){  *pP=fP(V,E,1);  *pT=fT(V,E,1);  *px=-1;  *pzone=1;  return 0;  }
      else if(tabZONE[1][pos_fic]==-1){  
        if(affichage==1){*pP=-10e9; *pT=0.;  *px=-1;   *pzone=-1;  return 0;}
        else if(affichage==0){printf("Error out of domain : V= %.15lf , E= %.15lf\n",V,E); return 1;}
      }
      else{
        zone=tabZONE[1][pos_fic];
        phaseB=zone%10;  phaseA=(zone-phaseB)/10;
 
        err=VE_mixte(Nmax, epsilon, phaseA, phaseB, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px);
        if(err){return 1;}
        *pP=fP(Vx,Ex,phaseA);  *pT=fT(Vx,Ex,phaseA);  *pzone=zone;

        return 0;
      } 
    }
    else if(zone==0){ 
      if(ftri(V, E)){  *pP=PT;  *pT=TT;  *px=-1;  *pzone=0;  return 0;  }
      else{
        zone=tabZONE[1][pos_fic];
        phaseB=zone%10;  phaseA=(zone-phaseB)/10;
 
        err=VE_mixte(Nmax, epsilon, phaseA, phaseB, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px);
        if(err){return 1;}
        *pP=fP(Vx,Ex,phaseA);  *pT=fT(Vx,Ex,phaseA);  *pzone=zone;

        return 0;
      } 
    }
    else if(zone=-1){  
      if(affichage==1){*pP=-10e9; *pT=0.; *px=-1; *pzone=-1; return 0;}
      else if(affichage==0){printf("Error out of domain : V= %.15lf , E= %.15lf\n",V,E); return 1;}
      else{printf("Erreur dans la valeur de affichage== %d",affichage); return 1;} 
    }
    else{printf(" Erreur valeur de zone= %d\n",zone );}
  }
  //  ##  3 zones ou plus  ##
  else{
    for(int ind=0; ind<nb_zone-1; ind++){
      zone=tabZONE[ind][pos_fic];
      //printf("   - zone= %d\n",zone );
      if(zone==3){
        if(intPOLY(V,E,NC,SC)){  *pP=fP(V,E,3);  *pT=fT(V,E,3);  *px=-1;  *pzone=3;  return 0;  }
      }
      else if(zone==2){
        if(intPOLY(V,E,NB,SB)){  *pP=fP(V,E,2);  *pT=fT(V,E,2);  *px=-1;  *pzone=2;  return 0;  }
      }
      else if(zone==1){
        if(intPOLY(V,E,NA,SA)){  *pP=fP(V,E,1);  *pT=fT(V,E,1);  *px=-1;  *pzone=1;  return 0;  }
      }
      else if(zone==23){
        if(intPOLY(V,E,NBC,SBC)){
          if(isout){Vmx=SBC[N23][0]; Emx=SBC[N23][1]; Vmy=SBC[N23-1][0]; Emy=SBC[23-1][1];}
          else{Vmx=VbT; Emx=EbT; Vmy=VcT; Emy=EcT;}
          err=VE_mixte(Nmax, epsilon, 2, 3, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px);
          if(err){return 1;}
          *pP=fP(Vx,Ex,2);  *pT=fT(Vx,Ex,2);  *pzone=23;

          return 0;
        }
      }
      else if(zone==13){
        if(intPOLY(V,E,NAC,SAC)){
          if(isout){Vmx=SAC[N13][0]; Emx=SAC[N13][1]; Vmy=SAC[N13-1][0]; Emy=SAC[N13-1][1];}
          else{Vmx=VaT; Emx=EaT; Vmy=VcT; Emy=EcT;}
          err=VE_mixte(Nmax, epsilon, 1, 3, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px);
          if(err){return 1;}
          *pP=fP(Vx,Ex,1);  *pT=fT(Vx,Ex,1); *pzone=13;

          return 0;
        }
      }  
      else if(zone==12){
        if(intPOLY(V,E,NAB,SAB)){
          if(isout){Vmx=SAB[N12][0]; Emx=SAB[N12][1]; Vmy=SAB[N12-1][0]; Emy=SAB[N12-1][1];}
          else{Vmx=VaT; Emx=EaT; Vmy=VbT; Emy=EbT;}
          err=VE_mixte(Nmax, epsilon, 1, 2, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px);
          if(err){return 1;}
          *pP=fP(Vx,Ex,1);  *pT=fT(Vx,Ex,1);  *pzone=12;

          return 0;
        }
      }
      else if(zone==0){   if(ftri(V, E)){  *pP=PT;  *pT=TT;  *px=-1;  *pzone=0;  return 0;  }   }  
      else if(zone=-1){  
        if(affichage==1){*pP=-10e9;  *pT=0.;  *px=-1;  *pzone=-1;  return 0;}
        else if(affichage==0){printf("Error out of domain : V= %.15lf , E= %.15lf\n",V,E);  return 1;}
        else{printf("Erreur dans la valeur de affichage= %d",affichage);  return 1;}
      }
    }
    
    zone=tabZONE[nb_zone-1][pos_fic];
    if(zone==3){  *pP=fP(V,E,3);  *pT=fT(V,E,3);  *px=-1;  *pzone=3;  return 0;  }
    else if(zone==2){  *pP=fP(V,E,2);  *pT=fT(V,E,2);  *px=-1;  *pzone=2;  return 0;  }
    else if(zone==1){   *pP=fP(V,E,1);  *pT=fT(V,E,1);  *px=-1;  *pzone=1;  return 0;  }
    else if(zone==23){
      if(isout){Vmx=SBC[N23][0]; Emx=SBC[N23][1]; Vmy=SBC[N23-1][0]; Emy=SBC[23-1][1];}
      else{Vmx=VbT; Emx=EbT; Vmy=VcT; Emy=EcT;}
      err=VE_mixte(Nmax, epsilon, 2, 3, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px);
      if(err){return 1;}
      *pP=fP(Vx,Ex,2);  *pT=fT(Vx,Ex,2);  *pzone=23;

      return 0;
    }
    else if(zone==13){
      if(isout){Vmx=SAC[N13][0]; Emx=SAC[N13][1]; Vmy=SAC[N13-1][0]; Emy=SAC[N13-1][1];}
      else{Vmx=VaT; Emx=EaT; Vmy=VcT; Emy=EcT;}
      err=VE_mixte(Nmax, epsilon, 1, 3, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px);
      if(err){return 1;}
      *pP=fP(Vx,Ex,1);  *pT=fT(Vx,Ex,1);  *pzone=13;

      return 0;
    }  
    else if(zone==12){
      if(isout){Vmx=SAB[N12][0]; Emx=SAB[N12][1]; Vmy=SAB[N12-1][0]; Emy=SAB[N12-1][1];}
      else{Vmx=VaT; Emx=EaT; Vmy=VbT; Emy=EbT;}
      err=VE_mixte(Nmax, epsilon, 1, 2, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px);
      if(err){return 1;}
      *pP=fP(Vx,Ex,1);  *pT=fT(Vx,Ex,1);  *pzone=12;

      return 0;
    }
    else if(zone==0){  *pP=PT;  *pT=TT;  *px=-1;  *pzone=0;  return 0;  }
    else if(zone=-1){  
      if(affichage==1){*pP=-10e9; *pT=0.;  *px=-1;   *pzone=-1;  return 0;}
      else if(affichage==0){printf("Error out of domain : V= %.15lf , E= %.15lf\n",V,E); return 1;}
      else{printf("Erreur dans la valeur de affichage== %d",affichage); return 1;} 
    }
    else{printf("Erreur dans tabZONE[%d][%d]= %d\n",nb_zone-1,pos_fic,tabZONE[nb_zone-1][pos_fic]); return 1;}
    
    
  }
  // Si on arrive la c'est que la bonne zone n'a pas été trouvé
  printf("Erreur choix zone (lec_benchmark_PT) ");
  printf("  nb_zone= %d",nb_zone );
  printf("  tabZONE[0][%d]= %d",pos_fic,tabZONE[0][pos_fic]);
  for(int ind=1; ind<nb_zone; ind++){
    printf("  tabZONE[%d][%d]= %d",ind,pos_fic,tabZONE[ind][pos_fic]); 
  }
  printf("\n V= %.10lf , E=%.10lf\n",V,E);
  return 1; 
  
}




// Script pour tester le Benchmark
int use_benchmark_PT(int Nv_test, int Ne_test, double Vmin_test, double Emin_test, double Vmax_test, double Emax_test){
  int err;
  double E,V;
  double P,T;
  double x;
  
  // Pour le calcul de la celerité du son
  double V0x, E0x, V0y, E0y;

  // Discrétisation du benchmark/tabulation
  int nb_points;
  int Nv_dis,Ne_dis;
  double Vmin_dis, Vmax_dis;
  double Emin_dis, Emax_dis;

  // Ecriture dans un fichier texte des paramètres de la discretisation
  /* pour pouvoir utiliser la discretisations juste en lisant les fichiers
  Organisation du fichier :
    nb_points
    Nv_dis , Ne_dis
    Vmin_dis , Vmax_dis
    Emin_dis , Emax_dis
  */
  FILE *fparams;
  if((fparams = fopen("fichiers/PARAMS_dis.txt", "r")) == NULL){printf("erreur ouverture fichier PARAMS_dis.txt \n"); return 1;}
  fscanf(fparams, "%d \n",&nb_points );
  fscanf(fparams, "%d %d\n",&Nv_dis,&Ne_dis );
  fscanf(fparams, "%lf %lf\n",&Vmin_dis,&Vmax_dis );
  fscanf(fparams, "%lf %lf\n",&Emin_dis,&Emax_dis );
  fclose(fparams);
  

  // Calcul de la célerité du son
  double Pp,Pm, Tp,Tm;
  double dE,c;
  double dV=1e-8;


  double hv_dis=(Vmax_dis-Vmin_dis)/(Nv_dis-1);
  double he_dis=(Emax_dis-Emin_dis)/(Ne_dis-1);

  double hv_test=(Vmax_test-Vmin_test)/(Nv_test-1);
  double he_test=(Emax_test-Emin_test)/(Ne_test-1);

  printf(" hv_test/hv_dis = %lf percent \n",100*hv_test/hv_dis );
  printf(" he_test/he_dis = %lf percent \n",100*he_test/he_dis );

  // Nombre de points sur les frontières
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3; 
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;

  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );


  double** SA=alloctabd(NA+1,2);
  double** SB=alloctabd(NB+1,2);
  double** SC=alloctabd(NC+1,2); // on ajoute S[N]=S[0]  --> utile pour l'algo
  double** SAB=alloctabd(NAB+1,2);
  double** SAC=alloctabd(NAC+1,2);
  double** SBC=alloctabd(NBC+1,2);

  // fichier binaire    couple (V,E) pour les zones mixte
  FILE *fileVmx, *fileEmx, *fileVmy, *fileEmy;
  if((fileVmx = fopen("fichiers/Vbench_mixte_x.bin", "rb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmx = fopen("fichiers/Ebench_mixte_x.bin", "rb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  if((fileVmy = fopen("fichiers/Vbench_mixte_y.bin", "rb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmy = fopen("fichiers/Ebench_mixte_y.bin", "rb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  
  double* tabVx=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double)); double* tabEx=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double));
  double* tabVy=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double)); double* tabEy=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double));

  // lecture des tableaux
  fread(tabVx,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileVmx);
  fread(tabEx,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileEmx);
  fread(tabVy,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileVmy);
  fread(tabEy,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileEmy);

  fclose(fileVmx);  fclose(fileEmx);
  fclose(fileVmy);  fclose(fileEmy);
  
  // Ouverture du fichier ZONE
  FILE *fileZONE_0,*fileZONE_1,*fileZONE_2,*fileZONE_3,*fileZONE_4;
  if((fileZONE_0 = fopen("fichiers/ZONEbench_0.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 0\n"); return 1;}
  if((fileZONE_1 = fopen("fichiers/ZONEbench_1.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 1\n"); return 1;}
  if((fileZONE_2 = fopen("fichiers/ZONEbench_2.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 2\n"); return 1;}
  if((fileZONE_3 = fopen("fichiers/ZONEbench_3.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 3\n"); return 1;}
  if((fileZONE_4 = fopen("fichiers/ZONEbench_4.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 4\n"); return 1;}
  
  int** tabZONE=alloctab(5,(Nv_dis-1)*(Ne_dis-1));
  
  // Ecriture des tableaux
  fread(tabZONE[0],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_0);
  fread(tabZONE[1],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_1);
  fread(tabZONE[2],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_2);
  fread(tabZONE[3],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_3);
  fread(tabZONE[4],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_4);

  fclose(fileZONE_0);  fclose(fileZONE_1);  fclose(fileZONE_2);
  fclose(fileZONE_3);  fclose(fileZONE_4);


  // fichier binaire NOMBRE DE ZONE
  FILE *file_Nb_ZONE;
  if((file_Nb_ZONE = fopen("fichiers/Nb_ZONEbench.bin", "rb")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  
  int* tab_Nb_ZONE=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(int)); 

  // Ecriture des tableaux
  fread(tab_Nb_ZONE,sizeof(int),(Nv_dis-1)*(Ne_dis-1),file_Nb_ZONE);

  fclose(file_Nb_ZONE);


  // fichier binaire POLYGONES Frontières de zones
  FILE *fileA_V,*fileB_V,*fileC_V,*fileAB_V,*fileAC_V,*fileBC_V;
  FILE *fileA_E,*fileB_E,*fileC_E,*fileAB_E,*fileAC_E,*fileBC_E;
  if((fileA_V = fopen("fichiers/A_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE A\n"); return 1;}     if((fileA_E = fopen("fichiers/A_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE A\n"); return 1;}
  if((fileB_V = fopen("fichiers/B_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE B\n"); return 1;}     if((fileB_E = fopen("fichiers/B_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE B\n"); return 1;}
  if((fileC_V = fopen("fichiers/C_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE C\n"); return 1;}     if((fileC_E = fopen("fichiers/C_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE C\n"); return 1;}
  if((fileAB_V = fopen("fichiers/AB_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AB\n"); return 1;}  if((fileAB_E = fopen("fichiers/AB_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AB\n"); return 1;}
  if((fileAC_V = fopen("fichiers/AC_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AC\n"); return 1;}  if((fileAC_E = fopen("fichiers/AC_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AC\n"); return 1;}
  if((fileBC_V = fopen("fichiers/BC_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE BC\n"); return 1;}  if((fileBC_E = fopen("fichiers/BC_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE BC\n"); return 1;}
  
  double* SA_V=malloc((NA+1)*sizeof(double));    double* SA_E=malloc((NA+1)*sizeof(double));
  double* SB_V=malloc((NB+1)*sizeof(double));    double* SB_E=malloc((NB+1)*sizeof(double));
  double* SC_V=malloc((NC+1)*sizeof(double));    double* SC_E=malloc((NC+1)*sizeof(double));
  double* SAB_V=malloc((NAB+1)*sizeof(double));  double* SAB_E=malloc((NAB+1)*sizeof(double));
  double* SAC_V=malloc((NAC+1)*sizeof(double));  double* SAC_E=malloc((NAC+1)*sizeof(double));
  double* SBC_V=malloc((NBC+1)*sizeof(double));  double* SBC_E=malloc((NBC+1)*sizeof(double));

  // Ecriture des tableaux
  fread(SA_V,sizeof(double),NA+1,fileA_V);     fread(SA_E,sizeof(double),NA+1,fileA_E);
  fread(SB_V,sizeof(double),NB+1,fileB_V);     fread(SB_E,sizeof(double),NB+1,fileB_E);
  fread(SC_V,sizeof(double),NC+1,fileC_V);     fread(SC_E,sizeof(double),NC+1,fileC_E);
  fread(SAB_V,sizeof(double),NAB+1,fileAB_V);  fread(SAB_E,sizeof(double),NAB+1,fileAB_E);
  fread(SAC_V,sizeof(double),NAC+1,fileAC_V);  fread(SAC_E,sizeof(double),NAC+1,fileAC_E);
  fread(SBC_V,sizeof(double),NBC+1,fileBC_V);  fread(SBC_E,sizeof(double),NBC+1,fileBC_E);

  for(int i=0; i<NA+1; i++){  SA[i][0]=SA_V[i];  SA[i][1]=SA_E[i];  }
  for(int i=0; i<NB+1; i++){  SB[i][0]=SB_V[i];  SB[i][1]=SB_E[i];  }
  for(int i=0; i<NC+1; i++){  SC[i][0]=SC_V[i];  SC[i][1]=SC_E[i];  }

  for(int i=0; i<NAB+1; i++){  SAB[i][0]=SAB_V[i];  SAB[i][1]=SAB_E[i];  }
  for(int i=0; i<NAC+1; i++){  SAC[i][0]=SAC_V[i];  SAC[i][1]=SAC_E[i];  }
  for(int i=0; i<NBC+1; i++){  SBC[i][0]=SBC_V[i];  SBC[i][1]=SBC_E[i];  }

  fclose(fileA_V);   fclose(fileB_V);   fclose(fileC_V);
  fclose(fileA_E);   fclose(fileB_E);   fclose(fileC_E);
  fclose(fileAB_V);  fclose(fileAC_V);  fclose(fileBC_V);
  fclose(fileAB_E);  fclose(fileAC_E);  fclose(fileBC_E);

  free(SA_V); free(SB_V); free(SC_V); free(SAB_V); free(SAC_V); free(SBC_V);
  free(SA_E); free(SB_E); free(SC_E); free(SAB_E); free(SAC_E); free(SBC_E);








 //////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////
  ///////////
  // Ecriture des sommets des poly pour vérif
  FILE *fSA, *fSB, *fSC, *fSAB, *fSAC, *fSBC;
  if((fSA = fopen("fichiers/pSA.txt", "w+")) == NULL){printf("erreur ouverture fichier pSA\n");return 1;}
  if((fSB = fopen("fichiers/pSB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSB\n");return 1;}
  if((fSC = fopen("fichiers/pSC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSC\n");return 1;}
  if((fSAB = fopen("fichiers/pSAB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAB\n");return 1;}
  if((fSAC = fopen("fichiers/pSAC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAC\n");return 1;}
  if((fSBC = fopen("fichiers/pSBC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSBC\n");return 1;}
  
  for(int i=0; i<NA+1; i++){   fprintf(fSA,"%.10lf %.10lf\n", SA[i][0], SA[i][1] );  }
  for(int i=0; i<NB+1; i++){   fprintf(fSB,"%.10lf %.10lf\n", SB[i][0], SB[i][1] );  }
  for(int i=0; i<NC+1; i++){   fprintf(fSC,"%.10lf %.10lf\n", SC[i][0], SC[i][1] );  }
  for(int i=0; i<NAB+1; i++){  fprintf(fSAB,"%.10lf %.10lf\n",SAB[i][0],SAB[i][1] );  }
  for(int i=0; i<NAC+1; i++){  fprintf(fSAC,"%.10lf %.10lf\n",SAC[i][0],SAC[i][1] );  }
  for(int i=0; i<NBC+1; i++){  fprintf(fSBC,"%.10lf %.10lf\n",SBC[i][0],SBC[i][1] );  }

  fclose(fSA);   fclose(fSB);   fclose(fSC);  fclose(fSAB);  fclose(fSAC);  fclose(fSBC);
 //////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////











  //  ****************************************
  // Pour tester les valeurs
  FILE *fPtest, *fTtest, *fXtest, *fCtest, *fGtest;
  if((fPtest = fopen("fichiers/Ptest.txt", "w+")) == NULL){printf("erreur ouverture fichier P test\n");return 1;}
  if((fTtest = fopen("fichiers/Ttest.txt", "w+")) == NULL){printf("erreur ouverture fichier T test\n");return 1;}
  if((fXtest = fopen("fichiers/Xtest.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  if((fCtest = fopen("fichiers/Ctest.txt", "w+")) == NULL){printf("erreur ouverture fichier C test\n");return 1;}
  if((fGtest = fopen("fichiers/Gtest.txt", "w+")) == NULL){printf("erreur ouverture fichier G test\n");return 1;}
  FILE *fVtest, *fEtest;
  if((fVtest = fopen("fichiers/Vtest.txt", "w+")) == NULL){printf("erreur ouverture fichier V test \n");return 1;}
  if((fEtest = fopen("fichiers/Etest.txt", "w+")) == NULL){printf("erreur ouverture fichier E test \n");return 1;}
  for(int j=0; j<Ne_test; j++){
      E=Emin_test+j*he_test;
      fprintf(fEtest, "%.15lf ",E);
  }
  for(int j=0; j<Nv_test; j++){
      V=Vmin_test+j*hv_test;
      fprintf(fVtest, "%.15lf ",V);
  }
  fclose(fVtest);  fclose(fEtest);


  int zone;
  FILE *fZONEtest;
  if((fZONEtest = fopen("fichiers/ZONEtest.txt", "w+")) == NULL){printf("erreur ouverture fichier ZONE test\n");return 1;}
  //  ****************************************
  
  //  $$$$$$$   Ecriture texte pour verif   $$$$$$$$$$
  FILE *file_Nb_ZONE_test;
  if((file_Nb_ZONE_test = fopen("fichiers/Nb_ZONEbench.txt", "w+")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  for(int i=0; i<(Nv_dis-1)*(Ne_dis-1); i++){
    fprintf(file_Nb_ZONE_test, "%d \n",tab_Nb_ZONE[i] );
  }
  fclose(file_Nb_ZONE_test);

  FILE *file_V_dis;
  if((file_V_dis = fopen("fichiers/Vbench_dis.txt", "w+")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  for(int i=0; i<(Nv_dis-1); i++){
    fprintf(file_V_dis, "%.15lf \n",Vmin_dis+i*hv_dis+hv_dis/2 );
  }
  fclose(file_V_dis);
  
  FILE *file_E_dis;
  if((file_E_dis = fopen("fichiers/Ebench_dis.txt", "w+")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  for(int i=0; i<(Ne_dis-1); i++){
    fprintf(file_E_dis, "%.15lf \n",Emin_dis+i*he_dis+he_dis/2 );
  }
  fclose(file_E_dis);
  //  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  double G;

  int Nmax=1e2;
  double epsilon=1e-8;
  // ***** Boucle *****
  for(int i=0; i<Nv_test; i++){
    printf("  - i= %d\n",i );
    V=Vmin_test+i*hv_test;
    for(int j=0; j<Ne_test; j++){
      //printf(" - i= %d, j= %d \n",i,j );
      E=Emin_test+j*he_test;
      err=lec_benchmark_PT(1,V, E, nb_points, SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE, 
                           Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, &P, &T, &x, &zone, &V0x, &E0x, &V0y, &E0y);
      if(err){return err;}
      
      if(zone==-1){c=-1e3; G=-1e1;}
      else{
        err=c_celerite_son(Nmax, epsilon,zone, P, V, E, dV, V0x, E0x, V0y, E0y, &c, &G);
        if(err){return err;}
      }


      fprintf(fPtest,"%.10lf\n",P);
      fprintf(fTtest,"%.10lf\n",T);
      fprintf(fXtest,"%.10lf\n",x);
      fprintf(fZONEtest,"%d\n",zone);
      fprintf(fCtest,"%.10lf\n",c);      
      fprintf(fGtest,"%.10lf\n",G); // dérivée fondamentale      
    }
  }

  fclose(fPtest);  fclose(fTtest);  fclose(fXtest); fclose(fCtest); fclose(fGtest);
  

  // Liberation de la mémoire
  // Liberation
  free(tabVx);   free(tabEx);
  free(tabVy);   free(tabEy);
  
  freetab(tabZONE);  free(tab_Nb_ZONE);

  freetab(SA);  freetab(SB);  freetab(SC);
  freetab(SAB);  freetab(SAC);  freetab(SBC);
  
  return 0;
}



//  Remplie les tableaux suivants pour les utiliser dans le main 
/*

*/
int lec_tabulation(int nb_points, int Nv_dis, int Ne_dis, int** tabZONE, int* tab_Nb_ZONE, double* tabVx, double* tabEx, double* tabVy, double* tabEy,
                   double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC){
  


  // Nombre de points sur les frontières
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3; 
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;

  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );
  
  //printf(" ok 1\n");
  // fichier binaire    couple (V,E) pour les zones mixte
  FILE *fileVmx, *fileEmx, *fileVmy, *fileEmy;
  if((fileVmx = fopen("fichiers/Vbench_mixte_x.bin", "rb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmx = fopen("fichiers/Ebench_mixte_x.bin", "rb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  if((fileVmy = fopen("fichiers/Vbench_mixte_y.bin", "rb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmy = fopen("fichiers/Ebench_mixte_y.bin", "rb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  
  // lecture des tableaux
  fread(tabVx,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileVmx);
  fread(tabEx,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileEmx);
  fread(tabVy,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileVmy);
  fread(tabEy,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileEmy);

  fclose(fileVmx);  fclose(fileEmx);
  fclose(fileVmy);  fclose(fileEmy);
  
  //printf(" ok 2\n"); 
  // Ouverture du fichier ZONE
  FILE *fileZONE_0,*fileZONE_1,*fileZONE_2,*fileZONE_3,*fileZONE_4;
  if((fileZONE_0 = fopen("fichiers/ZONEbench_0.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 0\n"); return 1;}
  if((fileZONE_1 = fopen("fichiers/ZONEbench_1.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 1\n"); return 1;}
  if((fileZONE_2 = fopen("fichiers/ZONEbench_2.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 2\n"); return 1;}
  if((fileZONE_3 = fopen("fichiers/ZONEbench_3.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 3\n"); return 1;}
  if((fileZONE_4 = fopen("fichiers/ZONEbench_4.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 4\n"); return 1;}
  
  // Ecriture des tableaux
  fread(tabZONE[0],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_0);
  fread(tabZONE[1],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_1);
  fread(tabZONE[2],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_2);
  fread(tabZONE[3],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_3);
  fread(tabZONE[4],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_4);

  fclose(fileZONE_0);  fclose(fileZONE_1);  fclose(fileZONE_2);
  fclose(fileZONE_3);  fclose(fileZONE_4);
  
  //printf(" ok 3\n");
  // fichier binaire NOMBRE DE ZONE
  FILE *file_Nb_ZONE;
  if((file_Nb_ZONE = fopen("fichiers/Nb_ZONEbench.bin", "rb")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  
  // Ecriture des tableaux
  fread(tab_Nb_ZONE,sizeof(int),(Nv_dis-1)*(Ne_dis-1),file_Nb_ZONE);

  fclose(file_Nb_ZONE);
  
  //printf(" ok 4\n");
  // fichier binaire POLYGONES Frontières de zones
  FILE *fileA_V,*fileB_V,*fileC_V,*fileAB_V,*fileAC_V,*fileBC_V;
  FILE *fileA_E,*fileB_E,*fileC_E,*fileAB_E,*fileAC_E,*fileBC_E;
  if((fileA_V = fopen("fichiers/A_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE A\n"); return 1;}     if((fileA_E = fopen("fichiers/A_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE A\n"); return 1;}
  if((fileB_V = fopen("fichiers/B_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE B\n"); return 1;}     if((fileB_E = fopen("fichiers/B_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE B\n"); return 1;}
  if((fileC_V = fopen("fichiers/C_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE C\n"); return 1;}     if((fileC_E = fopen("fichiers/C_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE C\n"); return 1;}
  if((fileAB_V = fopen("fichiers/AB_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AB\n"); return 1;}  if((fileAB_E = fopen("fichiers/AB_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AB\n"); return 1;}
  if((fileAC_V = fopen("fichiers/AC_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AC\n"); return 1;}  if((fileAC_E = fopen("fichiers/AC_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AC\n"); return 1;}
  if((fileBC_V = fopen("fichiers/BC_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE BC\n"); return 1;}  if((fileBC_E = fopen("fichiers/BC_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE BC\n"); return 1;}
  

  double* SA_V=malloc((NA+1)*sizeof(double));    double* SA_E=malloc((NA+1)*sizeof(double));
  double* SB_V=malloc((NB+1)*sizeof(double));    double* SB_E=malloc((NB+1)*sizeof(double));
  double* SC_V=malloc((NC+1)*sizeof(double));    double* SC_E=malloc((NC+1)*sizeof(double));
  double* SAB_V=malloc((NAB+1)*sizeof(double));  double* SAB_E=malloc((NAB+1)*sizeof(double));
  double* SAC_V=malloc((NAC+1)*sizeof(double));  double* SAC_E=malloc((NAC+1)*sizeof(double));
  double* SBC_V=malloc((NBC+1)*sizeof(double));  double* SBC_E=malloc((NBC+1)*sizeof(double));

  // Ecriture des tableaux
  fread(SA_V,sizeof(double),NA+1,fileA_V);     fread(SA_E,sizeof(double),NA+1,fileA_E);
  fread(SB_V,sizeof(double),NB+1,fileB_V);     fread(SB_E,sizeof(double),NB+1,fileB_E);
  fread(SC_V,sizeof(double),NC+1,fileC_V);     fread(SC_E,sizeof(double),NC+1,fileC_E);
  fread(SAB_V,sizeof(double),NAB+1,fileAB_V);  fread(SAB_E,sizeof(double),NAB+1,fileAB_E);
  fread(SAC_V,sizeof(double),NAC+1,fileAC_V);  fread(SAC_E,sizeof(double),NAC+1,fileAC_E);
  fread(SBC_V,sizeof(double),NBC+1,fileBC_V);  fread(SBC_E,sizeof(double),NBC+1,fileBC_E);

  for(int i=0; i<NA+1; i++){  SA[i][0]=SA_V[i];  SA[i][1]=SA_E[i];  }
  for(int i=0; i<NB+1; i++){  SB[i][0]=SB_V[i];  SB[i][1]=SB_E[i];  }
  for(int i=0; i<NC+1; i++){  SC[i][0]=SC_V[i];  SC[i][1]=SC_E[i];  }

  for(int i=0; i<NAB+1; i++){  SAB[i][0]=SAB_V[i];  SAB[i][1]=SAB_E[i];  }
  for(int i=0; i<NAC+1; i++){  SAC[i][0]=SAC_V[i];  SAC[i][1]=SAC_E[i];  }
  for(int i=0; i<NBC+1; i++){  SBC[i][0]=SBC_V[i];  SBC[i][1]=SBC_E[i];  }

  fclose(fileA_V);   fclose(fileB_V);   fclose(fileC_V);
  fclose(fileA_E);   fclose(fileB_E);   fclose(fileC_E);
  fclose(fileAB_V);  fclose(fileAC_V);  fclose(fileBC_V);
  fclose(fileAB_E);  fclose(fileAC_E);  fclose(fileBC_E);

  free(SA_V); free(SB_V); free(SC_V); free(SAB_V); free(SAC_V); free(SBC_V);
  free(SA_E); free(SB_E); free(SC_E); free(SAB_E); free(SAC_E); free(SBC_E);

  
  //printf(" ok fin lec\n");
  return 0;
}

// #######################  FONCTIONS POUR LE CALCUL DE P,T et C avec tabulation + calcul #############################



// #######################  FONCTIONS POUR LE CALCUL DE P,T et C avec calcul uniquement (sans tabulation) #############################

int graine_zone_mixte_modif(double V, double E, int N, double** Sp, double* pVmx, double* pEmx, double* pVmy, double* pEmy){
  int res;

  int Ntot=2*N-1;
  
  // reparage du quadrilatere
  double** S=alloctabd(5,2);
  int ind=0;
  int testNO=1;
  for(int i=0; i<=N-2; i++){
    S[0][0]=Sp[Ntot-i][0];     S[0][1]=Sp[Ntot-i][1];
    S[1][0]=Sp[Ntot-(i+1)][0]; S[1][1]=Sp[Ntot-(i+1)][1];
    S[2][0]=Sp[i+1][0];        S[2][1]=Sp[i+1][1];
    S[3][0]=Sp[i][0];          S[3][1]=Sp[i][1];
    S[4][0]=Sp[Ntot-i][0];     S[4][1]=Sp[Ntot-i][1];

    /*
    if( fabs(Sp[i][0] - Sp[i+1][0])>fabs(Sp[i][0] - V) || fabs(Sp[N+i][0] - Sp[N+i+1][0])>fabs(Sp[i][0] - V) || fabs(Sp[i][1] - Sp[i+1][1])>fabs(Sp[i][1] - E) || fabs(Sp[N+i][1] - Sp[N+i+1][1])>fabs(Sp[i][1] - E)   ){  
      printf(" V= %g; E= %g; \n", V,E);
      printf(" Sp[N+%d][0]= %g;  Sp[N+%d][1]= %g; \n",i,Sp[N+i][0],i,Sp[N+i][1]);
      printf(" Sp[N+%d+1][0]= %g;  Sp[N+%d+1][1]= %g; \n",i,Sp[N+i+1][0],i,Sp[N+i+1][1]);
      printf(" Sp[%d+1][0]= %g;  Sp[%d+1][1]= %g; \n",i,Sp[i+1][0],i,Sp[i+1][1]);
      printf(" Sp[%d][0]= %g;  Sp[%d][1]= %g; \n",i,Sp[i][0],i,Sp[i][1]);
    }
    */

    if(intPOLY(V, E, 4, S)){ind=i; testNO=0;   break;  } // ou i=N
  
  }

  if(testNO){
    printf("\n");
    printf("le point est dans aucun des quadrilatere (graine_zone_mixte)\n");
    printf(" ind= %d\n",ind );
    printf("V= %g , E= %g \n",V,E); 
    return 1;
  }

  *pVmx=Sp[N+ind][0]; *pEmx=Sp[N+ind][1];
  *pVmy=Sp[ind][0];   *pEmy=Sp[ind][1];

  
  free(S);

}


// #######################  FONCTIONS POUR LE CALCUL DE P,T et C avec tabulation uniquement #############################


/// Petite fonction pour le schema Godunov en energie totale
int epsilon_etain(double V, double P, double T, double* pE, int nb_points,
                  double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC ){
  double Vf,Ef;
  int phase;
  // nb_points est l'argument d'entree de la fonction
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3; 
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;
  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );
  
  double** SA_PT=alloctabd(NA+1,2);
  double** SB_PT=alloctabd(NB+1,2);
  double** SC_PT=alloctabd(NC+1,2);
  if(SA_PT==NULL){printf("Problème d'allocation (epsilon_etain)\n"); return 1;}
  if(SB_PT==NULL){printf("Problème d'allocation (epsilon_etain)\n"); return 1;}
  if(SC_PT==NULL){printf("Problème d'allocation (epsilon_etain)\n"); return 1;}

  for(int i=0; i<NA+1; i++){
    Vf=SA[i][0];Ef=SA[i][1];
    SA_PT[i][0]=fP(Vf,Ef,1);  
    SA_PT[i][1]=fT(Vf,Ef,1);  
  }
  for(int i=0; i<NB+1; i++){
    Vf=SB[i][0];Ef=SB[i][1];
    SB_PT[i][0]=fP(Vf,Ef,2);  
    SB_PT[i][1]=fT(Vf,Ef,2);  
  }
  for(int i=0; i<NC+1; i++){
    Vf=SC[i][0];Ef=SC[i][1];
    SC_PT[i][0]=fP(Vf,Ef,3);  
    SC_PT[i][1]=fT(Vf,Ef,3);  
  }

  /*
  printf(" SA[%d][0]= %.15lf , SA[%d][0]= %.15lf\n",0, SA[0][0]*1e6, NA, SA[NA][0]*1e6 );
  printf(" SA[%d][1]= %.15lf , SA[%d][1]= %.15lf\n",0, SA[0][1], NA, SA[NA][1] );

  printf(" SB[%d][0]= %.15lf , SB[%d][0]= %.15lf\n",0, SB[0][0]*1e6, NB, SB[NB][0]*1e6 );
  printf(" SB[%d][1]= %.15lf , SB[%d][1]= %.15lf\n",0, SB[0][1], NB, SB[NB][1] );

  printf(" SC[%d][0]= %.15lf , SC[%d][0]= %.15lf\n",0, SC[0][0]*1e6, NC, SC[NC][0]*1e6 );
  printf(" SC[%d][1]= %.15lf , SC[%d][1]= %.15lf\n",0, SC[0][1], NC, SC[NC][1] );


  printf(" SA_PT[%d][0]= %.15lf , SA_PT[%d][0]= %.15lf\n",0,SA_PT[0][0],NA,SA_PT[NA][0] );
  printf(" SA_PT[%d][1]= %.15lf , SA_PT[%d][1]= %.15lf\n",0,SA_PT[0][1],NA,SA_PT[NA][1] );

  printf(" SB_PT[%d][0]= %.15lf , SB_PT[%d][0]= %.15lf\n",0, SB_PT[0][0], NB, SB_PT[NB][0] );
  printf(" SB_PT[%d][1]= %.15lf , SB_PT[%d][1]= %.15lf\n",0, SB_PT[0][1], NB, SB_PT[NB][1] );

  printf(" SC_PT[%d][0]= %.15lf , SC_PT[%d][0]= %.15lf\n",0, SC_PT[0][0], NC, SC_PT[NC][0] );
  printf(" SC_PT[%d][1]= %.15lf , SC_PT[%d][1]= %g\n",0, SC_PT[0][1], NC, SC_PT[NC][1] );
  */

  if(intPOLY(P,T,NA,SA_PT)){phase=1;}
  else if(intPOLY(P,T,NB,SB_PT)){phase=2;}
  else if(intPOLY(P,T,NC,SC_PT)){phase=3;}
  else{printf("Error out of domain P= %lf Gpa , T= %lf K\n",P*1e-9,T ); return 1;}

  *pE=fE_VP(V,P,phase);
  return 0;

}



// ********                                                                  ************
///////////*************  CINEMATIQUE DE CHANGEMENT DE PHASE *************////////////////
// ********                                                                  ************


// *******  Newton pour la cinétique de chgmt de phase  ********

/*  Valeur d'entree :
      Vmel, Emel : 
      lambda : fraction massique
      Va,Ea, Vb, Eb : Vraiables


   Valeur sortie
     Vx,Ex, Vy,Ey, 

   double* A et InvA permet une seule allocation 

*/
void Newton_cin_phase(double* A, double* InvA, double* Delta, int phaseX, int phaseY, double Vmel, double Emel, double lambda, double Vx, double Ex, double Vy, double Ey, double* pVX, double* pEX, double* pVY, double* pEY, double* dX){
  int dim=4;

  // Grandeurs thermodynamiques
  double Px, Tx, dPvx, dPex, dTvx, dTex;
  double Py, Ty, dPvy, dPey, dTvy, dTey;
  // phase X
  Px=fP(Vx,Ex, phaseX); Tx=fT(Vx,Ex, phaseX);
  dPvx=dPdV(Vx,phaseX); dPex=dPdE(phaseX);
  dTvx=dTdV(Vx,phaseX); dTex=dTdE(phaseX);
  // phase Y
  Py=fP(Vy,Ey, phaseY); Ty=fT(Vy,Ey, phaseY);
  dPvy=dPdV(Vy,phaseY); dPey=dPdE(phaseY);
  dTvy=dTdV(Vy,phaseY); dTey=dTdE(phaseY);
  
  
  //printf("Px= %g, Tx= %g\n",Px,Tx );
  //printf("Py= %g, Ty= %g\n",Py,Ty );
  
  
  Delta[0]=Py-Px;   Delta[1]=Ty-Tx;  
  Delta[2]=Vmel - (1.-lambda)*Vx - lambda*Vy;
  Delta[3]=Emel - (1.-lambda)*Ex - lambda*Ey;
   
  /*
  for(int j=0; j<dim; j++){
    printf("Delta[%d]= %g  ",j, Delta[j] );
  }
  printf("\n");
  */
  
  
  // on met le second menbre dans la dernière colonne
  //matrice A
  A[dim*0+0]=dPvx;       A[dim*0+1]=dPex;       A[dim*0+2]=-dPvy;   A[dim*0+3]=-dPey; 
  
  A[dim*1+0]=dTvx;       A[dim*1+1]=dTex;       A[dim*1+2]=-dTvy;   A[dim*1+3]=-dTey;
  
  A[dim*2+0]=1.-lambda;  A[dim*2+1]=0.;         A[dim*2+2]=lambda;  A[dim*2+3]=0.;     
  
  A[dim*3+0]=0.;         A[dim*3+1]=1.-lambda;  A[dim*3+2]=0.;      A[dim*3+3]=lambda;    
  
  /*
  printf("  * A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      printf("A[%d][%d]= %g  ",i,j, A[dim*i+j] );
    }
    printf("\n");
  }
  printf("\n");  
  */


  double det_a=inverse_matrice_pivot(A, dim, InvA);
  //printf("detA= %g\n",det_a );
  if ( fabs(det_a)<1e-8){printf("A non inversible : detA= %g\n", det_a);}
  
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
  
  
  // Mutiplcation matricielle de A et InvA  TEST
  /*
  double res;
  //printf("  * InvA*A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      res=0;
      for(int k=0; k<dim; k++){
        res+=InvA[dim*i+k]*A[dim*k+j];
      }
      //printf("A*InvA[%d][%d]= %g ",i,j,res );
    }
    //printf("\n");
  }
  //printf("\n");
  */


  // Mutiplcation matricielle
  
  for(int i=0; i<dim; i++){
    dX[i]=0;
    for(int j=0; j<dim; j++){
      dX[i]+=InvA[dim*i+j]*Delta[j];
    }
    //printf("  dX[%d]= %g\n",i,dX[i] );
  }
  //printf("\n");
  


  double correc=1.;

  *pVX=Vx+correc*dX[0];
  *pEX=Ex+correc*dX[1];
  *pVY=Vy+correc*dX[2];
  *pEY=Ey+correc*dX[3];

  /*
  printf("  Vx= %.10lf, Ex= %.10lf   Vy= %.10lf, Ey= %.10lf    x= %.10lf\n",Vx,Ex,Vy,Ey,x );
  printf("  (Vmel-Vx)/(Vy-Vx)= %.10lf  ,    (Emel-Ex)/(Ey-Ex)= %.10lf\n", (Vmel-Vx)/(Vy-Vx), (Emel-Ex)/(Ey-Ex));
  printf("  (1-x)*Vx + x*Vy - Vmel= %.10lf  ,    (1-x)*Ex + x*Ey - Emel= %.10lf\n", (1-x)*Vx + x*Vy - Vmel, (1-x)*Ex + x*Ey - Emel);
  */

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
int VE_mixte_cin_phase(int Nmax, double epsilon, int phaseX, int phaseY, double lambda, double V0x, double E0x, double V0y, double E0y, double Vmel, double Emel, double* pVx, double* pEx, double* pVy, double* pEy){
  int n=1;
  int dim=4;
   
  // printf(" epsilon= %g\n",epsilon);
  /*
  printf("  phaseX= %d , phaseY= %d\n",phaseX,phaseY );
  printf("  Vmel= %g , Emel= %g\n",Vmel, Emel );
  */


  double* A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));
  double* Delta=malloc(dim*sizeof(double));
  double* dX=malloc(dim*sizeof(double));
  
  double Px, Tx;
  double Py, Ty;

  double Vx, Ex, Vy, Ey;
  double VX, EX, VY, EY;

  Vx=V0x; Ex=E0x; 
  Vy=V0y; Ey=E0y; 

  // phase X
  Px=fP(Vx,Ex, phaseX); Tx=fT(Vx,Ex, phaseX);
  // phase Y
  Py=fP(Vy,Ey, phaseY); Ty=fT(Vy,Ey, phaseY);


  double critere;
  
  critere=1.;
  //printf("  critere= %.14lf\n",critere );
    
  n=0;
  while(critere>=epsilon && n<Nmax){
    
    Newton_cin_phase(A, InvA, Delta, phaseX, phaseY, Vmel, Emel, lambda, Vx, Ex, Vy, Ey, &VX, &EX, &VY, &EY, dX);
    

    Vx=VX; Ex=EX;
    Vy=VY; Ey=EY;
    
    /*
    // phase X
    Px=fP(Vx,Ex, phaseX); Tx=fT(Vx,Ex, phaseX);
    // phase Y
    Py=fP(Vy,Ey, phaseY); Ty=fT(Vy,Ey, phaseY);
    */

    //critere=fabs(Gx-Gy)+fabs(Px-Py)+fabs(Ty-Tx);
    critere=fabs(dX[0]/Vx)+fabs(dX[1]/Ex)+fabs(dX[2]/Vy)+fabs(dX[3]/Ey);
    //critere=fabs(Delta[0])+fabs(Delta[1])+fabs(Delta[2])+fabs(Delta[3]);
    //printf("  n= %d    critere= %g\n",n,critere );

    
    n++;
  }
  //printf(" --n= %d, critere= %g (VE_mixte_cin_phase)\n",n,critere );

  if(n==Nmax){printf("Newton_mixte_cin_phase non convergent : -n= %d , critere= %g (VE_mixte_cin_phase)\n",n,critere ); return 1;}

  *pVx=Vx; *pEx=Ex;
  *pVy=Vy; *pEy=Ey;

  //printf("  ok fin boucle i\n");
  free(A); 
  free(InvA);
  free(Delta);
  free(dX);
  
  return 0;
}



int Newton_cin_3_phase_all_dom(double* A, double* InvA, double* Delta, double Vmel, double Emel, double* lambda, double Vx, double Ex, double Vy, double Ey, double Vz, double Ez, double* pVX, double* pEX, double* pVY, double* pEY, double* pVZ, double* pEZ, double* dX){
  int dim=6;

  double lambdaX,lambdaY,lambdaZ;
  lambdaX=lambda[0];
  lambdaY=lambda[1];
  lambdaZ=lambda[2];

  int phaseX=1, phaseY=2, phaseZ=3;

  // Grandeurs thermodynamiques
  double Px, Tx, dPvx, dPex, dTvx, dTex;
  double Py, Ty, dPvy, dPey, dTvy, dTey;
  double Pz, Tz, dPvz, dPez, dTvz, dTez;
  // phase X
  Px=fP(Vx,Ex, phaseX); Tx=fT(Vx,Ex, phaseX);
  dPvx=dPdV(Vx,phaseX); dPex=dPdE(phaseX);
  dTvx=dTdV(Vx,phaseX); dTex=dTdE(phaseX);
  // phase Y
  Py=fP(Vy,Ey, phaseY); Ty=fT(Vy,Ey, phaseY);
  dPvy=dPdV(Vy,phaseY); dPey=dPdE(phaseY);
  dTvy=dTdV(Vy,phaseY); dTey=dTdE(phaseY);
  // phase Z
  Pz=fP(Vz,Ez, phaseZ); Tz=fT(Vz,Ez, phaseZ);
  dPvz=dPdV(Vz,phaseZ); dPez=dPdE(phaseZ);
  dTvz=dTdV(Vz,phaseZ); dTez=dTdE(phaseZ);
  
  
  //printf("Px= %g, Tx= %g\n",Px,Tx );
  //printf("Py= %g, Ty= %g\n",Py,Ty );
  
  
  Delta[0]=Py-Px;   Delta[1]=Pz-Px;
  Delta[2]=Ty-Tx;   Delta[3]=Tz-Tx;
  Delta[4]=Vmel - lambdaX*Vx - lambdaY*Vy - lambdaZ*Vz;
  Delta[5]=Emel - lambdaX*Ex - lambdaY*Ey - lambdaZ*Ez;
   
  /*
  for(int j=0; j<dim; j++){
    printf("Delta[%d]= %g  ",j, Delta[j] );
  }
  printf("\n");
  */
  
  
  // on met le second menbre dans la dernière colonne
  //matrice A

  A[dim*0+0]=dPvx;     A[dim*0+1]=dPex;     A[dim*0+2]=-dPvy;    A[dim*0+3]=-dPey;    A[dim*0+4]=0.;       A[dim*0+5]=0.;  

  A[dim*1+0]=dPvx;     A[dim*1+1]=dPex;     A[dim*1+2]=0.;       A[dim*1+3]=0.;       A[dim*1+4]=-dPvz;    A[dim*1+5]=-dPez; 

  A[dim*2+0]=dTvx;     A[dim*2+1]=dTex;     A[dim*2+2]=-dTvy;    A[dim*2+3]=-dTey;    A[dim*2+4]=0.;       A[dim*2+5]=0.; 

  A[dim*3+0]=dTvx;     A[dim*3+1]=dTex;     A[dim*3+2]=0.;       A[dim*3+3]=0.;       A[dim*3+4]=-dTvz;    A[dim*3+5]=-dTez; 

  A[dim*4+0]=lambdaX;  A[dim*4+1]=0.;       A[dim*4+2]=lambdaY;  A[dim*4+3]=0.;       A[dim*4+4]=lambdaZ;  A[dim*4+5]=0.; 

  A[dim*5+0]=0.;       A[dim*5+1]=lambdaX;  A[dim*5+2]=0.;       A[dim*5+3]=lambdaY;  A[dim*5+4]=0.;       A[dim*5+5]=lambdaZ;
  
 
  
  /*
  printf("  * A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      printf("A[%d][%d]= %g  ",i,j, A[dim*i+j] );
    }
    printf("\n");
  }
  printf("\n");  
  */


  double det_a=inverse_matrice_pivot(A, dim, InvA);
  //printf("detA= %g\n",det_a );
  if ( fabs(det_a)<1e-8){printf("A non inversible : detA= %g\n", det_a); return 1;}
  
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
  
  
  // Mutiplcation matricielle de A et InvA  TEST
  /*
  double res;
  //printf("  * InvA*A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      res=0;
      for(int k=0; k<dim; k++){
        res+=InvA[dim*i+k]*A[dim*k+j];
      }
      //printf("A*InvA[%d][%d]= %g ",i,j,res );
    }
    //printf("\n");
  }
  //printf("\n");
  */


  // Mutiplcation matricielle
  
  for(int i=0; i<dim; i++){
    dX[i]=0;
    for(int j=0; j<dim; j++){
      dX[i]+=InvA[dim*i+j]*Delta[j];
    }
    //printf("  dX[%d]= %g\n",i,dX[i] );
  }
  //printf("\n");
  


  double correc=1.;

  *pVX=Vx+correc*dX[0];
  *pEX=Ex+correc*dX[1];
  *pVY=Vy+correc*dX[2];
  *pEY=Ey+correc*dX[3];
  *pVZ=Vz+correc*dX[4];
  *pEZ=Ez+correc*dX[5];

  /*
  printf("  Vx= %.10lf, Ex= %.10lf   Vy= %.10lf, Ey= %.10lf    x= %.10lf\n",Vx,Ex,Vy,Ey,x );
  printf("  (Vmel-Vx)/(Vy-Vx)= %.10lf  ,    (Emel-Ex)/(Ey-Ex)= %.10lf\n", (Vmel-Vx)/(Vy-Vx), (Emel-Ex)/(Ey-Ex));
  printf("  (1-x)*Vx + x*Vy - Vmel= %.10lf  ,    (1-x)*Ex + x*Ey - Emel= %.10lf\n", (1-x)*Vx + x*Vy - Vmel, (1-x)*Ex + x*Ey - Emel);
  */

  return 0;

}

int VE_mixte_cin_3_phase_all_dom(int Nmax, double epsilon, double* lambda, double V0x, double E0x, double V0y, double E0y, double V0z, double E0z, double Vmel, double Emel, double* pVx, double* pEx, double* pVy, double* pEy, double* pVz, double* pEz){
  int err;
  int n=1;
  int dim=6;
   
  // printf(" epsilon= %g\n",epsilon);
  /*
  printf("  phaseX= %d , phaseY= %d\n",phaseX,phaseY );
  printf("  Vmel= %g , Emel= %g\n",Vmel, Emel );
  */
  

  double* A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));
  double* Delta=malloc(dim*sizeof(double));
  double* dX=malloc(dim*sizeof(double));
  
  double Px, Tx;
  double Py, Ty;
  double Pz, Tz;

  double Vx, Ex, Vy, Ey, Vz, Ez;
  double VX, EX, VY, EY, VZ, EZ;

  Vx=V0x; Ex=E0x; 
  Vy=V0y; Ey=E0y; 
  Vz=V0z; Ez=E0z; 

  // phase X
  Px=fP(Vx,Ex, 1); Tx=fT(Vx,Ex, 1);
  // phase Y
  Py=fP(Vy,Ey, 2); Ty=fT(Vy,Ey, 2);
  // phase Z
  Pz=fP(Vz,Ez, 3); Tz=fT(Vz,Ez, 3);


  double critere;
  
  critere=1.;
  //printf("  critere= %.14lf\n",critere );
    
  n=0;
  while(critere>=epsilon && n<Nmax){
    
    err= Newton_cin_3_phase_all_dom(A, InvA, Delta, Vmel, Emel, lambda, Vx, Ex, Vy, Ey, Vz, Ez, &VX, &EX, &VY, &EY, &VZ, &EZ, dX);
    if(err){return 1;}

    Vx=VX; Ex=EX;
    Vy=VY; Ey=EY;
    Vz=VZ; Ez=EZ;
    
    /*
    // phase X
    Px=fP(Vx,Ex, phaseX); Tx=fT(Vx,Ex, phaseX);
    // phase Y
    Py=fP(Vy,Ey, phaseY); Ty=fT(Vy,Ey, phaseY);
    */

    //critere=fabs(Gx-Gy)+fabs(Px-Py)+fabs(Ty-Tx);
    critere=fabs(dX[0]/Vx)+fabs(dX[1]/Ex)+fabs(dX[2]/Vy)+fabs(dX[3]/Ey)+fabs(dX[4]/Vz)+fabs(dX[5]/Ez);
    //critere=fabs(Delta[0])+fabs(Delta[1])+fabs(Delta[2])+fabs(Delta[3]);
    //printf("  n= %d    critere= %g\n",n,critere );

    
    n++;
  }
  //printf(" --n= %d, critere= %g (VE_mixte_cin_3_phase_all_dom)\n",n,critere );

  if(n==Nmax){printf("Newton_mixte_cin_phase non convergent : -n= %d , critere= %g (VE_mixte_cin_phase)\n",n,critere ); return 1;}

  *pVx=Vx; *pEx=Ex;
  *pVy=Vy; *pEy=Ey;
  *pVz=Vz; *pEz=Ez;

  //printf("  ok fin boucle i\n");
  free(A); 
  free(InvA);
  free(Delta);
  free(dX);
  
  return 0;
}




// Renvoie la zone et les fractoins massiques de chaque phases
/*  Arguments d'entrées
     - Xfm : tableau des 3 fractions massiques  

*/
int lec_benchmark_zone(int affichage, double V, double E, int nb_points, double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC, 
                       double* tabVx, double* tabEx, double* tabVy, double* tabEy, int** tabZONE, int* tab_Nb_ZONE, 
                       int Nv, int Ne, double Vmin, double Vmax, double Emin, double Emax, 
                       double* Xfm, int* pzone, double* pV0x, double *pE0x, double* pV0y, double* pE0y){

  int Nmax=3e2;
  double epsilon=2e-8;

  int err; 
  double P,T;
  double Vmx, Emx, Vmy, Emy;
  double Vx, Ex, Vy, Ey;

  double x;

  // nb_points est un argument de la fonction
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3; 
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;

  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );


  // Point triple
  double PT= 4510421781.609182357788085;
  double TT= 624.63982538492518;

  if(V<Vmin || Vmax<V){printf("V hors domaine : Vmin= %g  Vmax= %g  et V= %g (lec_benchmark_zone)\n",Vmin,Vmax,V);   return 1;}  
  if(E<Emin || Emax<E){printf("E hors domaine : Emin= %g  Emax= %g  et E= %g (lec_benchmark_zone)\n",Emin,Emax,E);   return 1;}
 
  double hv=(Vmax-Vmin)/(Nv-1);
  double he=(Emax-Emin)/(Ne-1);
  int i=(V-Vmin)/hv;
  int j=(E-Emin)/he;


  int pos_fic=i*(Ne-1) + j;
  //printf("  --i_dis= %d, j_dis= %d , pos_fic= %d \n",i,j,pos_fic );


  Vmx=tabVx[pos_fic];  Emx=tabEx[pos_fic];
  Vmy=tabVy[pos_fic];  Emy=tabEy[pos_fic];
  //printf("Vmx= %lf , Emx= %lf \n",Vmx,Emx );
  //printf("Vmy= %lf , Emy= %lf \n",Vmy,Emy );
 

  //*pV0x=Vmx; *pE0x=Emx;
  //*pV0y=Vmy; *pE0y=Emy;

  // Valeur d'initialisation
  double VaT=129.905671875844291e-6;  double EaT=72.168190127265447e3;
  double VbT=127.666112185063994e-6;  double EbT=99.480256119236515e3;
  double VcT=132.538478805442338e-6;  double EcT=136.154346228414312e3;
  
  // pour le calcul des fractions massiques dans la zone trinagle --> calcul des coord barycentriques dans le triangle
  double* VE=malloc(2*sizeof(double));
  VE[0]=V; VE[1]=E;
  double** TRI=alloctabd(3,2);
  TRI[0][0]=VaT; TRI[0][1]=EaT;
  TRI[1][0]=VbT; TRI[1][1]=EbT;
  TRI[2][0]=VcT; TRI[2][1]=EcT;
  double* coorBar=malloc(3*sizeof(double));

  int nb_zone=tab_Nb_ZONE[pos_fic];
  
  
  //printf("   nb zone= %d  : ",nb_zone );
  //for(int i=0; i<nb_zone; i++){  printf("  tabZONE[%d][pos_fic]= %d",i,tabZONE[i][pos_fic]);  }  printf("\n");
  

  int zone;
  int phaseA, phaseB;

  
  // On regarde si lorsque qu'il ya 3 zones ou plus, on a -1 dans les zones
  // Dans ce cas on mets comme graines segV[N] 
  // sinon c'est qu'on est au abord du point triple et dans ce cas on met comme graines les sommets du point triple 
  int isout=0;
  if(nb_zone>=3){  
    for(int i=0; i<nb_zone; i++){
      if(tabZONE[i][pos_fic]==-1){isout=1;}
    }
  }
  // on se sert du isout dans la boucle ou nb_zone>=3

  //printf("   nb_zone= %d\n",nb_zone );
  
  //  ##  1 seule zone ##
  if(nb_zone==1){
    zone=tabZONE[0][pos_fic];
    //printf("   - zone= %d\n",zone );

    if(zone==-1){ 
      if(affichage==1){  Xfm[0]=-1; Xfm[1]=-1; Xfm[2]=-1;   *pzone=-1;  return 0;}
      else if(affichage==0){printf("Error out of domain : V= %g , E= %g (lec_benchmark_zone)\n",V,E); return 1;}
      else{printf("Erreur dans la valeur de affichage== %d (lec_benchmark_zone)\n",affichage); return 1;} 
    }
    else if(zone==3 || zone==2 || zone==1){  Xfm[zone-1]=1.;  Xfm[(zone-1+1)%3]=0.;  Xfm[(zone-1+2)%3]=0.;  *pzone=zone;  return 0;  }
    else if(zone==23 || zone==13 || zone==12){
      phaseB=zone%10;  phaseA=(zone-phaseB)/10;
      err=VE_mixte(Nmax, epsilon, phaseA, phaseB, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x);  if(err){return 1;}
      Xfm[0]=0.;  Xfm[1]=0.;  Xfm[2]=0.;
      Xfm[phaseA-1]=1.-x;  Xfm[phaseB-1]=x;
      *pzone=zone;
      return 0;
    }
    else if(zone==0){ 
      // fonction triangle 
      err=coorbary(VE, TRI, &Xfm); if(err){return 1;} 
      *pzone=0;  return 0; 
    }
  }
  //  ##  2 zones   ##   la première zone est forcément la zone pur (ou triple)
  else if(nb_zone==2){
    zone=tabZONE[0][pos_fic];
    //printf("   - zone= %d\n",zone );
    if(zone==3){
      if(intPOLY(V,E,NC,SC)){ Xfm[zone-1]=1.;  Xfm[(zone-1+1)%3]=0.;  Xfm[(zone-1+2)%3]=0.;  *pzone=zone;  return 0;  }
      else if(tabZONE[1][pos_fic]==-1){  
        if(affichage==1){  Xfm[0]=-1; Xfm[1]=-1; Xfm[2]=-1; *pzone=-1;  return 0;}
        else if(affichage==0){printf("Error out of domain : V= %g , E= %g (lec_benchmark_zone)\n",V,E); return 1;}
        else{printf("Erreur dans la valeur de affichage== %d (lec_benchmark_zone)\n",affichage); return 1;} 
      }
      else{
        zone=tabZONE[1][pos_fic];
        phaseB=zone%10;  phaseA=(zone-phaseB)/10;

        err=VE_mixte(Nmax, epsilon, phaseA, phaseB, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x);
        if(err){return 1;}

        Xfm[0]=0.;  Xfm[1]=0.;  Xfm[2]=0.;
        Xfm[phaseA-1]=1.-x;  Xfm[phaseB-1]=x;
        *pzone=zone;
        return 0;
      } 
    }
    else if(zone==2){
      if(intPOLY(V,E,NB,SB)){  Xfm[zone-1]=1.;  Xfm[(zone-1+1)%3]=0.;  Xfm[(zone-1+2)%3]=0.;  *pzone=zone;  return 0;  }
      else if(tabZONE[1][pos_fic]==-1){ 
        if(affichage==1){  Xfm[0]=-1; Xfm[1]=-1; Xfm[2]=-1; *pzone=-1;  return 0;}
        else if(affichage==0){printf("Error out of domain : V= %g , E= %g (lec_benchmark_zone)\n",V,E); return 1;}
        else{printf("Erreur dans la valeur de affichage== %d(lec_benchmark_zone)\n",affichage); return 1;} 
      }
      else{
        zone=tabZONE[1][pos_fic];
        phaseB=zone%10;  phaseA=(zone-phaseB)/10;

        err=VE_mixte(Nmax, epsilon, phaseA, phaseB, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x);
        if(err){return 1;}
        Xfm[0]=0.;  Xfm[1]=0.;  Xfm[2]=0.;
        Xfm[phaseA-1]=1.-x;  Xfm[phaseB-1]=x;
        *pzone=zone;
        
        return 0;
      } 
    }
    else if(zone==1){
      if(intPOLY(V,E,NA,SA)){  Xfm[zone-1]=1.;  Xfm[(zone-1+1)%3]=0.;  Xfm[(zone-1+2)%3]=0.;  *pzone=zone;  return 0;  }
      else if(tabZONE[1][pos_fic]==-1){ 
        if(affichage==1){  Xfm[0]=-1; Xfm[1]=-1; Xfm[2]=-1; *pzone=-1;  return 0;}
        else if(affichage==0){printf("Error out of domain : V= %g , E= %g (lec_benchmark_zone)\n",V,E); return 1;}
        else{printf("Erreur dans la valeur de affichage== %d(lec_benchmark_zone)\n",affichage); return 1;} 
      }
      else{
        zone=tabZONE[1][pos_fic];
        phaseB=zone%10;  phaseA=(zone-phaseB)/10;
 
        err=VE_mixte(Nmax, epsilon, phaseA, phaseB, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x);
        if(err){return 1;}
        Xfm[0]=0.;  Xfm[1]=0.;  Xfm[2]=0.;
        Xfm[phaseA-1]=1.-x;  Xfm[phaseB-1]=x;
        *pzone=zone;

        return 0;
      } 
    }
    else if(zone==0){ 
      if(ftri(V, E)){  err=coorbary(VE, TRI, &Xfm); if(err){return 1;}   *pzone=0;  return 0;   }
      else{
        zone=tabZONE[1][pos_fic];
        phaseB=zone%10;  phaseA=(zone-phaseB)/10;
 
        err=VE_mixte(Nmax, epsilon, phaseA, phaseB, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x);
        if(err){return 1;}
        Xfm[0]=0.;  Xfm[1]=0.;  Xfm[2]=0.;
        Xfm[phaseA-1]=1.-x;  Xfm[phaseB-1]=x;
        *pzone=zone;

        return 0;
      } 
    }
    else if(zone=-1){  
      if(affichage==1){  Xfm[0]=-1; Xfm[1]=-1; Xfm[2]=-1; *pzone=-1;  return 0;}
      else if(affichage==0){printf("Error out of domain : V= %g , E= %g (lec_benchmark_zone)\n",V,E); return 1;}
      else{printf("Erreur dans la valeur de affichage== %d (lec_benchmark_zone)\n",affichage); return 1;} 
    }
    else{printf(" Erreur valeur de zone= %d\n",zone );}
  }
  //  ##  3 zones ou plus  ##
  else{
    for(int ind=0; ind<nb_zone-1; ind++){
      zone=tabZONE[ind][pos_fic];
      //printf("   - zone= %d\n",zone );
      if(zone==3){
        if(intPOLY(V,E,NC,SC)){  Xfm[zone-1]=1.;  Xfm[(zone-1+1)%3]=0.;  Xfm[(zone-1+2)%3]=0.;  *pzone=zone;  return 0;  }
      }
      else if(zone==2){
        if(intPOLY(V,E,NB,SB)){  Xfm[zone-1]=1.;  Xfm[(zone-1+1)%3]=0.;  Xfm[(zone-1+2)%3]=0.;  *pzone=zone;  return 0;  }
      }
      else if(zone==1){
        if(intPOLY(V,E,NA,SA)){  Xfm[zone-1]=1.;  Xfm[(zone-1+1)%3]=0.;  Xfm[(zone-1+2)%3]=0.;  *pzone=zone;  return 0;  }
      }
      else if(zone==23){
        if(intPOLY(V,E,NBC,SBC)){
          if(isout){Vmx=SBC[N23][0]; Emx=SBC[N23][1]; Vmy=SBC[N23-1][0]; Emy=SBC[23-1][1];}
          else{Vmx=VbT; Emx=EbT; Vmy=VcT; Emy=EcT;}
          err=VE_mixte(Nmax, epsilon, 2, 3, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x);
          if(err){return 1;}
          Xfm[0]=0.;  Xfm[1]=1.-x;  Xfm[2]=x;   *pzone=23;

          return 0;
        }
      }
      else if(zone==13){
        if(intPOLY(V,E,NAC,SAC)){
          if(isout){Vmx=SAC[N13][0]; Emx=SAC[N13][1]; Vmy=SAC[N13-1][0]; Emy=SAC[N13-1][1];}
          else{Vmx=VaT; Emx=EaT; Vmy=VcT; Emy=EcT;}
          err=VE_mixte(Nmax, epsilon, 1, 3, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x);
          if(err){return 1;}
          Xfm[0]=1.-x;  Xfm[1]=0;  Xfm[2]=x;   *pzone=13;

          return 0;
        }
      }  
      else if(zone==12){
        if(intPOLY(V,E,NAB,SAB)){
          if(isout){Vmx=SAB[N12][0]; Emx=SAB[N12][1]; Vmy=SAB[N12-1][0]; Emy=SAB[N12-1][1];}
          else{Vmx=VaT; Emx=EaT; Vmy=VbT; Emy=EbT;}
          err=VE_mixte(Nmax, epsilon, 1, 2, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x);
          if(err){return 1;}
          Xfm[0]=1.-x;  Xfm[1]=x;  Xfm[2]=0.;  *pzone=12;

          return 0;
        }
      }
      else if(zone==0){   if(ftri(V, E)){  err=coorbary(VE, TRI, &Xfm); if(err){return 1;}   *pzone=0;  return 0;  }   }  
      else if(zone=-1){  
        if(affichage==1){  Xfm[0]=-1; Xfm[1]=-1; Xfm[2]=-1; *pzone=-1;  return 0;}
        else if(affichage==0){printf("Error out of domain : V= %g , E= %g (lec_benchmark_zone)\n",V,E); return 1;}
        else{printf("Erreur dans la valeur de affichage== %d (lec_benchmark_zone)\n",affichage); return 1;} 
      }
    }
    
    zone=tabZONE[nb_zone-1][pos_fic];
    if(zone==3 || zone==2 || zone==1){  Xfm[zone-1]=1.;  Xfm[(zone-1+1)%3]=0.;  Xfm[(zone-1+2)%3]=0.;  *pzone=zone;  return 0;  }
    else if(zone==23){
      if(isout){Vmx=SBC[N23][0]; Emx=SBC[N23][1]; Vmy=SBC[N23-1][0]; Emy=SBC[23-1][1];}
      else{Vmx=VbT; Emx=EbT; Vmy=VcT; Emy=EcT;}
      err=VE_mixte(Nmax, epsilon, 2, 3, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x);
      if(err){return 1;}
      Xfm[0]=0.;  Xfm[1]=1.0-x;  Xfm[2]=x;  *pzone=23;

      return 0;
    }
    else if(zone==13){
      if(isout){Vmx=SAC[N13][0]; Emx=SAC[N13][1]; Vmy=SAC[N13-1][0]; Emy=SAC[N13-1][1];}
      else{Vmx=VaT; Emx=EaT; Vmy=VcT; Emy=EcT;}
      err=VE_mixte(Nmax, epsilon, 1, 3, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x);
      if(err){return 1;}
      Xfm[0]=1.-x;  Xfm[1]=0.;  Xfm[2]=x;  *pzone=13;

      return 0;
    }  
    else if(zone==12){
      if(isout){Vmx=SAB[N12][0]; Emx=SAB[N12][1]; Vmy=SAB[N12-1][0]; Emy=SAB[N12-1][1];}
      else{Vmx=VaT; Emx=EaT; Vmy=VbT; Emy=EbT;}
      err=VE_mixte(Nmax, epsilon, 1, 2, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x);
      if(err){return 1;}
      Xfm[0]=1.-x;  Xfm[1]=x;  Xfm[2]=0.;   *pzone=12;

      return 0;
    }
    else if(zone==0){  err=coorbary(VE, TRI, &Xfm); if(err){return 1;}   *pzone=0;  return 0;  }
    else if(zone=-1){  
      if(affichage==1){  Xfm[0]=-1; Xfm[1]=-1; Xfm[2]=-1; *pzone=-1;  return 0;}
      else if(affichage==0){printf("Error out of domain : V= %g , E= %g (lec_benchmark_zone)\n",V,E); return 1;}
      else{printf("Erreur dans la valeur de affichage== %d (lec_benchmark_zone)\n",affichage); return 1;} 
    }
    else{printf("Erreur dans tabZONE[%d][%d]= %d (lec_benchmark_zone)\n\n",nb_zone-1,pos_fic,tabZONE[nb_zone-1][pos_fic]); return 1;}
    
  }
  // Si on arrive la c'est que la bonne zone n'a pas été trouvé
  printf("Erreur choix zone (lec_benchmark_zone)\n");
  printf("  nb_zone= %d",nb_zone );
  printf("  tabZONE[0][%d]= %d",pos_fic,tabZONE[0][pos_fic]);
  for(int ind=1; ind<nb_zone; ind++){
    printf("  tabZONE[%d][%d]= %d",ind,pos_fic,tabZONE[ind][pos_fic]); 
  }
  printf("\n V= %.10lf , E= %.10lf\n",V,E);
  return 1; 
  
}


// Calcul de la célérite du son
/*  Pp=P(V+dV,E+dE(+dV) )  
    Pm=P(V-dV,E+dE(-dV) )  

    Argument de sortie :kugh
     - double* pc : pointeur sur la celerité du son
*/
int c_celerite_son_cin_phase(int Nmax, double epsilon, int zone, double P, double V, double E, double dV, 
                             double V0x, double E0x, double V0y, double E0y, double* lambda, double *pc){

  int err;

  double Pp,Pm;
  double dE=-P*dV;

  double Vx,Ex,Vy,Ey,x;
  double lambdaB;

  int phaseA,phaseB;


  if(zone==-1){printf("Erreur out of domain : zone= %d (c_celerite_son_cin_phase)\n",zone); return 1;  }
  else if(zone==1 || zone==2 || zone==3){  
    Pp=fP(V+dV,E+dE,zone); 
    Pm=fP(V-dV,E-dE,zone); 
  }
  else if(zone==12 || zone==13 || zone==23){
    phaseB=zone%10;
    phaseA=(zone-phaseB)/10;
    lambdaB=lambda[phaseB-1];

    err=VE_mixte_cin_phase(Nmax, epsilon, phaseA, phaseB, lambdaB, V0x, E0x, V0y, E0y, V+dV, E+dE, &Vx, &Ex, &Vy, &Ey);  if(err){return 1;}
    //err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V+dV, E+dE, &Vx, &Ex, &Vy, &Ey, &x);  if(err){return 1;}
    Pp=fP(Vx,Ex,phaseA);

    err=VE_mixte_cin_phase(Nmax, epsilon, phaseA, phaseB, lambdaB, V0x, E0x, V0y, E0y, V-dV, E-dE, &Vx, &Ex, &Vy, &Ey);  if(err){return 1;}
    //err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V-dV, E-dE, &Vx, &Ex, &Vy, &Ey, &x);  if(err){return 1;}
    Pm=fP(Vx,Ex,phaseA);
  }
  else if(zone==0){
    *pc=0;
    return 0;
  }

  double dPpdV=(Pp-P)/dV;
  double dPmdV=(Pm-P)/(-dV);
  double dPdV=(dPpdV + dPmdV)/2;
  double c2=-V*V*dPdV; 

  
  if(c2<0){
    printf("Erreur c2= %g negatif\n",c2);
    printf("  V= %g , E= %g , zone= %d \n",V,E,zone );  
    printf("  dV/V= %g , dE/E= %g \n",dV/V,dE/E );
    printf("  P= %g , Pp= %g , Pm= %g\n",P,Pp,Pm );
    printf("  dPp/dV= %g , dPm/dV= %g , dPdV= %g \n",dPpdV,dPmdV,dPdV );
    printf("  lambdaB= %g , c^2= %g \n",lambdaB, c2 );
    //*pc=3.3e2;
    return 1;
  }
  *pc=sqrt(c2);
  return 0;
}

int c_celerite_son_cin_3_phase_all_dom(int Nmax, double epsilon, int zone, double P, double V, double E, double dV, 
                             double Va0, double Ea0, double Vb0, double Eb0, double Vc0, double Ec0, double* lambda, double *pc){

  int err;

  double Pp,Pm;
  double dE=-P*dV;

  double Vx,Ex,Vy,Ey,Vz,Ez;



  err= VE_mixte_cin_3_phase_all_dom(Nmax, epsilon, lambda, Va0, Ea0, Vb0, Eb0, Vc0, Ec0, V+dV, E+dE, &Vx, &Ex, &Vy, &Ey, &Vz, &Ez);
  if(err){printf("erreur VE_mixte_cin_3_phase_all_dom\n"); return 1;}
  Pp=fP(Vx,Ex,1);

  err= VE_mixte_cin_3_phase_all_dom(Nmax, epsilon, lambda, Va0, Ea0, Vb0, Eb0, Vc0, Ec0, V-dV, E-dE, &Vx, &Ex, &Vy, &Ey, &Vz, &Ez);
  if(err){printf("erreur VE_mixte_cin_3_phase_all_dom\n"); return 1;}
  Pm=fP(Vx,Ex,1);

  double dPpdV=(Pp-P)/dV;
  double dPmdV=(Pm-P)/(-dV);
  double dPdV=(dPpdV + dPmdV)/2;
  double c2=-V*V*dPdV; 

  
  if(c2<0){
    printf("Erreur c2= %g negatif\n",c2);
    printf("  V= %g , E= %g , zone= %d \n",V,E,zone );  
    printf("  dV/V= %g , dE/E= %g \n",dV/V,dE/E );
    printf("  P= %g , Pp= %g , Pm= %g\n",P,Pp,Pm );
    printf("  dPp/dV= %g , dPm/dV= %g , dPdV= %g \n",dPpdV,dPmdV,dPdV );
    printf("  lambdaA= %g , lambdaB= %g , lambdaC= %g , c^2= %g \n",lambda[0],lambda[1],lambda[2], c2 );
    //*pc=3.3e2;
    return 1;
  }
  *pc=sqrt(c2);
  return 0;
}

// Renvoie la valeur de P,T et C pour un couple V,E en prenant en compte la cinétique de chgmt de phase
/*  Argments d'entrees :
      - CIN_PAHSE :  if CIN_PHASE==0 : pas de cinétique de phase
                     if CIN_PHASE==1 : prise en compte de la cinétique de phase

   Metre un param entier tels que 
    meth==0 : calcul + tabulation
    meth==1 : calcul pure
    meth==2 : tabulation pure
    meth==3 : maillage ??

    Pour le moment uniquement avec cal_tab
*/
int fPTC_cin_phase_old_wrong(int affichage, int CIN_PHASE, double V, double E, double* pP, double *pc, double* pT, double* lambda_in, double* lambda_out, int Ntau, double dx, double dt, double Ttau, int* pzone,
                   int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, 
                   double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC, 
                   double* tabVx, double* tabEx, double* tabVy, double* tabEy, int** tabZONE, int* tab_Nb_ZONE ){
  
  //      Ttau est calculé avec c(Vn, En)   
  
  if(CIN_PHASE==2){  Ttau=Ntau*dx/(*pc);   }
  //double Ttau=1e-9;   //Ntau*dx/(*pc);


  int Nmax=3e2;
  double epsilon=2e-9;
  int err;
  // Pour le calcul de la celerité du son
  double V0x, E0x, V0y, E0y;  
  double Vx, Ex, Vy, Ey;  
  double dV=1e-8;

  int zone;

  double* Xfm=malloc(3*sizeof(double)); // tableau des fractions massiques
  int phaseX,phaseY;

  // Point triple
  double PT= 4510421781.609182357788085;
  double TT= 624.63982538492518;

  //err= PT_tab_cal(V, E, nb_points, SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE, 
  //                Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, pP, pT, px, &zone, &V0x, &E0x, &V0y, &E0y);

  // Focntion renvoie la zone, la fraction massique et les graines pour les Newton_mixte
  
  err= lec_benchmark_zone(affichage, V, E, nb_points, SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE, 
                          Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, Xfm, &zone, &V0x, &E0x, &V0y, &E0y);
  if(err){return err;}
   
  *pzone=zone;

  //printf(" ZONE= %d\n",zone );
  
  //for(int i=0; i<3; i++){ if(Xfm[i]<0 ){Xfm[i]=0;}  if(Xfm[i]>1 ){Xfm[i]=1;} }
  
  //printf("  V= %g , E= %g , zone= %d\n",V,E,zone );
  //printf(" V= %g, E= %g, Xfm[0]= %g, Xfm[2]= %g, Xfm[2]= %g \n",V,E,Xfm[0],Xfm[1],Xfm[2] );  

  if(zone==-1){
    if(affichage==1){*pP=-10e9; *pT=0.;  lambda_out[0]=-1; lambda_out[1]=-1; lambda_out[2]=-1;  return 0;}
    else if(affichage==0){printf("Error out of domain : V= %g , E= %g\n",V,E); return 1;}
    else{printf("Erreur affichage= %d",affichage); return 1;} 
  }
  if(CIN_PHASE==0){
    for(int i=0; i<3; i++){  lambda_out[i]=Xfm[i];  }
  }
  else if(CIN_PHASE==1 || CIN_PHASE==2){
    // Calcul de x^{n+1} 
    if(zone==1 || zone==2 || zone==3){ for(int i=0; i<3; i++){  lambda_out[i]=Xfm[i]; } }
    else if(zone==0){
      //printf("DlambdaA= %g , DlambdaB= %g \n", (Ttau*(lambda[0])+Xfm[0])/(Ttau+1),(Ttau*(lambda[1])+Xfm[1])/(Ttau+1));
      lambda_out[0]=((dt/Ttau)*Xfm[0]+lambda_in[0])/((dt/Ttau)+1.);
      lambda_out[1]=((dt/Ttau)*Xfm[1]+lambda_in[1])/((dt/Ttau)+1.);
      lambda_out[2]=1.-lambda_out[1]-lambda_out[0];

    }
    else if(zone==12 ){
      //printf("DlambdaA= %g \n", (Ttau*(lambda[0])+Xfm[0])/(Ttau+1));
      lambda_out[0]=((dt/Ttau)*Xfm[0]+lambda_in[0])/((dt/Ttau)+1.);
      lambda_out[1]=1.-lambda_out[0]; // (autre phase)
      lambda_out[2]=0.;
      //printf(" ok 12  lambda[0]= %g , Xfm[0]= %g\n",lambda[0],Xfm[0]);
    }
    else if( zone==13 ){
      //printf("DlambdaA= %g \n", (Ttau*(lambda[0])+Xfm[0])/(Ttau+1));
      lambda_out[0]=((dt/Ttau)*Xfm[0]+lambda_in[0])/((dt/Ttau)+1.);
      lambda_out[1]=0;
      lambda_out[2]=1.-lambda_out[0]; // (autre phase)
      
    }
    else if( zone==23 ){
      //printf("DlambdaB= %g \n", (Ttau*(lambda[0])+Xfm[0])/(Ttau+1));
      lambda_out[0]=0;
      lambda_out[1]=((dt/Ttau)*Xfm[1]+lambda_in[1])/((dt/Ttau)+1.);
      lambda_out[2]=1.-lambda_out[1]; // (autre phase)
      
    }
  }
  else{printf("Erreur CIN_PHASE= %d\n",CIN_PHASE ); return 1;}


  //printf(" ZONE= %d | In : xA= %g , xB= %g, xC= %g   Out : xA= %g , xB= %g, xC= %g\n",zone,lambda_in[0],lambda_in[1],lambda_in[2],lambda_out[0],lambda_out[1],lambda_out[2] );  


  // HYDRO
  //printf("   debut hydro  :  zone= %d\n",zone);

  if(zone==0){  
    *pP=PT;
    *pT=TT;
    *pc=0;

    return 0;
  }
  else if(zone==12 || zone==13 || zone==23){
    phaseY=zone%10;
    phaseX=(zone-phaseY)/10;
    // Calcul de p,T et C grace aux nouvelles valeurs de fractions massiques
    // ON IMPOSE Lambda MAIS IL A DEJA ETE CALCULER DANS LE CALCUL DE LA ZONE DONC OK
    //lambda[phaseY-1]*=(1+5e-1);
    //if(lambda[phaseY-1]>1){ lambda[phaseY-1]=1;}
    err= VE_mixte_cin_phase(Nmax, epsilon, phaseX, phaseY, lambda_out[phaseY-1], V0x, E0x, V0y, E0y, V, E, &Vx, &Ex, &Vy, &Ey);
    if(err){return 1;}
    *pP=fP(Vx,Ex,phaseX);
    *pT=fT(Vx,Ex,phaseX);
  }
  else if(zone==1 || zone==2 || zone==3){
    *pP=fP(V,E,zone);
    *pT=fT(V,E,zone);
  }

  /*
  if(*pP>1e11 || *pT<-1e3){
    printf(" Vx= %g , Ex= %g , Vy= %g , Ey= %g\n",Vx,Ex,Vy,Ey);
    printf(" V= %g , E= %g , zone= %d, xA= %g , xB= %g , xC= %g\n",V,E,zone,lambda[0],lambda[01],lambda[2]);
    printf(" Px= %g , Tx= %g , Py= %g , Ty= %g \n",fP(Vx,Ex,phaseX),fT(Vx,Ex,phaseX),fP(Vy,Ey,phaseY),fT(Vy,Ey,phaseY));
    return 1;
  }
  */

  // Problème de cérité du son , 
  // On ne peut calculer P de la même manière, il faut utiliser VE_mixte_cin_phase si on est en zone de mélange mixte
  // Faire une nouvelle fonction celerite avec un test dedans
  err=c_celerite_son_cin_phase(Nmax, epsilon, zone, *pP, V, E, dV, V0x, E0x, V0y, E0y, lambda_out, pc);  
  if(err){return 1;}

  return 0;
}
/*  Implemenation du schema temporelle des fraction massiques
  

*/

int lambda_cin_phase(int sch_cin_phase, double* lambda_out, double* lambda_in, double* lambda_eq, double dt, double Ttau ){
  int sgn;
  if(sch_cin_phase==0){ // Schéma implcite d'odre 1 classique,
    lambda_out[0]=((dt/Ttau)*lambda_eq[0]+lambda_in[0])/((dt/Ttau)+1.);
    lambda_out[1]=((dt/Ttau)*lambda_eq[1]+lambda_in[1])/((dt/Ttau)+1.);
    lambda_out[2]=((dt/Ttau)*lambda_eq[2]+lambda_in[2])/((dt/Ttau)+1.);
  }
  //else if(sch==1){ // Schéma implcite d'odre 1 Fickett-Rivard

  return 0;
}



int fPTC_cin_phase_old_ok(int affichage, int CIN_PHASE, double V, double E, double* pP, double *pc, double* pT, double* lambda_in, double* lambda_out, int Ntau, double dx, double dt, double Ttau, int* pzone,
                   int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, 
                   double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC, 
                   double* tabVx, double* tabEx, double* tabVy, double* tabEy, int** tabZONE, int* tab_Nb_ZONE ){
  
  //      Ttau est calculé avec c(Vn, En)   
  
  if(CIN_PHASE==2){  Ttau=Ntau*dx/(*pc);   }
  //double Ttau=1e-9;   //Ntau*dx/(*pc);


  int Nmax=3e2;
  double epsilon=2e-9;
  int err;
  // Pour le calcul de la celerité du son
  double V0x, E0x, V0y, E0y;  
  double Vx, Ex, Vy, Ey;  
  double dV=1e-8;

  int zone;

  double* Xfm=malloc(3*sizeof(double)); // tableau des fractions massiques
  int phaseX,phaseY;

  // Point triple
  double PT= 4510421781.609182357788085;
  double TT= 624.63982538492518;

  double Va0, Ea0, Vb0, Eb0, Vc0, Ec0; 
  Va0=125e-6;  Ea0=2.5e4;
  Vb0=108e-6;  Eb0=3.4e5;
  Vc0=112e-6;  Ec0=6.8e5;
  double facINIT=1.1;

  //err= PT_tab_cal(V, E, nb_points, SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE, 
  //                Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, pP, pT, px, &zone, &V0x, &E0x, &V0y, &E0y);

  // Focntion renvoie la zone, la fraction massique et les graines pour les Newton_mixte
  err= lec_benchmark_zone(affichage, V, E, nb_points, SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE, 
                          Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, Xfm, &zone, &V0x, &E0x, &V0y, &E0y);
  if(err){return err;}
  
  *pzone=zone;

  //printf(" ZONE= %d\n",zone );
  
  //for(int i=0; i<3; i++){ if(Xfm[i]<0 ){Xfm[i]=0;}  if(Xfm[i]>1 ){Xfm[i]=1;} }
  
  //printf("  V= %g , E= %g , zone= %d\n",V,E,zone );
  //printf(" V= %g, E= %g, Xfm[0]= %g, Xfm[2]= %g, Xfm[2]= %g \n",V,E,Xfm[0],Xfm[1],Xfm[2] );  
  
  if(zone==-1){
    if(affichage==1){*pP=-10e9; *pT=0.;  lambda_out[0]=-1; lambda_out[1]=-1; lambda_out[2]=-1;  return 0;}
    else if(affichage==0){printf("Error out of domain : V= %g , E= %g\n",V,E); return 1;}
    else{printf("Erreur affichage= %d",affichage); return 1;} 
  }
  if(CIN_PHASE==0){
    for(int i=0; i<3; i++){  lambda_out[i]=Xfm[i];  }
  }
  else if(CIN_PHASE==1 || CIN_PHASE==2){
    //printf("DlambdaA= %g , DlambdaB= %g \n", (Ttau*(lambda[0])+Xfm[0])/(Ttau+1),(Ttau*(lambda[1])+Xfm[1])/(Ttau+1));

    //int lambda_cin_phase(int sch, double* lambda, double* lambda_eq, double dt, double Ttau );
    lambda_out[0]=((dt/Ttau)*Xfm[0]+lambda_in[0])/((dt/Ttau)+1.);
    lambda_out[1]=((dt/Ttau)*Xfm[1]+lambda_in[1])/((dt/Ttau)+1.);
    lambda_out[2]=((dt/Ttau)*Xfm[2]+lambda_in[2])/((dt/Ttau)+1.);

  }
  else{printf("Erreur CIN_PHASE= %d\n",CIN_PHASE ); return 1;}
  

  //printf(" ZONE= %d | xeqA= %g , xeqB= %g, xeqC= %g   Out : xA= %g , xB= %g, xC= %g\n",zone,Xfm[0],Xfm[1],Xfm[2],lambda_out[0],lambda_out[1],lambda_out[2] );  


  // HYDRO
  //printf("   debut hydro  :  zone= %d\n",zone);
  if(CIN_PHASE==0){
    if(zone==0){  
      *pP=PT;
      *pT=TT;
      *pc=0;

      return 0;
    }
    else if(zone==12 || zone==13 || zone==23){
      phaseY=zone%10;
      phaseX=(zone-phaseY)/10;
      // Calcul de p,T et C grace aux nouvelles valeurs de fractions massiques
      // ON IMPOSE Lambda MAIS IL A DEJA ETE CALCULER DANS LE CALCUL DE LA ZONE DONC OK
      //lambda[phaseY-1]*=(1+5e-1);
      //if(lambda[phaseY-1]>1){ lambda[phaseY-1]=1;}
      err= VE_mixte_cin_phase(Nmax, epsilon, phaseX, phaseY, lambda_out[phaseY-1], V0x, E0x, V0y, E0y, V, E, &Vx, &Ex, &Vy, &Ey);
      if(err){return 1;}
      *pP=fP(Vx,Ex,phaseX);
      *pT=fT(Vx,Ex,phaseX);
    }
    else if(zone==1 || zone==2 || zone==3){
      *pP=fP(V,E,zone);
      *pT=fT(V,E,zone);
    }

  }
  else if(CIN_PHASE!=0){
    if(zone==0){  // ????????   RESTE A TRAITER
      *pP=PT;
      *pT=TT;
      *pc=0;

      return 0;
    }
    else if(zone==12 || zone==13 || zone==23){
      phaseY=zone%10;
      phaseX=(zone-phaseY)/10;
      // Calcul de p,T et C grace aux nouvelles valeurs de fractions massiques
      // ON IMPOSE Lambda MAIS IL A DEJA ETE CALCULER DANS LE CALCUL DE LA ZONE DONC OK
      //lambda[phaseY-1]*=(1+5e-1);
      //if(lambda[phaseY-1]>1){ lambda[phaseY-1]=1;}
      err= VE_mixte_cin_phase(Nmax, epsilon, phaseX, phaseY, lambda_out[phaseY-1], V0x, E0x, V0y, E0y, V, E, &Vx, &Ex, &Vy, &Ey);
      if(err){return 1;}
      *pP=fP(Vx,Ex,phaseX);
      *pT=fT(Vx,Ex,phaseX);
    }

    else if(zone==1 || zone==2 || zone==3){
      if(zone==1 && lambda_out[1]==0 && lambda_out[2]==0){
        *pP=fP(V,E,zone);
        *pT=fT(V,E,zone);
      }
      else if(zone==2 && lambda_out[0]==0 && lambda_out[2]==0){
        *pP=fP(V,E,zone);
        *pT=fT(V,E,zone);
      }
      else if(zone==3 && lambda_out[0]==0 && lambda_out[1]==0){
        *pP=fP(V,E,zone);
        *pT=fT(V,E,zone);
      }
      else if(lambda_out[0]==0 || lambda_out[2]>0){
        if(zone==3){  V0y=facINIT*V; E0y=facINIT*E; V0x=Vb0; E0x=Eb0; }
        else if(zone==2){  Vx=V; Ex=E; Vy=Vc0; Ey=Ec0; }
        else if(zone==1){  printf("Erreur zone= %d, et lambdaA= %g\n",zone, lambda_out[0] ); }
        phaseY=3;
        phaseX=2;

        err= VE_mixte_cin_phase(Nmax, epsilon, phaseX, phaseY, lambda_out[phaseY-1], V0x, E0x, V0y, E0y, V, E, &Vx, &Ex, &Vy, &Ey);
        if(err){return 1;}

        *pP=fP(Vx,Ex,phaseX);
        *pT=fT(Vx,Ex,phaseX);
      }
      else if(lambda_out[1]==0){
        if(zone==3){  V0y=facINIT*V; E0y=facINIT*E; V0x=Va0; E0x=Ea0; }
        else if(zone==1){  Vx=V; Ex=E; Vy=Vc0; Ey=Ec0; }
        else if(zone==2){  printf("Erreur zone= %d, et lambdaA= %g\n",zone, lambda_out[1] ); }
        phaseY=3;
        phaseX=1;

        err= VE_mixte_cin_phase(Nmax, epsilon, phaseX, phaseY, lambda_out[phaseY-1], V0x, E0x, V0y, E0y, V, E, &Vx, &Ex, &Vy, &Ey);
        if(err){return 1;}

        *pP=fP(Vx,Ex,phaseX);
        *pT=fT(Vx,Ex,phaseX);

      }
      else if(lambda_out[2]==0){
        if(zone==2){  V0y=facINIT*V; E0y=facINIT*E; V0x=Va0; E0x=Ea0; }
        else if(zone==1){  Vx=V; Ex=E; Vy=Vb0; Ey=Eb0; }
        else if(zone==3){  printf("Erreur zone= %d, et lambdaA= %g\n",zone, lambda_out[2] ); }
        phaseY=2;
        phaseX=1;

        err= VE_mixte_cin_phase(Nmax, epsilon, phaseX, phaseY, lambda_out[phaseY-1], V0x, E0x, V0y, E0y, V, E, &Vx, &Ex, &Vy, &Ey);
        if(err){return 1;}

        *pP=fP(Vx,Ex,phaseX);
        *pT=fT(Vx,Ex,phaseX);
      }
      else{printf(" Erreur zone= %d, xA= %g, xB= %g, xC= %g\n",zone, lambda_out[0],  lambda_out[1], lambda_out[2]);}
      

    }
  }

  /*
  if(*pP>1e11 || *pT<-1e3){
    printf(" Vx= %g , Ex= %g , Vy= %g , Ey= %g\n",Vx,Ex,Vy,Ey);
    printf(" V= %g , E= %g , zone= %d, xA= %g , xB= %g , xC= %g\n",V,E,zone,lambda[0],lambda[01],lambda[2]);
    printf(" Px= %g , Tx= %g , Py= %g , Ty= %g \n",fP(Vx,Ex,phaseX),fT(Vx,Ex,phaseX),fP(Vy,Ey,phaseY),fT(Vy,Ey,phaseY));
    return 1;
  }
  */

  // Problème de cérité du son , 
  // On ne peut calculer P de la même manière, il faut utiliser VE_mixte_cin_phase si on est en zone de mélange mixte
  // Faire une nouvelle fonction celerite avec un test dedans
  err=c_celerite_son_cin_phase(Nmax, epsilon, zone, *pP, V, E, dV, V0x, E0x, V0y, E0y, lambda_out, pc);  
  if(err){return 1;}

  return 0;
}


int fPTC_cin_phase(int affichage, int CIN_PHASE, int sch_cin_phase, double V, double E, double* pP, double *pc, double* pT, double* lambda_in, double* lambda_out, int Ntau, double dx, double dt, double Ttau, int* pzone,
                   int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, 
                   double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC, 
                   double* tabVx, double* tabEx, double* tabVy, double* tabEy, int** tabZONE, int* tab_Nb_ZONE ){
  
  //      Ttau est calculé avec c(Vn, En)   
  
  if(CIN_PHASE==2){  Ttau=Ntau*dx/(*pc);   }
  //double Ttau=1e-9;   //Ntau*dx/(*pc);


  int Nmax=3e2;
  double epsilon=1e-14;
  int err;
  // Pour le calcul de la celerité du son
  double V0x, E0x, V0y, E0y, V0z, E0z;  
  double Vx, Ex, Vy, Ey, Vz, Ez;  
  double dV=1e-8;

  int zone;

  double* Xfm=malloc(3*sizeof(double)); // tableau des fractions massiques
  int phaseX,phaseY;

  // Point triple
  double PT= 4510421781.609182357788085;
  double TT= 624.63982538492518;

  double Va0, Ea0, Vb0, Eb0, Vc0, Ec0; 
  Va0=125e-6;  Ea0=2.5e4;
  Vb0=108e-6;  Eb0=3.4e5;
  Vc0=112e-6;  Ec0=6.8e5;
  double facINIT=1.1;

  //err= PT_tab_cal(V, E, nb_points, SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE, 
  //                Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, pP, pT, px, &zone, &V0x, &E0x, &V0y, &E0y);

  // Focntion renvoie la zone, la fraction massique et les graines pour les Newton_mixte
  err= lec_benchmark_zone(affichage, V, E, nb_points, SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE, 
                          Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, Xfm, &zone, &V0x, &E0x, &V0y, &E0y);
  if(err){return err;}
  
  *pzone=zone;

  //printf(" ZONE= %d\n",zone );
  
  for(int i=0; i<3; i++){ if(Xfm[i]<0 ){Xfm[i]=0;}  if(Xfm[i]>1 ){Xfm[i]=1;} }
  
  //printf("  V= %g , E= %g , zone= %d, cin_phase= %d \n",V,E,zone,CIN_PHASE );
  //printf(" V= %g, E= %g, Xfm[0]= %g, Xfm[2]= %g, Xfm[2]= %g \n",V,E,Xfm[0],Xfm[1],Xfm[2] );  
  
  if(zone==-1){
    if(affichage==1){*pP=-10e9; *pT=0.;  lambda_out[0]=-1; lambda_out[1]=-1; lambda_out[2]=-1;  return 0;}
    else if(affichage==0){printf("Error out of domain : V= %g , E= %g\n",V,E); return 1;}
    else{printf("Erreur affichage= %d",affichage); return 1;} 
  }
  if(CIN_PHASE==0){ 
    for(int i=0; i<3; i++){  lambda_out[i]=Xfm[i];  }
  }
  else if(CIN_PHASE==1 || CIN_PHASE==2){
    //printf("DlambdaA= %g , DlambdaB= %g \n", (Ttau*(lambda[0])+Xfm[0])/(Ttau+1),(Ttau*(lambda[1])+Xfm[1])/(Ttau+1));

    //int lambda_cin_phase(int sch, double* lambda, double* lambda_eq, double dt, double Ttau );
    /*
    lambda_out[0]=((dt/Ttau)*Xfm[0]+lambda_in[0])/((dt/Ttau)+1.);
    lambda_out[1]=((dt/Ttau)*Xfm[1]+lambda_in[1])/((dt/Ttau)+1.);
    lambda_out[2]=((dt/Ttau)*Xfm[2]+lambda_in[2])/((dt/Ttau)+1.);
    */
    err= lambda_cin_phase(sch_cin_phase, lambda_out, lambda_in, Xfm, dt, Ttau );
    if (err){ printf("Erreur lambda_cin_phase\n"); return 1; }
  }
  else{printf("Erreur CIN_PHASE= %d\n",CIN_PHASE ); return 1;}
  

  //printf(" ZONE= %d | xeqA= %g , xeqB= %g, xeqC= %g   Out : xA= %g , xB= %g, xC= %g\n",zone,Xfm[0],Xfm[1],Xfm[2],lambda_out[0],lambda_out[1],lambda_out[2] );  


  // HYDRO
  //printf("   debut hydro  :  zone= %d\n",zone);
  
  if(CIN_PHASE==0){
    if(zone==0){  
      *pP=PT;
      *pT=TT;
      *pc=0;

      return 0;
    }
    else if(zone==12 || zone==13 || zone==23){
      phaseY=zone%10;
      phaseX=(zone-phaseY)/10;
      // Calcul de p,T et C grace aux nouvelles valeurs de fractions massiques
      // ON IMPOSE Lambda MAIS IL A DEJA ETE CALCULER DANS LE CALCUL DE LA ZONE DONC OK
      //lambda[phaseY-1]*=(1+5e-1);
      //if(lambda[phaseY-1]>1){ lambda[phaseY-1]=1;}
      err= VE_mixte_cin_phase(Nmax, epsilon, phaseX, phaseY, lambda_out[phaseY-1], V0x, E0x, V0y, E0y, V, E, &Vx, &Ex, &Vy, &Ey);
      if(err){return 1;}
      *pP=fP(Vx,Ex,phaseX);
      *pT=fT(Vx,Ex,phaseX);
    }
    else if(zone==1 || zone==2 || zone==3){
      *pP=fP(V,E,zone);
      *pT=fT(V,E,zone);
    }

  }
  else if(CIN_PHASE!=0){
    if(zone==0){  
      *pP=PT;
      *pT=TT;
      *pc=0;

      return 0;
    }
    else {

      err= VE_mixte_cin_3_phase_all_dom(Nmax, epsilon, lambda_out, Va0, Ea0, Vb0, Eb0, Vc0, Ec0, V, E, &Vx, &Ex, &Vy, &Ey, &Vz, &Ez);
      if(err){printf("erreur VE_mixte_cin_3_phase_all_dom\n"); return 1;}
      *pP=fP(Vx,Ex,1);
      *pT=fT(Vx,Ex,1);
    }
  }

  /*
  if(*pP>1e11 || *pT<-1e3){
    printf(" Vx= %g , Ex= %g , Vy= %g , Ey= %g\n",Vx,Ex,Vy,Ey);
    printf(" V= %g , E= %g , zone= %d, xA= %g , xB= %g , xC= %g\n",V,E,zone,lambda[0],lambda[01],lambda[2]);
    printf(" Px= %g , Tx= %g , Py= %g , Ty= %g \n",fP(Vx,Ex,phaseX),fT(Vx,Ex,phaseX),fP(Vy,Ey,phaseY),fT(Vy,Ey,phaseY));
    return 1;
  }
  */

  // Problème de cérité du son , 
  // On ne peut calculer P de la même manière, il faut utiliser VE_mixte_cin_phase si on est en zone de mélange mixte
  // Faire une nouvelle fonction celerite avec un test dedans
  err=c_celerite_son_cin_3_phase_all_dom(Nmax, epsilon, zone, *pP, V, E, dV, Va0, Ea0, Vb0, Eb0, Vc0, Ec0, lambda_out, pc);  
  if(err){return 1;}

  return 0;
}


//   fonction ENERGIE INTERNE SPECIFIQUE epsilon(tau, p)
// Alogorthime de Newton pour calculer EPSn+1 sans conaitre Pn+1
// Necessaire pour le schema vNR
int fPn1_etain_cin_phase(int tst, int cin_phase, int sch_cin_phase, double TAUn, double TAUn1, double Pn, double Qn, double EPSn, double* pEPSn1, double* pPn1, double* lambda_in, int Ntau, double dX, double dt, double Ttau,
                         int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, 
                         double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC, 
                         double* tabVx, double* tabEx, double* tabVy, double* tabEy, int** tabZONE, int* tab_Nb_ZONE ){

  int Nmax=6e2;
  double eps;
  if(tst/100==1){  eps=1e-12;  }
  else if(tst/100==0){  eps=1e-10;  }

  double epsilon=1e-7;

  int err;
  double res;
  double dTAU=TAUn1-TAUn; 
  double dE;
  double E=EPSn;
  double C;

  //printf(" V= %g , E= %g , TAUn1= %g , Pn= %g , xA= %g , xB= %g , xC= %g\n",TAUn,EPSn,TAUn1,Pn,lambda[0],lambda[1],lambda[2] );
  
  // Calcul de la derivée
  double dE_inf=EPSn*1e-5;
  double P,T,x;
  int zone;
  double Pp,Pm;
  int phaseA,phaseB;
  double dPpdE, dPmdE, dPdE;
  double V0x, E0x, V0y, E0y, Vx, Ex, Vy, Ey;

  double hv=(Vmax_dis-Vmin_dis)/(Nv_dis-1);
  double he=(Emax_dis-Emin_dis)/(Ne_dis-1);
  int i=(TAUn1-Vmin_dis)/hv;
  int j=(EPSn-Emin_dis)/he;

  int pos_fic=i*(Ne_dis-1) + j;
  //printf("  --i_dis= %d, j_dis= %d \n",i,j );

  V0x=tabVx[pos_fic];  E0x=tabEx[pos_fic];
  V0y=tabVy[pos_fic];  E0y=tabEy[pos_fic];

  // Calcul de la derivee
  double lambdaB;
  double* lambda_out=malloc(3*sizeof(double)); // récupére la modif de lambda dans fPTC_cin_phase (inutile)
  double* lambda_test=malloc(3*sizeof(double)); // récupére la modif de lambda dans fPTC_cin_phase (inutile)
  lambda_test[0]=lambda_in[0];
  lambda_test[1]=lambda_in[1];
  lambda_test[2]=lambda_in[2];


  double f,f_prime;
  double critere=1.;
  int n=0;
  while(critere>eps && n<Nmax){
    ///////// Calcul de la dérive de P(.,TAUn)
    /*
    err= lec_benchmark_PT(0,TAUn1, E, nb_points, SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE, 
                          Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, &P, &T, &x, &zone, &V0x, &E0x, &V0y, &E0y);
    */

    /*
    err=fPTC_cin_phase(0, cin_phase, TAUn1, E, &P, &C, &T, lambda_in, lambda_out, Ntau, dX, dt, Ttau, &zone,
                       nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis,
                       SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
    */
    err=fPTC_cin_phase(0, 0, sch_cin_phase, TAUn1, E, &P, &C, &T, lambda_test, lambda_test, Ntau, dX, dt, Ttau, &zone,
                       nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis,
                       SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
    if(err){return err;}
    
    /*
    err=fPTC_cin_phase(0, cin_phase, TAUn1, E+dE_inf, &Pp, &C, &T, lambda, Ntau, dX, dt, Ttau,
                       nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis,
                       SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
    if(err){return err;}

    err=fPTC_cin_phase(0, cin_phase, TAUn1, E-dE_inf, &Pm, &C, &T, lambda, Ntau, dX, dt, Ttau,
                       nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis,
                       SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
    if(err){return err;}
    */

    if(zone==-1){printf("Erreur out of domain : zone= %d (fPn1_etain_cin_phase)\n",zone); return 1;  }
    else if(zone==1 || zone==2 || zone==3){
      Pp=fP(TAUn1,E+dE_inf,zone);
      Pm=fP(TAUn1,E-dE_inf,zone);
    }
    else if(zone==12 || zone==13 || zone==23){
      j=(E-Emin_dis)/he;
      pos_fic=i*(Ne_dis-1) + j;
      V0x=tabVx[pos_fic];  E0x=tabEx[pos_fic];
      V0y=tabVy[pos_fic];  E0y=tabEy[pos_fic];

      phaseB=zone%10;
      phaseA=(zone-phaseB)/10;
      lambdaB=lambda_out[phaseB-1];

      err=VE_mixte_cin_phase(Nmax, epsilon, phaseA, phaseB, lambdaB, V0x, E0x, V0y, E0y, TAUn1, E+dE_inf, &Vx, &Ex, &Vy, &Ey);  if(err){return 1;}
      //err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V+dV, E+dE, &Vx, &Ex, &Vy, &Ey, &x);  if(err){return 1;}
      Pp=fP(Vx,Ex,phaseA);

      err=VE_mixte_cin_phase(Nmax, epsilon, phaseA, phaseB, lambdaB, V0x, E0x, V0y, E0y, TAUn1, E-dE_inf, &Vx, &Ex, &Vy, &Ey);  if(err){return 1;}
      //err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V-dV, E-dE, &Vx, &Ex, &Vy, &Ey, &x);  if(err){return 1;}
      Pm=fP(Vx,Ex,phaseA);
    }
    else if(zone==0){
      Pp=P;
      Pm=P;
    }

    dPpdE=(Pp-P)/dE_inf;
    dPmdE=(Pm-P)/(-dE_inf);
    dPdE=(dPpdE + dPmdE)/2;

    f= E-EPSn + ((P+Pn)/2 + Qn)*dTAU;
    f_prime= 1 + dPdE*dTAU/2;
    dE=(1.)*(-f/f_prime)/E;
   
    /*
    printf(" n= %d : f= %g , f_prime= %g , dPdE= %g , P= %g , Pp= %g , Pm= %g , TAUn1= %g , E= %g , zone= %d \n",n,f,f_prime,dPdE,P,Pp,Pm,TAUn1,E,zone ); 
    printf("         dE= %g, E-EPSn= %g, P+Pn= %g, Qn= %g, dTAU= %g \n",dE,E-EPSn,P+Pn,Qn,dTAU );
    printf("         Vx= %g , Ex= %g , phaseA= %d , V0x= %g , E0x= %g , V0y= %g , E0y= %g, nb_zone= %d , zone1= %d , zone2= %d , zone3= %d \n\n",Vx,Ex,phaseA,V0x,E0x,V0y,E0y,tab_Nb_ZONE[pos_fic],tabZONE[0][pos_fic],tabZONE[1][pos_fic],tabZONE[2][pos_fic] ); 
    */
    /*
    if(isnan(dE)){
      printf(" f= %g , f_prime= %g , dPdE= %g , P= %g , Pp= %g , Pm= %g , TAUn1= %g , E= %g , zone= %d , pos_fic= %d , i= %d , j= %d \n",f,f_prime,dPdE,P,Pp,Pm,TAUn1,E,zone,pos_fic,i,j ); 
      printf(" Vx= %g , Ex= %g , phaseA= %d , V0x= %g , E0x= %g , V0y= %g , E0y= %g, nb_zone= %d , zone1= %d , zone2= %d , zone3= %d \n",Vx,Ex,phaseA,V0x,E0x,V0y,E0y,tab_Nb_ZONE[pos_fic],tabZONE[0][pos_fic],tabZONE[1][pos_fic],tabZONE[2][pos_fic] ); 
      return 1;
    }
    */

    E+=dE;

    //critere=fabs(E-EPSn + ((P-Pn)/2 + Qn)*dTAU);
    critere=fabs(dE);
    
    /*
    printf("  P= %g , Pp= %g , Pm= %g\n",P,Pp,Pm );
    printf("  P-Pn= %g\n",(P-Pn)*1e-9);
    printf("  dPdE= %g\n",dPdE);
    printf("  n= %d , dE= %g, critere= %g , |f|= %g \n",n,dE,critere,fabs(f) );
    */

    n++;
  }
  //printf("  n= %d , dE= %g, critere= %g , |f|= %g , zone= %d \n",n,dE,critere,fabs(f),zone );
  if(n>=Nmax){printf("  n= %d , critere= %g , |f|= %g , zone= %d \n",n,critere,fabs(f),zone ); return 1;}

  *pEPSn1=E;
  *pPn1=P;

  free(lambda_out);

  return 0;
}



// Script pour tester le Benchmark
int use_benchmark_PT_cin_phase(int Nv_test, int Ne_test, double Vmin_test, double Emin_test, double Vmax_test, double Emax_test){
  int err;
  double E,V;
  double P,T;
  double x;
  double dt;

  int sch_cin_phase;
  
  // Pour le calcul de la celerité du son
  double V0x, E0x, V0y, E0y;

  // Discrétisation du benchmark/tabulation
  int nb_points;
  int Nv_dis,Ne_dis;
  double Vmin_dis, Vmax_dis;
  double Emin_dis, Emax_dis;

  double* lambda=malloc(3*sizeof(double));
  double Ttau;

  // Ecriture dans un fichier texte des paramètres de la discretisation
  /* pour pouvoir utiliser la discretisations juste en lisant les fichiers
  Organisation du fichier :
    nb_points
    Nv_dis , Ne_dis
    Vmin_dis , Vmax_dis
    Emin_dis , Emax_dis
  */
  FILE *fparams;
  if((fparams = fopen("fichiers/PARAMS_dis.txt", "r")) == NULL){printf("erreur ouverture fichier PARAMS_dis.txt \n"); return 1;}
  fscanf(fparams, "%d \n",&nb_points );
  fscanf(fparams, "%d %d\n",&Nv_dis,&Ne_dis );
  fscanf(fparams, "%lf %lf\n",&Vmin_dis,&Vmax_dis );
  fscanf(fparams, "%lf %lf\n",&Emin_dis,&Emax_dis );
  fclose(fparams);
  

  // Calcul de la célerité du son
  double Pp,Pm, Tp,Tm;
  double dE,c;
  double dV=1e-8;


  double hv_dis=(Vmax_dis-Vmin_dis)/(Nv_dis-1);
  double he_dis=(Emax_dis-Emin_dis)/(Ne_dis-1);

  double hv_test=(Vmax_test-Vmin_test)/(Nv_test-1);
  double he_test=(Emax_test-Emin_test)/(Ne_test-1);

  printf(" hv_test/hv_dis = %lf percent \n",100*hv_test/hv_dis );
  printf(" he_test/he_dis = %lf percent \n",100*he_test/he_dis );

  // Nombre de points sur les frontières
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3; 
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;

  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );


  double** SA=alloctabd(NA+1,2);
  double** SB=alloctabd(NB+1,2);
  double** SC=alloctabd(NC+1,2); // on ajoute S[N]=S[0]  --> utile pour l'algo
  double** SAB=alloctabd(NAB+1,2);
  double** SAC=alloctabd(NAC+1,2);
  double** SBC=alloctabd(NBC+1,2);

  // fichier binaire    couple (V,E) pour les zones mixte
  FILE *fileVmx, *fileEmx, *fileVmy, *fileEmy;
  if((fileVmx = fopen("fichiers/Vbench_mixte_x.bin", "rb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmx = fopen("fichiers/Ebench_mixte_x.bin", "rb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  if((fileVmy = fopen("fichiers/Vbench_mixte_y.bin", "rb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmy = fopen("fichiers/Ebench_mixte_y.bin", "rb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  
  double* tabVx=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double)); double* tabEx=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double));
  double* tabVy=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double)); double* tabEy=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double));

  // lecture des tableaux
  fread(tabVx,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileVmx);
  fread(tabEx,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileEmx);
  fread(tabVy,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileVmy);
  fread(tabEy,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileEmy);

  fclose(fileVmx);  fclose(fileEmx);
  fclose(fileVmy);  fclose(fileEmy);
  
  // Ouverture du fichier ZONE
  FILE *fileZONE_0,*fileZONE_1,*fileZONE_2,*fileZONE_3,*fileZONE_4;
  if((fileZONE_0 = fopen("fichiers/ZONEbench_0.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 0\n"); return 1;}
  if((fileZONE_1 = fopen("fichiers/ZONEbench_1.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 1\n"); return 1;}
  if((fileZONE_2 = fopen("fichiers/ZONEbench_2.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 2\n"); return 1;}
  if((fileZONE_3 = fopen("fichiers/ZONEbench_3.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 3\n"); return 1;}
  if((fileZONE_4 = fopen("fichiers/ZONEbench_4.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 4\n"); return 1;}
  
  int** tabZONE=alloctab(5,(Nv_dis-1)*(Ne_dis-1));
  
  // Ecriture des tableaux
  fread(tabZONE[0],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_0);
  fread(tabZONE[1],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_1);
  fread(tabZONE[2],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_2);
  fread(tabZONE[3],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_3);
  fread(tabZONE[4],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_4);

  fclose(fileZONE_0);  fclose(fileZONE_1);  fclose(fileZONE_2);
  fclose(fileZONE_3);  fclose(fileZONE_4);


  // fichier binaire NOMBRE DE ZONE
  FILE *file_Nb_ZONE;
  if((file_Nb_ZONE = fopen("fichiers/Nb_ZONEbench.bin", "rb")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  
  int* tab_Nb_ZONE=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(int)); 

  // Ecriture des tableaux
  fread(tab_Nb_ZONE,sizeof(int),(Nv_dis-1)*(Ne_dis-1),file_Nb_ZONE);

  fclose(file_Nb_ZONE);


  // fichier binaire POLYGONES Frontières de zones
  FILE *fileA_V,*fileB_V,*fileC_V,*fileAB_V,*fileAC_V,*fileBC_V;
  FILE *fileA_E,*fileB_E,*fileC_E,*fileAB_E,*fileAC_E,*fileBC_E;
  if((fileA_V = fopen("fichiers/A_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE A\n"); return 1;}     if((fileA_E = fopen("fichiers/A_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE A\n"); return 1;}
  if((fileB_V = fopen("fichiers/B_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE B\n"); return 1;}     if((fileB_E = fopen("fichiers/B_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE B\n"); return 1;}
  if((fileC_V = fopen("fichiers/C_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE C\n"); return 1;}     if((fileC_E = fopen("fichiers/C_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE C\n"); return 1;}
  if((fileAB_V = fopen("fichiers/AB_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AB\n"); return 1;}  if((fileAB_E = fopen("fichiers/AB_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AB\n"); return 1;}
  if((fileAC_V = fopen("fichiers/AC_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AC\n"); return 1;}  if((fileAC_E = fopen("fichiers/AC_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AC\n"); return 1;}
  if((fileBC_V = fopen("fichiers/BC_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE BC\n"); return 1;}  if((fileBC_E = fopen("fichiers/BC_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE BC\n"); return 1;}
  
  double* SA_V=malloc((NA+1)*sizeof(double));    double* SA_E=malloc((NA+1)*sizeof(double));
  double* SB_V=malloc((NB+1)*sizeof(double));    double* SB_E=malloc((NB+1)*sizeof(double));
  double* SC_V=malloc((NC+1)*sizeof(double));    double* SC_E=malloc((NC+1)*sizeof(double));
  double* SAB_V=malloc((NAB+1)*sizeof(double));  double* SAB_E=malloc((NAB+1)*sizeof(double));
  double* SAC_V=malloc((NAC+1)*sizeof(double));  double* SAC_E=malloc((NAC+1)*sizeof(double));
  double* SBC_V=malloc((NBC+1)*sizeof(double));  double* SBC_E=malloc((NBC+1)*sizeof(double));

  // Ecriture des tableaux
  fread(SA_V,sizeof(double),NA+1,fileA_V);     fread(SA_E,sizeof(double),NA+1,fileA_E);
  fread(SB_V,sizeof(double),NB+1,fileB_V);     fread(SB_E,sizeof(double),NB+1,fileB_E);
  fread(SC_V,sizeof(double),NC+1,fileC_V);     fread(SC_E,sizeof(double),NC+1,fileC_E);
  fread(SAB_V,sizeof(double),NAB+1,fileAB_V);  fread(SAB_E,sizeof(double),NAB+1,fileAB_E);
  fread(SAC_V,sizeof(double),NAC+1,fileAC_V);  fread(SAC_E,sizeof(double),NAC+1,fileAC_E);
  fread(SBC_V,sizeof(double),NBC+1,fileBC_V);  fread(SBC_E,sizeof(double),NBC+1,fileBC_E);

  for(int i=0; i<NA+1; i++){  SA[i][0]=SA_V[i];  SA[i][1]=SA_E[i];  }
  for(int i=0; i<NB+1; i++){  SB[i][0]=SB_V[i];  SB[i][1]=SB_E[i];  }
  for(int i=0; i<NC+1; i++){  SC[i][0]=SC_V[i];  SC[i][1]=SC_E[i];  }

  for(int i=0; i<NAB+1; i++){  SAB[i][0]=SAB_V[i];  SAB[i][1]=SAB_E[i];  }
  for(int i=0; i<NAC+1; i++){  SAC[i][0]=SAC_V[i];  SAC[i][1]=SAC_E[i];  }
  for(int i=0; i<NBC+1; i++){  SBC[i][0]=SBC_V[i];  SBC[i][1]=SBC_E[i];  }

  fclose(fileA_V);   fclose(fileB_V);   fclose(fileC_V);
  fclose(fileA_E);   fclose(fileB_E);   fclose(fileC_E);
  fclose(fileAB_V);  fclose(fileAC_V);  fclose(fileBC_V);
  fclose(fileAB_E);  fclose(fileAC_E);  fclose(fileBC_E);

  free(SA_V); free(SB_V); free(SC_V); free(SAB_V); free(SAC_V); free(SBC_V);
  free(SA_E); free(SB_E); free(SC_E); free(SAB_E); free(SAC_E); free(SBC_E);





 //////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////
  ///////////
  // Ecriture des sommets des poly pour vérif
  FILE *fSA, *fSB, *fSC, *fSAB, *fSAC, *fSBC;
  if((fSA = fopen("fichiers/pSA.txt", "w+")) == NULL){printf("erreur ouverture fichier pSA\n");return 1;}
  if((fSB = fopen("fichiers/pSB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSB\n");return 1;}
  if((fSC = fopen("fichiers/pSC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSC\n");return 1;}
  if((fSAB = fopen("fichiers/pSAB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAB\n");return 1;}
  if((fSAC = fopen("fichiers/pSAC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAC\n");return 1;}
  if((fSBC = fopen("fichiers/pSBC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSBC\n");return 1;}
  
  for(int i=0; i<NA+1; i++){   fprintf(fSA,"%.10lf %.10lf\n", SA[i][0], SA[i][1] );  }
  for(int i=0; i<NB+1; i++){   fprintf(fSB,"%.10lf %.10lf\n", SB[i][0], SB[i][1] );  }
  for(int i=0; i<NC+1; i++){   fprintf(fSC,"%.10lf %.10lf\n", SC[i][0], SC[i][1] );  }
  for(int i=0; i<NAB+1; i++){  fprintf(fSAB,"%.10lf %.10lf\n",SAB[i][0],SAB[i][1] );  }
  for(int i=0; i<NAC+1; i++){  fprintf(fSAC,"%.10lf %.10lf\n",SAC[i][0],SAC[i][1] );  }
  for(int i=0; i<NBC+1; i++){  fprintf(fSBC,"%.10lf %.10lf\n",SBC[i][0],SBC[i][1] );  }

  fclose(fSA);   fclose(fSB);   fclose(fSC);  fclose(fSAB);  fclose(fSAC);  fclose(fSBC);
 //////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////





  //  ****************************************
  // Pour tester les valeurs
  FILE *fPtest, *fTtest, *fXAtest, *fXBtest, *fXCtest, *fCtest;
  if((fPtest = fopen("fichiers/Ptest.txt", "w+")) == NULL){printf("erreur ouverture fichier P test\n");return 1;}
  if((fTtest = fopen("fichiers/Ttest.txt", "w+")) == NULL){printf("erreur ouverture fichier T test\n");return 1;}
  if((fXAtest = fopen("fichiers/XAtest.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  if((fXBtest = fopen("fichiers/XBtest.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  if((fXCtest = fopen("fichiers/XCtest.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  if((fCtest = fopen("fichiers/Ctest.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  FILE *fVtest, *fEtest;
  if((fVtest = fopen("fichiers/Vtest.txt", "w+")) == NULL){printf("erreur ouverture fichier V test \n");return 1;}
  if((fEtest = fopen("fichiers/Etest.txt", "w+")) == NULL){printf("erreur ouverture fichier E test \n");return 1;}
  for(int j=0; j<Ne_test; j++){
      E=Emin_test+j*he_test;
      fprintf(fEtest, "%.15lf ",E);
  }
  for(int j=0; j<Nv_test; j++){
      V=Vmin_test+j*hv_test;
      fprintf(fVtest, "%.15lf ",V);
  }
  fclose(fVtest);  fclose(fEtest);


  int zone;
  FILE *fZONEtest;
  if((fZONEtest = fopen("fichiers/ZONEtest.txt", "w+")) == NULL){printf("erreur ouverture fichier ZONE test\n");return 1;}
  //  ****************************************
  
  //  $$$$$$$   Ecriture texte pour verif   $$$$$$$$$$
  FILE *file_Nb_ZONE_test;
  if((file_Nb_ZONE_test = fopen("fichiers/Nb_ZONEbench.txt", "w+")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  for(int i=0; i<(Nv_dis-1)*(Ne_dis-1); i++){
    fprintf(file_Nb_ZONE_test, "%d \n",tab_Nb_ZONE[i] );
  }
  fclose(file_Nb_ZONE_test);

  FILE *file_V_dis;
  if((file_V_dis = fopen("fichiers/Vbench_dis.txt", "w+")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  for(int i=0; i<(Nv_dis-1); i++){
    fprintf(file_V_dis, "%.15lf \n",Vmin_dis+i*hv_dis+hv_dis/2 );
  }
  fclose(file_V_dis);
  
  FILE *file_E_dis;
  if((file_E_dis = fopen("fichiers/Ebench_dis.txt", "w+")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  for(int i=0; i<(Ne_dis-1); i++){
    fprintf(file_E_dis, "%.15lf \n",Emin_dis+i*he_dis+he_dis/2 );
  }
  fclose(file_E_dis);
  //  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  double Ntau=1, dx=1;


  // ***** Boucle *****
  for(int i=0; i<Nv_test; i++){
    printf("  - i= %d\n",i );
    V=Vmin_test+i*hv_test;
    for(int j=0; j<Ne_test; j++){
      //printf(" - i= %d, j= %d \n",i,j );
      //printf("   V= %g  E= %g\n",V,E );
      E=Emin_test+j*he_test;
      //err=lec_benchmark_PT(1,V, E, nb_points, SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE, 
      //                     Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, &P, &T, &x, &zone, &V0x, &E0x, &V0y, &E0y);
  
      err=fPTC_cin_phase(1, 0, sch_cin_phase, V, E, &P, &c, &T, lambda, lambda, Ntau, dx, dt, Ttau, &zone,
                         nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis,
                         SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
      if(err){return err;}
      
      //

      fprintf(fPtest,"%.10lf\n",P);
      fprintf(fTtest,"%.10lf\n",T);
      fprintf(fXAtest,"%.10lf\n",lambda[0]);
      fprintf(fXBtest,"%.10lf\n",lambda[1]);
      fprintf(fXCtest,"%.10lf\n",lambda[2]);
      fprintf(fZONEtest,"%d\n",zone);
      fprintf(fCtest,"%.10lf\n",c);      
    }
  }

  fclose(fPtest);  fclose(fTtest);  fclose(fXAtest); fclose(fXBtest); fclose(fXCtest); fclose(fCtest);
  

  // Liberation de la mémoire
  // Liberation
  free(tabVx);   free(tabEx);
  free(tabVy);   free(tabEy);
  
  freetab(tabZONE);  free(tab_Nb_ZONE);

  freetab(SA);  freetab(SB);  freetab(SC);
  freetab(SAB);  freetab(SAC);  freetab(SBC);
  
  return 0;
}



int debug_Newton(){
  printf("\n\n  DEBUG NEWTON  \n");
  int Nmax=1e2;
  double epsilon=1e-5;
  int err;
   
  double x;
  double V,E;
  V=0.0001; E=831111;

  double Vx, Ex, Vy, Ey;
  double V0x, E0x, V0y, E0y;
  V0x= 9.92845e-05; E0x= 806926; V0y= 0.000101004; E0y= 866875;
  //V0y*=(1.0+1e-3);
  //E0y*=(1.0+1e-3);
  printf(" V0x= %g , E0x= %g , V0y= %g , E0y= %g\n",V0x,E0x,V0y,E0y );

  err=VE_mixte(Nmax, epsilon, 2, 3, V0x, E0x, V0y, E0y, V, E, &Vx, &Ex, &Vy, &Ey, &x);
  if(err){return 1;}
  printf(" Vx= %g , Ex= %g , Vy= %g , Ey= %g\n",Vx,Ex,Vy,Ey );
  printf(" x= %g ,\n",x );
  
  double lambda3=x;


  //V0x= Vx; E0x= Ex; V0y= Vy; E0y= Ey;
  lambda3*=1.3;
  printf("lambda= %g\n",lambda3 );
  err= VE_mixte_cin_phase(Nmax, epsilon, 2, 3, lambda3, V0x, E0x, V0y, E0y, V, E, &Vx, &Ex, &Vy, &Ey);
  if(err){return 1;}

  printf("Vx= %g , Ex= %g , Vy= %g , Ey= %g\n",Vx,Ex,Vy,Ey );
  printf("P= %g  T= %g \n",fP(Vx,Ex,2), fT(Vx,Ex,3));

}




// IMPLEMENTATION DES METHODES TABULATIONS OU CALCUL PUR


// Fonction pour avoir une benchmark de (P,T) sur (V,E)
/*

*/
/*
int ecr_benchmark_pur_METH2(int nb_points, int Nv ,int Ne, double Vmin, double Vmax, double Emin, double Emax){

  // Ecriture dans un fichier texte des paramètres de la discretisation
  // pour pouvoir utiliser la discretisations juste en lisant les fichiers
  FILE *fparams;
  if((fparams = fopen("fichiers/PARAMS_dis_METH2.txt", "w+")) == NULL){printf("erreur ouverture fichier PARAMS_dis.txt \n"); return 1;}
  fprintf(fparams, "%d \n",nb_points );
  fprintf(fparams, "%d %d\n",Nv,Ne );
  fprintf(fparams, "%.15lf %.15lf\n",Vmin,Vmax );
  fprintf(fparams, "%.15lf %.15lf\n",Emin,Emax );
  fprintf(fparams, "\n\n Organisation du fichier\n");
  fprintf(fparams, "nb_points\n");
  fprintf(fparams, "Nv_dis  Ne_dis\n");
  fprintf(fparams, "Vmin_dis Vmax_dis\n");
  fprintf(fparams, "Emin_dis Emax_dis\n");


  double he=(Emax-Emin)/(Ne-1);
  double hv=(Vmax-Vmin)/(Nv-1);


  int Nmax=1e4;
  double epsilon=1e-6;
  double Vx, Ex, Vy, Ey, x;
  

  int err,res;
  

  double E,V;

  double Vmx, Emx;
  double Vmy, Emy;


  // Varaibles pour discretisation de la maille
  double V_sub, E_sub;
  int N_rec=10;
  double he_sub=he/(N_rec-1);
  double hv_sub=hv/(N_rec-1);
  double** VE_front_rec=alloctabd(4*(N_rec-1),2);
  int count;
  

  // nb_points est un argument de la fonction
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3; 
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;
  
  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );
  

  double** SA=alloctabd(NA+1,2);
  double** SB=alloctabd(NB+1,2);
  double** SC=alloctabd(NC+1,2); // on ajoute S[N]=S[0]  --> utile pour l'algo
  double** SAB=alloctabd(NAB+1,2);
  double** SAC=alloctabd(NAC+1,2);
  double** SBC=alloctabd(NBC+1,2);
  
   
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
  


  // CALCUL DES SOMMETS  
  printf(" Sommets_polygone \n");
  err=sommets_polygone(NA, NB, NC, NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3, N12, N13, N23, SA, SB, SC, SAB, SAC, SBC, segP12, segT12, segP13, segT13, segP23, segT23, 
                       segVa12, segEa12, segVb12, segEb12, segVa13, segEa13, segVc13, segEc13, segVb23, segEb23, segVc23, segEc23 );
  if(err){printf("err= %d (sommets_polygone)\n",err ); return err;}


  
  // Pour tester les valeurs (coord des points du maillage)
  FILE *fileV, *fileE;
  if((fileV = fopen("fichiers/Vbench_METH2.txt", "w+")) == NULL){printf("erreur ouverture fichier V bench \n"); return 1;}
  if((fileE = fopen("fichiers/Ebench_METH2.txt", "w+")) == NULL){printf("erreur ouverture fichier E bench \n"); return 1;}
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
  

  double* tabVx=malloc((Nv-1)*(Ne-1)*sizeof(double)); double* tabEx=malloc((Nv-1)*(Ne-1)*sizeof(double));
  double* tabVy=malloc((Nv-1)*(Ne-1)*sizeof(double)); double* tabEy=malloc((Nv-1)*(Ne-1)*sizeof(double));
  
  int* tabZONE=malloc((Nv-1)*(Ne-1)*sizeof(double));  // mémorise la zone dans le rectangle


  int test;

  int k=0;
  for(int i=0; i<Nv-1; i++){
    printf("  --i= %d\n",i );
    V=Vmin+i*hv + hv/2; // pour prendre la valeur au centre du rectangle
    for(int j=0; j<Ne-1; j++){
      //printf("  -- i= %d, j= %d\n",i,j );

      E=Emin+j*he + he/2;  
      

      // Attribution d'une zone aux points milieux des rectangles
      if(intPOLY(V,E,NC,SC)){ 
        tabZONE[k]=3;
      }
      else if(intPOLY(V,E,NB,SB)){
        tabZONE[k]=2;
      }      
      else if(intPOLY(V,E,NA,SA)){
        tabZONE[k]=1;
      }
      else if(ftri(V,E)){
        tabZONE[k]=0;
      }
      else if(intPOLY(V,E,NAB,SAB)){
        err=graine_zone_mixte(V, E, N12, segVa12, segVb12, segEa12, segEb12, &Vmx, &Emx, &Vmy, &Emy);
        if(err){printf("err= %d (benchmark_PT)\n",err ); return err;}
        tabVx[k]=Vmx; tabEx[k]=Emx;
        tabVy[k]=Vmy; tabEy[k]=Emy;
        tabZONE[k]=12;
      }
      else if(intPOLY(V,E,NAC,SAC)){
        err=graine_zone_mixte(V, E, N13, segVa13, segVc13, segEa13, segEc13, &Vmx, &Emx, &Vmy, &Emy);
        if(err){printf("err= %d (benchmark_PT)\n",err ); return err;}
        tabVx[k]=Vmx; tabEx[k]=Emx;
        tabVy[k]=Vmy; tabEy[k]=Emy;
        tabZONE[k]=13;
      }
      else if(intPOLY(V,E,NBC,SBC)){
        err=graine_zone_mixte(V, E, N23, segVb23, segVc23, segEb23, segEc23, &Vmx, &Emx, &Vmy, &Emy);
        if(err){printf("err= %d (benchmark_PT)\n",err ); return err;}
        tabVx[k]=Vmx; tabEx[k]=Emx;
        tabVy[k]=Vmy; tabEy[k]=Emy;
        tabZONE[k]=23;
      }
      else{  // hors domaine --> on ne se soucis des pb de chevauchement des rectangles sur les zones
        tabZONE[k]=-1;
      }
    }
  }


  // fichier binaire   couple (V,E) pour les zones mixte
  FILE *fileVmx, *fileEmx;
  if((fileVmx = fopen("fichiers/Vbench_mixte_x_METH2.bin", "wb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmx = fopen("fichiers/Ebench_mixte_x_METH2.bin", "wb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  FILE *fileVmy, *fileEmy;
  if((fileVmy = fopen("fichiers/Vbench_mixte_y_METH2.bin", "wb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmy = fopen("fichiers/Ebench_mixte_y_METH2.bin", "wb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  
  // Ecriture des tableaux
  fwrite(tabVx,sizeof(double),(Nv-1)*(Ne-1),fileVmx);
  fwrite(tabEx,sizeof(double),(Nv-1)*(Ne-1),fileEmx);
  fwrite(tabVy,sizeof(double),(Nv-1)*(Ne-1),fileVmy);
  fwrite(tabEy,sizeof(double),(Nv-1)*(Ne-1),fileEmy);

  fclose(fileVmx);  fclose(fileEmx);
  fclose(fileVmy);  fclose(fileEmy);
  
  // fichier binaire ZONE
  FILE *fileZONE;
  if((fileZONE = fopen("fichiers/ZONEbench_METH2.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE bench \n"); return 1;}
  
  // Ecriture des tableaux
  fwrite(tabZONE,sizeof(int),(Nv-1)*(Ne-1),fileZONE);

  fclose(fileZONE_0);


  // fichier binaire POLYGONES Frontières de zones
  FILE *fileA_V,*fileB_V,*fileC_V,*fileAB_V,*fileAC_V,*fileBC_V;
  FILE *fileA_E,*fileB_E,*fileC_E,*fileAB_E,*fileAC_E,*fileBC_E;
  if((fileA_V = fopen("fichiers/A_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE A (ecr_benchmark_PT)\n"); return 1;}     if((fileA_E = fopen("fichiers/A_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE A (ecr_benchmark_PT)\n"); return 1;}
  if((fileB_V = fopen("fichiers/B_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE B (ecr_benchmark_PT)\n"); return 1;}     if((fileB_E = fopen("fichiers/B_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE B (ecr_benchmark_PT)\n"); return 1;}
  if((fileC_V = fopen("fichiers/C_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE C (ecr_benchmark_PT)\n"); return 1;}     if((fileC_E = fopen("fichiers/C_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE C (ecr_benchmark_PT)\n"); return 1;}
  if((fileAB_V = fopen("fichiers/AB_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE AB (ecr_benchmark_PT)\n"); return 1;}  if((fileAB_E = fopen("fichiers/AB_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE AB (ecr_benchmark_PT)\n"); return 1;}
  if((fileAC_V = fopen("fichiers/AC_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE AC (ecr_benchmark_PT)\n"); return 1;}  if((fileAC_E = fopen("fichiers/AC_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE AC (ecr_benchmark_PT)\n"); return 1;}
  if((fileBC_V = fopen("fichiers/BC_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE BC (ecr_benchmark_PT)\n"); return 1;}  if((fileBC_E = fopen("fichiers/BC_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE BC (ecr_benchmark_PT)\n"); return 1;}
  
  double* SA_V=malloc((NA+1)*sizeof(double));    double* SA_E=malloc((NA+1)*sizeof(double));
  double* SB_V=malloc((NB+1)*sizeof(double));    double* SB_E=malloc((NB+1)*sizeof(double));
  double* SC_V=malloc((NC+1)*sizeof(double));    double* SC_E=malloc((NC+1)*sizeof(double));
  double* SAB_V=malloc((NAB+1)*sizeof(double));  double* SAB_E=malloc((NAB+1)*sizeof(double));
  double* SAC_V=malloc((NAC+1)*sizeof(double));  double* SAC_E=malloc((NAC+1)*sizeof(double));
  double* SBC_V=malloc((NBC+1)*sizeof(double));  double* SBC_E=malloc((NBC+1)*sizeof(double));

  for(int i=0; i<NA+1; i++){  SA_V[i]=SA[i][0];  SA_E[i]=SA[i][1];  }
  for(int i=0; i<NB+1; i++){  SB_V[i]=SB[i][0];  SB_E[i]=SB[i][1];  }
  for(int i=0; i<NC+1; i++){  SC_V[i]=SC[i][0];  SC_E[i]=SC[i][1];  }

  for(int i=0; i<NAB+1; i++){  SAB_V[i]=SAB[i][0];  SAB_E[i]=SAB[i][1];  }
  for(int i=0; i<NAC+1; i++){  SAC_V[i]=SAC[i][0];  SAC_E[i]=SAC[i][1];  }
  for(int i=0; i<NBC+1; i++){  SBC_V[i]=SBC[i][0];  SBC_E[i]=SBC[i][1];  }

  // Ecriture des tableaux
  fwrite(SA_V,sizeof(double),NA+1,fileA_V);     fwrite(SA_E,sizeof(double),NA+1,fileA_E);
  fwrite(SB_V,sizeof(double),NB+1,fileB_V);     fwrite(SB_E,sizeof(double),NB+1,fileB_E);
  fwrite(SC_V,sizeof(double),NC+1,fileC_V);     fwrite(SC_E,sizeof(double),NC+1,fileC_E);
  fwrite(SAB_V,sizeof(double),NAB+1,fileAB_V);  fwrite(SAB_E,sizeof(double),NAB+1,fileAB_E);
  fwrite(SAC_V,sizeof(double),NAC+1,fileAC_V);  fwrite(SAC_E,sizeof(double),NAC+1,fileAC_E);
  fwrite(SBC_V,sizeof(double),NBC+1,fileBC_V);  fwrite(SBC_E,sizeof(double),NBC+1,fileBC_E);


  fclose(fileA_V);   fclose(fileB_V);   fclose(fileC_V);
  fclose(fileA_E);   fclose(fileB_E);   fclose(fileC_E);
  fclose(fileAB_V);  fclose(fileAC_V);  fclose(fileBC_V);
  fclose(fileAB_E);  fclose(fileAC_E);  fclose(fileBC_E);

  free(SA_V); free(SB_V); free(SC_V); free(SAB_V); free(SAC_V); free(SBC_V);
  free(SA_E); free(SB_E); free(SC_E); free(SAB_E); free(SAC_E); free(SBC_E);



  // Liberation
  free(tabVx);   free(tabEx);
  free(tabVy);   free(tabEy);
  
  freetab(tabZONE);
  freetab(REC);

  free(tabNbZONE);
  freetab(tabINTER);
  freetab(tabInd_Mixte);
  freetab(numZONE);
  freetab(VE_front_rec);


  free(segP12);   free(segT12);
  free(segVa12);  free(segEa12);
  free(segVb12);  free(segEb12);
  
  free(segP13);   free(segT13);
  free(segVa13);  free(segEa13);
  free(segVc13);  free(segEc13);
  
  free(segP23);   free(segT23);
  free(segVb23);  free(segEb23);
  free(segVc23);  free(segEc23);

  freetab(SA);   freetab(SB);   freetab(SC);
  freetab(SAB);  freetab(SAC);  freetab(SBC);
  freetab(Sa12); freetab(Sb12);
  freetab(Sa13); freetab(Sc13);
  freetab(Sb23); freetab(Sc23);

  freetab(STRI_12);  freetab(STRI_13);  freetab(STRI_23);

  fprintf(fparams, "\nEcriture bien terminé\n");
  fclose(fparams);


  return 0;

}
*/

/*  Arguments d'entrée :
     - V,E
     - nb_points : nombre de points de discretisation des frontières
     - SA,SB,SC etc.. contiennent les discretisations des frontières
    Arguments de sorties :
     - Xfm,pzone : pointeurs sur les fraction massique et la zone de (V,E)
     - pV0x, pE0x, etc...  pointeur sur la valeur initiale choisi pour lancer une zone mixte (On peut y avoir acces à la fonction avant grace à tabVx)
*/

int cal_zone(int affichage, double V, double E, int nb_points, double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC, 
                double* Xfm, int* pzone, double* pV0x, double* pE0x, double* pV0y, double* pE0y){
 
  int Nmax=1e2;
  double epsilon=1e-8;

  int err; 
  double P,T;
  double Vmx, Emx, Vmy, Emy;
  double Vx, Ex, Vy, Ey;

  double x;


  // pour le calcul des fractions massiques dans la zone trinagle --> calcul des coord barycentriques dans le triangle
  double VaT=129.905671875844291e-6;  double EaT=72.168190127265447e3;
  double VbT=127.666112185063994e-6;  double EbT=99.480256119236515e3;
  double VcT=132.538478805442338e-6;  double EcT=136.154346228414312e3;
  double* VE=malloc(2*sizeof(double));
  VE[0]=V; VE[1]=E;
  double** TRI=alloctabd(3,2);
  TRI[0][0]=VaT; TRI[0][1]=EaT;
  TRI[1][0]=VbT; TRI[1][1]=EbT;
  TRI[2][0]=VcT; TRI[2][1]=EcT;


  // nb_points est un argument de la fonction
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3; 
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;

  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );

  
  for(int i=0; i<3; i++){ Xfm[i]=0; }

  //printf(" nb_zone= %d\n",nb_zone );
  if(intPOLY(V,E,NC,SC)){   Xfm[2]=1.0;  *pzone=3;  return 0;  }
  else if(intPOLY(V,E,NB,SB)){ Xfm[1]=1.0;  *pzone=2;  return 0;  }
  else if(intPOLY(V,E,NA,SA)){ Xfm[0]=1.0;  *pzone=1;  return 0;  }
  else if(intPOLY(V,E,NBC,SBC)){
    err= graine_zone_mixte_modif(V, E, N23, SBC, &Vmx, &Emx, &Vmy, &Emy);
    if(err){printf("err= %d BC (graine_zone_mixte_modif)\n",err ); return err;}
    err=VE_mixte(Nmax, epsilon, 2, 3, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x);
    if(err){printf("err= %d (VE_mixte)\n",err ); return 1;}
    Xfm[2]=x;  Xfm[1]=1.0-x;  *pzone=23;

    return 0;
  }
  else if(intPOLY(V,E,NAC,SAC)){
    //err=graine_zone_mixte(V, E, N13, segVa13, segVc13, segEa13, segEc13, &Vmx, &Emx, &Vmy, &Emy);
    err= graine_zone_mixte_modif(V, E, N13, SAC, &Vmx, &Emx, &Vmy, &Emy);
    if(err){printf("err= %d AC (graine_zone_mixte_modif)\n",err ); return err;}
    err=VE_mixte(Nmax, epsilon, 1, 3, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x);
    if(err){printf("err= %d (VE_mixte)\n",err ); return 1;}
    Xfm[2]=x;  Xfm[0]=1.0-x;  *pzone=13;

    return 0;
  }
  else if(intPOLY(V,E,NAB,SAB)){
    //err=graine_zone_mixte(V, E, N12, segVa12, segVb12, segEa12, segEb12, &Vmx, &Emx, &Vmy, &Emy);
    err= graine_zone_mixte_modif(V, E, N12, SAB, &Vmx, &Emx, &Vmy, &Emy);
    if(err){printf("err= %d AB (graine_zone_mixte_modif)\n",err ); return err;}
    err=VE_mixte(Nmax, epsilon, 1, 2, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x);
    if(err){printf("err= %d (VE_mixte)\n",err ); return 1;}
    Xfm[1]=x;  Xfm[0]=1.0-x;  *pzone=12;

    return 0;
  }
  else if(ftri(V,E)){
    err=coorbary(VE, TRI, &Xfm); if(err){return 1;} 
    *pzone=-1;  return 0;
  }
  else{
   if(affichage==1){Xfm[0]=-1.0;  Xfm[1]=-1.0;  Xfm[2]=-1.0;  *pzone=-1;  return 0;}
   else if(affichage==0){printf("Error out of domain : V= %g , E= %g\n",V,E); return 1;}
   else{printf("Erreur dans la valeur de affichage= %d",affichage); return 1;} 
  }

  return 1;
  
}



// Renvoie la valeur de P,T et C pour un couple V,E en prenant en compte la cinétique de chgmt de phase
/*  Argments d'entrees :
      - CIN_PAHSE :  if CIN_PHASE==0 : pas de cinétique de phase
                     if CIN_PHASE==1 : prise en compte de la cinétique de phase

   Metre un param entier tels que 
    meth==0 : calcul + tabulation
    meth==1 : calcul pure
    meth==2 : tabulation pure
    meth==3 : maillage ??

    Pour le moment uniquement avec cal_tab
*/
/*
int fPTC_METH_cin_phase(int methode, int affichage, int CIN_PHASE, double V, double E, double* pP, double *pc, double* pT, double* lambda_in, double* lambda_out, int Ntau, double dx, double dt, double Ttau, int* pzone,
                   int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, 
                   double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC, 
                   double* tabVx, double* tabEx, double* tabVy, double* tabEy, int** tabZONE, int* tab_Nb_ZONE ){
  
  //      Ttau est calculé avec c(Vn, En)   
  
  if(CIN_PHASE==2){  Ttau=Ntau*dx/(*pc);   }
  //double Ttau=1e-9;   //Ntau*dx/(*pc);

  int Nmax=3e2;
  double epsilon=2e-9;
  int err;
  // Pour le calcul de la celerité du son
  double V0x, E0x, V0y, E0y;  
  double Vx, Ex, Vy, Ey;  
  double dV=1e-8;

  int zone;

  double* Xfm=malloc(3*sizeof(double)); // tableau des fractions massiques
  int phaseX,phaseY;

  // Point triple
  double PT= 4510421781.609182357788085;
  double TT= 624.63982538492518;

  //err= PT_tab_cal(V, E, nb_points, SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE, 
  //                Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, pP, pT, px, &zone, &V0x, &E0x, &V0y, &E0y);

  // Focntion renvoie la zone, la fraction massique et les graines pour les Newton_mixte
  
  if(methode==0){ // methode mixte
    err= lec_benchmark_zone(affichage, V, E, nb_points, SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE, 
                            Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, Xfm, &zone, &V0x, &E0x, &V0y, &E0y);
  }
  else if(methode==1){ // calcul pur
    err= cal_zone(affichage, V, E, nb_points, SA, SB, SC, SAB, SAC, SBC, Xfm, &zone, &V0x, &E0x, &V0y, &E0y);
  }
  else if(methode==1){ // tabulation pure
  if(err){return err;}
  
  *pzone=zone;


  //printf(" ZONE= %d\n",zone );
  
  //for(int i=0; i<3; i++){ if(Xfm[i]<0 ){Xfm[i]=0;}  if(Xfm[i]>1 ){Xfm[i]=1;} }
  
  //printf("  V= %g , E= %g , zone= %d\n",V,E,zone );
  //printf(" V= %g, E= %g, Xfm[0]= %g, Xfm[2]= %g, Xfm[2]= %g \n",V,E,Xfm[0],Xfm[1],Xfm[2] );  

  if(zone==-1){
    if(affichage==1){*pP=-10e9; *pT=0.;  lambda_out[0]=-1; lambda_out[1]=-1; lambda_out[2]=-1;  return 0;}
    else if(affichage==0){printf("Error out of domain : V= %g , E= %g\n",V,E); return 1;}
    else{printf("Erreur affichage= %d",affichage); return 1;} 
  }
  if(CIN_PHASE==0){
    for(int i=0; i<3; i++){  lambda_out[i]=Xfm[i];  }
    // if mehode==2 (tabulation pure) return 0   car on a deja toutes les valeurs à l'équilibre 
    // en revanche si on active la cinétique de chgmt de phase, on doit recalculer P,T et C si on est en zone de cgmt de phase 
  }
  else if(CIN_PHASE==1 || CIN_PHASE==2){
    // Calcul de x^{n+1} 
    if(zone==1 || zone==2 || zone==3){ for(int i=0; i<3; i++){  lambda_out[i]=Xfm[i]; } }
    else if(zone==0){
      //printf("DlambdaA= %g , DlambdaB= %g \n", (Ttau*(lambda[0])+Xfm[0])/(Ttau+1),(Ttau*(lambda[1])+Xfm[1])/(Ttau+1));
      lambda_out[0]=((dt/Ttau)*Xfm[0]+lambda_in[0])/((dt/Ttau)+1.);
      lambda_out[1]=((dt/Ttau)*Xfm[1]+lambda_in[1])/((dt/Ttau)+1.);
      lambda_out[2]=1.-lambda_out[1]-lambda_out[0];

    }
    else if(zone==12 ){
      //printf("DlambdaA= %g \n", (Ttau*(lambda[0])+Xfm[0])/(Ttau+1));
      lambda_out[0]=((dt/Ttau)*Xfm[0]+lambda_in[0])/((dt/Ttau)+1.);
      lambda_out[1]=1.-lambda_out[0]; // (autre phase)
      lambda_out[2]=0.;
      //printf(" ok 12  lambda[0]= %g , Xfm[0]= %g\n",lambda[0],Xfm[0]);
    }
    else if( zone==13 ){
      //printf("DlambdaA= %g \n", (Ttau*(lambda[0])+Xfm[0])/(Ttau+1));
      lambda_out[0]=((dt/Ttau)*Xfm[0]+lambda_in[0])/((dt/Ttau)+1.);
      lambda_out[1]=0;
      lambda_out[2]=1.-lambda_out[0]; // (autre phase)
      
    }
    else if( zone==23 ){
      //printf("DlambdaB= %g \n", (Ttau*(lambda[0])+Xfm[0])/(Ttau+1));
      lambda_out[0]=0;
      lambda_out[1]=((dt/Ttau)*Xfm[1]+lambda_in[1])/((dt/Ttau)+1.);
      lambda_out[2]=1.-lambda_out[1]; // (autre phase)
      
    }
  }
  else{printf("Erreur CIN_PHASE= %d\n",CIN_PHASE ); return 1;}


  //printf(" ZONE= %d | In : xA= %g , xB= %g, xC= %g   Out : xA= %g , xB= %g, xC= %g\n",zone,lambda_in[0],lambda_in[1],lambda_in[2],lambda_out[0],lambda_out[1],lambda_out[2] );  
  
 

  /*

   Si methode==2 + cinétique, on ne recalcul P,T et C que si on est en zone bi-phasique sinon on resort ici
   ( apres la modif des fractions massiques )

   if(methode==2){return 0;}

   Si la cinétique est activé alors on lance aussi le Newton et le calcul de C

  */



/*
  // HYDRO
  //printf("   debut hydro  :  zone= %d\n",zone);
  if(zone==0){  
    *pP=PT;
    *pT=TT;
    *pc=0;

    return 0;
  }
  else if(zone==12 || zone==13 || zone==23){
    phaseY=zone%10;
    phaseX=(zone-phaseY)/10;
    // Calcul de p,T et C grace aux nouvelles valeurs de fractions massiques
    if(CIN_PHASE!=0){
      err= VE_mixte_cin_phase(Nmax, epsilon, phaseX, phaseY, lambda_out[phaseY-1], V0x, E0x, V0y, E0y, V, E, &Vx, &Ex, &Vy, &Ey);
      if(err){return 1;}
      *pP=fP(Vx,Ex,phaseX);
      *pT=fT(Vx,Ex,phaseX);
    }
    else if(CIN_PHASE==0){
      *pP=fP(V0x,E0x,phaseX);
      *pT=fT(V0x,E0x,phaseX);
    } 
  }
  else if(zone==1 || zone==2 || zone==3){
    *pP=fP(V,E,zone);
    *pT=fT(V,E,zone);
  }

  /*
  if(*pP>1e11 || *pT<-1e3){
    printf(" Vx= %g , Ex= %g , Vy= %g , Ey= %g\n",Vx,Ex,Vy,Ey);
    printf(" V= %g , E= %g , zone= %d, xA= %g , xB= %g , xC= %g\n",V,E,zone,lambda[0],lambda[01],lambda[2]);
    printf(" Px= %g , Tx= %g , Py= %g , Ty= %g \n",fP(Vx,Ex,phaseX),fT(Vx,Ex,phaseX),fP(Vy,Ey,phaseY),fT(Vy,Ey,phaseY));
    return 1;
  }
  */
  
  /*
  // Problème de cérité du son 
  // On ne peut calculer P de la même manière, il faut utiliser VE_mixte_cin_phase si on est en zone de mélange mixte
  // Faire une nouvelle fonction celerite avec un test dedans
  err=c_celerite_son_cin_phase(Nmax, epsilon, zone, *pP, V, E, dV, V0x, E0x, V0y, E0y, lambda_out, pc);  
  if(err){return 1;}

  return 0;
}
*/

/*
// Script pour tester le Benchmark
int use_fPTC_METH_cin_phase(int methode, int Nv_test, int Ne_test, double Vmin_test, double Emin_test, double Vmax_test, double Emax_test){
  int err;
  double E,V;
  double P,T;
  double x;
  double dt;
  
  // Pour le calcul de la celerité du son
  double V0x, E0x, V0y, E0y;

  // Discrétisation du benchmark/tabulation
  int nb_points;
  int Nv_dis,Ne_dis;
  double Vmin_dis, Vmax_dis;
  double Emin_dis, Emax_dis;

  double* lambda=malloc(3*sizeof(double));
  double Ttau;

  // Ecriture dans un fichier texte des paramètres de la discretisation
  // pour pouvoir utiliser la discretisations juste en lisant les fichiers
  //Organisation du fichier :
  //  nb_points
  //  Nv_dis , Ne_dis
  //  Vmin_dis , Vmax_dis
  //  Emin_dis , Emax_dis
  
  FILE *fparams;
  if((fparams = fopen("fichiers/PARAMS_dis.txt", "r")) == NULL){printf("erreur ouverture fichier PARAMS_dis.txt \n"); return 1;}
  fscanf(fparams, "%d \n",&nb_points );
  fscanf(fparams, "%d %d\n",&Nv_dis,&Ne_dis );
  fscanf(fparams, "%lf %lf\n",&Vmin_dis,&Vmax_dis );
  fscanf(fparams, "%lf %lf\n",&Emin_dis,&Emax_dis );
  fclose(fparams);
  

  // Calcul de la célerité du son
  double Pp,Pm, Tp,Tm;
  double dE,c;
  double dV=1e-8;


  double hv_dis=(Vmax_dis-Vmin_dis)/(Nv_dis-1);
  double he_dis=(Emax_dis-Emin_dis)/(Ne_dis-1);

  double hv_test=(Vmax_test-Vmin_test)/(Nv_test-1);
  double he_test=(Emax_test-Emin_test)/(Ne_test-1);

  printf(" hv_test/hv_dis = %lf percent \n",100*hv_test/hv_dis );
  printf(" he_test/he_dis = %lf percent \n",100*he_test/he_dis );

  // Nombre de points sur les frontières
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3; 
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;

  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );

  double** SA=alloctabd(NA+1,2);
  double** SB=alloctabd(NB+1,2);
  double** SC=alloctabd(NC+1,2); // on ajoute S[N]=S[0]  --> utile pour l'algo
  double** SAB=alloctabd(NAB+1,2);
  double** SAC=alloctabd(NAC+1,2);
  double** SBC=alloctabd(NBC+1,2);

  // fichier binaire    couple (V,E) pour les zones mixte
  FILE *fileVmx, *fileEmx, *fileVmy, *fileEmy;
  if((fileVmx = fopen("fichiers/Vbench_mixte_x.bin", "rb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmx = fopen("fichiers/Ebench_mixte_x.bin", "rb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  if((fileVmy = fopen("fichiers/Vbench_mixte_y.bin", "rb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmy = fopen("fichiers/Ebench_mixte_y.bin", "rb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  
  double* tabVx=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double)); double* tabEx=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double));
  double* tabVy=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double)); double* tabEy=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double));

  // lecture des tableaux
  fread(tabVx,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileVmx);
  fread(tabEx,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileEmx);
  fread(tabVy,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileVmy);
  fread(tabEy,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileEmy);

  fclose(fileVmx);  fclose(fileEmx);
  fclose(fileVmy);  fclose(fileEmy);
  
  // Ouverture du fichier ZONE
  FILE *fileZONE_0,*fileZONE_1,*fileZONE_2,*fileZONE_3,*fileZONE_4;
  if((fileZONE_0 = fopen("fichiers/ZONEbench_0.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 0\n"); return 1;}
  if((fileZONE_1 = fopen("fichiers/ZONEbench_1.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 1\n"); return 1;}
  if((fileZONE_2 = fopen("fichiers/ZONEbench_2.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 2\n"); return 1;}
  if((fileZONE_3 = fopen("fichiers/ZONEbench_3.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 3\n"); return 1;}
  if((fileZONE_4 = fopen("fichiers/ZONEbench_4.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 4\n"); return 1;}
  
  int** tabZONE=alloctab(5,(Nv_dis-1)*(Ne_dis-1));
  
  // Ecriture des tableaux
  fread(tabZONE[0],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_0);
  fread(tabZONE[1],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_1);
  fread(tabZONE[2],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_2);
  fread(tabZONE[3],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_3);
  fread(tabZONE[4],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_4);

  fclose(fileZONE_0);  fclose(fileZONE_1);  fclose(fileZONE_2);
  fclose(fileZONE_3);  fclose(fileZONE_4);


  // fichier binaire NOMBRE DE ZONE
  FILE *file_Nb_ZONE;
  if((file_Nb_ZONE = fopen("fichiers/Nb_ZONEbench.bin", "rb")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  
  int* tab_Nb_ZONE=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(int)); 

  // Ecriture des tableaux
  fread(tab_Nb_ZONE,sizeof(int),(Nv_dis-1)*(Ne_dis-1),file_Nb_ZONE);

  fclose(file_Nb_ZONE);


  // fichier binaire POLYGONES Frontières de zones
  FILE *fileA_V,*fileB_V,*fileC_V,*fileAB_V,*fileAC_V,*fileBC_V;
  FILE *fileA_E,*fileB_E,*fileC_E,*fileAB_E,*fileAC_E,*fileBC_E;
  if((fileA_V = fopen("fichiers/A_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE A\n"); return 1;}     if((fileA_E = fopen("fichiers/A_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE A\n"); return 1;}
  if((fileB_V = fopen("fichiers/B_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE B\n"); return 1;}     if((fileB_E = fopen("fichiers/B_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE B\n"); return 1;}
  if((fileC_V = fopen("fichiers/C_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE C\n"); return 1;}     if((fileC_E = fopen("fichiers/C_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE C\n"); return 1;}
  if((fileAB_V = fopen("fichiers/AB_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AB\n"); return 1;}  if((fileAB_E = fopen("fichiers/AB_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AB\n"); return 1;}
  if((fileAC_V = fopen("fichiers/AC_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AC\n"); return 1;}  if((fileAC_E = fopen("fichiers/AC_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AC\n"); return 1;}
  if((fileBC_V = fopen("fichiers/BC_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE BC\n"); return 1;}  if((fileBC_E = fopen("fichiers/BC_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE BC\n"); return 1;}
  
  double* SA_V=malloc((NA+1)*sizeof(double));    double* SA_E=malloc((NA+1)*sizeof(double));
  double* SB_V=malloc((NB+1)*sizeof(double));    double* SB_E=malloc((NB+1)*sizeof(double));
  double* SC_V=malloc((NC+1)*sizeof(double));    double* SC_E=malloc((NC+1)*sizeof(double));
  double* SAB_V=malloc((NAB+1)*sizeof(double));  double* SAB_E=malloc((NAB+1)*sizeof(double));
  double* SAC_V=malloc((NAC+1)*sizeof(double));  double* SAC_E=malloc((NAC+1)*sizeof(double));
  double* SBC_V=malloc((NBC+1)*sizeof(double));  double* SBC_E=malloc((NBC+1)*sizeof(double));

  // Ecriture des tableaux
  fread(SA_V,sizeof(double),NA+1,fileA_V);     fread(SA_E,sizeof(double),NA+1,fileA_E);
  fread(SB_V,sizeof(double),NB+1,fileB_V);     fread(SB_E,sizeof(double),NB+1,fileB_E);
  fread(SC_V,sizeof(double),NC+1,fileC_V);     fread(SC_E,sizeof(double),NC+1,fileC_E);
  fread(SAB_V,sizeof(double),NAB+1,fileAB_V);  fread(SAB_E,sizeof(double),NAB+1,fileAB_E);
  fread(SAC_V,sizeof(double),NAC+1,fileAC_V);  fread(SAC_E,sizeof(double),NAC+1,fileAC_E);
  fread(SBC_V,sizeof(double),NBC+1,fileBC_V);  fread(SBC_E,sizeof(double),NBC+1,fileBC_E);

  for(int i=0; i<NA+1; i++){  SA[i][0]=SA_V[i];  SA[i][1]=SA_E[i];  }
  for(int i=0; i<NB+1; i++){  SB[i][0]=SB_V[i];  SB[i][1]=SB_E[i];  }
  for(int i=0; i<NC+1; i++){  SC[i][0]=SC_V[i];  SC[i][1]=SC_E[i];  }

  for(int i=0; i<NAB+1; i++){  SAB[i][0]=SAB_V[i];  SAB[i][1]=SAB_E[i];  }
  for(int i=0; i<NAC+1; i++){  SAC[i][0]=SAC_V[i];  SAC[i][1]=SAC_E[i];  }
  for(int i=0; i<NBC+1; i++){  SBC[i][0]=SBC_V[i];  SBC[i][1]=SBC_E[i];  }

  fclose(fileA_V);   fclose(fileB_V);   fclose(fileC_V);
  fclose(fileA_E);   fclose(fileB_E);   fclose(fileC_E);
  fclose(fileAB_V);  fclose(fileAC_V);  fclose(fileBC_V);
  fclose(fileAB_E);  fclose(fileAC_E);  fclose(fileBC_E);

  free(SA_V); free(SB_V); free(SC_V); free(SAB_V); free(SAC_V); free(SBC_V);
  free(SA_E); free(SB_E); free(SC_E); free(SAB_E); free(SAC_E); free(SBC_E);





 //////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////
  ///////////
  // Ecriture des sommets des poly pour vérif
  FILE *fSA, *fSB, *fSC, *fSAB, *fSAC, *fSBC;
  if((fSA = fopen("fichiers/pSA.txt", "w+")) == NULL){printf("erreur ouverture fichier pSA\n");return 1;}
  if((fSB = fopen("fichiers/pSB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSB\n");return 1;}
  if((fSC = fopen("fichiers/pSC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSC\n");return 1;}
  if((fSAB = fopen("fichiers/pSAB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAB\n");return 1;}
  if((fSAC = fopen("fichiers/pSAC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAC\n");return 1;}
  if((fSBC = fopen("fichiers/pSBC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSBC\n");return 1;}
  
  for(int i=0; i<NA+1; i++){   fprintf(fSA,"%.10lf %.10lf\n", SA[i][0], SA[i][1] );  }
  for(int i=0; i<NB+1; i++){   fprintf(fSB,"%.10lf %.10lf\n", SB[i][0], SB[i][1] );  }
  for(int i=0; i<NC+1; i++){   fprintf(fSC,"%.10lf %.10lf\n", SC[i][0], SC[i][1] );  }
  for(int i=0; i<NAB+1; i++){  fprintf(fSAB,"%.10lf %.10lf\n",SAB[i][0],SAB[i][1] );  }
  for(int i=0; i<NAC+1; i++){  fprintf(fSAC,"%.10lf %.10lf\n",SAC[i][0],SAC[i][1] );  }
  for(int i=0; i<NBC+1; i++){  fprintf(fSBC,"%.10lf %.10lf\n",SBC[i][0],SBC[i][1] );  }

  fclose(fSA);   fclose(fSB);   fclose(fSC);  fclose(fSAB);  fclose(fSAC);  fclose(fSBC);
 //////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////





  //  ****************************************
  // Pour tester les valeurs
  FILE *fPtest, *fTtest, *fXAtest, *fXBtest, *fXCtest, *fCtest;
  if((fPtest = fopen("fichiers/Ptest.txt", "w+")) == NULL){printf("erreur ouverture fichier P test\n");return 1;}
  if((fTtest = fopen("fichiers/Ttest.txt", "w+")) == NULL){printf("erreur ouverture fichier T test\n");return 1;}
  if((fXAtest = fopen("fichiers/XAtest.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  if((fXBtest = fopen("fichiers/XBtest.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  if((fXCtest = fopen("fichiers/XCtest.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  if((fCtest = fopen("fichiers/Ctest.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  FILE *fVtest, *fEtest;
  if((fVtest = fopen("fichiers/Vtest.txt", "w+")) == NULL){printf("erreur ouverture fichier V test \n");return 1;}
  if((fEtest = fopen("fichiers/Etest.txt", "w+")) == NULL){printf("erreur ouverture fichier E test \n");return 1;}
  for(int j=0; j<Ne_test; j++){
      E=Emin_test+j*he_test;
      fprintf(fEtest, "%.15lf ",E);
  }
  for(int j=0; j<Nv_test; j++){
      V=Vmin_test+j*hv_test;
      fprintf(fVtest, "%.15lf ",V);
  }
  fclose(fVtest);  fclose(fEtest);

  int zone;
  FILE *fZONEtest;
  if((fZONEtest = fopen("fichiers/ZONEtest.txt", "w+")) == NULL){printf("erreur ouverture fichier ZONE test\n");return 1;}
  //  ****************************************
  
  //  $$$$$$$   Ecriture texte pour verif   $$$$$$$$$$
  FILE *file_Nb_ZONE_test;
  if((file_Nb_ZONE_test = fopen("fichiers/Nb_ZONEbench.txt", "w+")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  for(int i=0; i<(Nv_dis-1)*(Ne_dis-1); i++){
    fprintf(file_Nb_ZONE_test, "%d \n",tab_Nb_ZONE[i] );
  }
  fclose(file_Nb_ZONE_test);

  FILE *file_V_dis;
  if((file_V_dis = fopen("fichiers/Vbench_dis.txt", "w+")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  for(int i=0; i<(Nv_dis-1); i++){
    fprintf(file_V_dis, "%.15lf \n",Vmin_dis+i*hv_dis+hv_dis/2 );
  }
  fclose(file_V_dis);
  
  FILE *file_E_dis;
  if((file_E_dis = fopen("fichiers/Ebench_dis.txt", "w+")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  for(int i=0; i<(Ne_dis-1); i++){
    fprintf(file_E_dis, "%.15lf \n",Emin_dis+i*he_dis+he_dis/2 );
  }
  fclose(file_E_dis);
  //  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  double Ntau=1, dx=1;


  // ***** Boucle *****
  for(int i=0; i<Nv_test; i++){
    printf("  - i= %d\n",i );
    V=Vmin_test+i*hv_test;
    for(int j=0; j<Ne_test; j++){
      //printf(" - i= %d, j= %d \n",i,j );
      //printf("   V= %g  E= %g\n",V,E );
      E=Emin_test+j*he_test;
      //err=lec_benchmark_PT(1,V, E, nb_points, SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE, 
      //                     Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, &P, &T, &x, &zone, &V0x, &E0x, &V0y, &E0y);
  
      err=fPTC_METH_cin_phase(methode, 1, 0, V, E, &P, &c, &T, lambda, lambda, Ntau, dx, dt, Ttau, &zone,
                         nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis,
                         SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
      if(err){return err;}
      
      //

      fprintf(fPtest,"%.10lf\n",P);
      fprintf(fTtest,"%.10lf\n",T);
      fprintf(fXAtest,"%.10lf\n",lambda[0]);
      fprintf(fXBtest,"%.10lf\n",lambda[1]);
      fprintf(fXCtest,"%.10lf\n",lambda[2]);
      fprintf(fZONEtest,"%d\n",zone);
      fprintf(fCtest,"%.10lf\n",c);      
    }
  }

  fclose(fPtest);  fclose(fTtest);  fclose(fXAtest); fclose(fXBtest); fclose(fXCtest); fclose(fZONEtest); fclose(fCtest);
  

  // Liberation de la mémoire
  // Liberation
  free(tabVx);   free(tabEx);
  free(tabVy);   free(tabEy);
  
  freetab(tabZONE);  free(tab_Nb_ZONE);

  freetab(SA);  freetab(SB);  freetab(SC);
  freetab(SAB);  freetab(SAC);  freetab(SBC);
  
  return 0;
}
*/
