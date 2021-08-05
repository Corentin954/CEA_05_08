#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Matrice.h"
#include "EOS_etain.h"

/*  Valeur d'entree :
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



int trace_phase(int Nv, int Ne, double Vmin, double Vmax, double Emin, double Emax){

  FILE *fileT1, *fileT2, *fileT3;
  FILE *fileP1, *fileP2, *fileP3;
  FILE *fileS1, *fileS2, *fileS3;
  FILE *fileG1, *fileG2, *fileG3;
  FILE *fileV;
  FILE *fileE;

  if((fileT1 = fopen("T1.txt", "w+")) == NULL){printf("erreur ouverture fichier T\n");return 1;}
  if((fileT2 = fopen("T2.txt", "w+")) == NULL){printf("erreur ouverture fichier T\n");return 1;}
  if((fileT3 = fopen("T3.txt", "w+")) == NULL){printf("erreur ouverture fichier T\n");return 1;}

  if((fileP1 = fopen("P1.txt", "w+")) == NULL){printf("erreur ouverture fichier P\n");return 1;}
  if((fileP2 = fopen("P2.txt", "w+")) == NULL){printf("erreur ouverture fichier P\n");return 1;}
  if((fileP3 = fopen("P3.txt", "w+")) == NULL){printf("erreur ouverture fichier P\n");return 1;}
  
  if((fileS1 = fopen("S1.txt", "w+")) == NULL){printf("erreur ouverture fichier S\n");return 1;}
  if((fileS2 = fopen("S2.txt", "w+")) == NULL){printf("erreur ouverture fichier S\n");return 1;}
  if((fileS3 = fopen("S3.txt", "w+")) == NULL){printf("erreur ouverture fichier S\n");return 1;}
  
  if((fileG1 = fopen("G1.txt", "w+")) == NULL){printf("erreur ouverture fichier G\n");return 1;}
  if((fileG2 = fopen("G2.txt", "w+")) == NULL){printf("erreur ouverture fichier G\n");return 1;}
  if((fileG3 = fopen("G3.txt", "w+")) == NULL){printf("erreur ouverture fichier G\n");return 1;}
  
  if((fileV = fopen("V1.txt", "w+")) == NULL){printf("erreur ouverture fichier V\n");return 1;}

  if((fileE = fopen("E1.txt", "w+")) == NULL){printf("erreur ouverture fichier E\n");return 1;}
  
  
  double he,hv;
  
  he=(Emax-Emin)/(Ne-1);
  hv=(Vmax-Vmin)/(Nv-1);

  double E, V;
  double P1, P2, P3;
  double T1, T2, T3;
  double S1, S2, S3;
  double G1, G2, G3;

  for(int i=0; i<Nv; i++){
    V=Vmin+i*hv;
    fprintf(fileV, "%.15lf ",V );
    for(int j=0; j<Ne; j++){

      E=Emin+i*he;
      fprintf(fileE, "%.15lf ",E );
      
      //phase 1
      P1=fP(V,E, 1); T1=fT(V,E, 1);
      S1=fS(V,E, 1); G1=E+P1*V-T1*S1;

      fprintf(fileT1,"%.15lf ",T1);
      fprintf(fileP1,"%.15lf ",P1);
      fprintf(fileS1,"%.15lf ",S1);
      fprintf(fileG1,"%.15lf ",G1);

      // phase 2
      P2=fP(V,E, 2); T2=fT(V,E, 2);
      S2=fS(V,E, 2); G2=E+P2*V-T2*S2;

      fprintf(fileT2,"%.15lf ",T2);
      fprintf(fileP2,"%.15lf ",P2);
      fprintf(fileS2,"%.15lf ",S2);
      fprintf(fileG2,"%.15lf ",G2);
      
      // phase 3
      P3=fP(V,E, 3); T3=fT(V,E, 3);
      S3=fS(V,E, 3); G3=E+P3*V-T3*S3;

      fprintf(fileT3,"%.15lf ",T3);
      fprintf(fileP3,"%.15lf ",P3);
      fprintf(fileS3,"%.15lf ",S3);
      fprintf(fileG3,"%.15lf ",G3);

    }
    fprintf(fileT1,"\n"); fprintf(fileT2,"\n"); fprintf(fileT3,"\n");
    fprintf(fileP1,"\n"); fprintf(fileP2,"\n"); fprintf(fileP3,"\n");
    fprintf(fileS1,"\n"); fprintf(fileS2,"\n"); fprintf(fileS3,"\n");
    fprintf(fileG1,"\n"); fprintf(fileG2,"\n"); fprintf(fileG3,"\n");
  }

  fclose(fileT1);   fclose(fileT2);   fclose(fileT3);
  fclose(fileP1);   fclose(fileP2);   fclose(fileP3);
  fclose(fileS1);   fclose(fileS2);   fclose(fileS3);
  fclose(fileG1);   fclose(fileG2);   fclose(fileG3);
  
  return 0;

}

/*  Valeur d'entree :
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

  A[dim*1+0]=0.;    A[dim*1+1]=0.;    A[dim*1+2]=-dPvb;    A[dim*1+3]=-dPeb;      

  A[dim*2+0]=dTva;  A[dim*2+1]=dTea;  A[dim*2+2]=-dTvb; A[dim*2+3]=-dTeb;  

  A[dim*3+0]=dGva;  A[dim*3+1]=dGea;  A[dim*3+2]=-dGvb;    A[dim*3+3]=-dGeb; 
  
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
    printf("\n****  i= %d ******\n",i );
    segP[i]=Pdeb+i*hp;
    Pcst=segP[i];
    //printf("  Pcst= %.15lf, hp= %.15lf, Pdeb= %.15lf\n", Pcst, hp, Pdeb );
    
    printf("  ok 1\n");
    // phase A
    Pa=fP(Va,Ea, phaseA); Ta=fT(Va,Ea, phaseA);
    Sa=fS(Va,Ea, phaseA); Ga=Ea+Pa*Va-Ta*Sa;
    //phase B
    printf("  ok 2\n");
    Pb=fP(Vb,Eb, phaseB); Tb=fT(Vb,Eb, phaseB);
    Sb=fS(Vb,Eb, phaseB); Gb=Eb+Pb*Vb-Tb*Sb;
    
    printf("  ok 3\n");

    critere=fabs(Ga-Gb)+fabs(Pa-Pcst)+fabs(Pb-Pcst);
    //printf("  critere= %.14lf\n",critere );
    
    n=1;
    while(critere>=epsilon && n<Nmax){
      printf("  --n= %d\n",n );

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
    printf("  n fin =%d \n",n);
    printf("  ok 4\n");

    segT[i]=Tb;
    segVa[i]=Va;
    segEa[i]=Ea;
    segVb[i]=Vb;
    segEb[i]=Eb;

    printf("  ok 5\n");
    
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

  printf("  ok fin boucle i\n");
  free(A); 
  free(InvA);

}

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

// focntion d'inversion 
// Obtention de (V,E) à partir de (P,T)
// NEwton classique 
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



  ///////////////////////////////////////////////////
  //     Calcul du Point courbe frontières
  ///////////////////////////////////////////////////
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
  double* segVb13=malloc(Np13*sizeof(double));
  double* segEb13=malloc(Np13*sizeof(double));

  double* segP23=malloc(Np23*sizeof(double));
  double* segT23=malloc(Np23*sizeof(double));
  double* segVa23=malloc(Np23*sizeof(double));
  double* segEa23=malloc(Np23*sizeof(double));
  double* segVb23=malloc(Np23*sizeof(double));
  double* segEb23=malloc(Np23*sizeof(double));

  double* segPHAUT=malloc(Nh*sizeof(double));
  double* segPBAS=malloc((Nb1+Nb2)*sizeof(double));
  double* segTGAUCHE=malloc((Ng1+Ng2)*sizeof(double));
  double* segTDROITE=malloc((Nd1+Nd2)*sizeof(double));

  // frontiere 12
  printf("PHASE 12 :\n");
  ligne_double(Nmax, epsilon, phaseA12, phaseB12, Pdeb, Pfin12, Np12, Va012, Ea012, Vb012, Eb012, segP12, segT12, segVa12, segEa12, segVb12, segEb12);

  // frontiere 13
  printf("PHASE 13 :\n");
  ligne_double(Nmax, epsilon, phaseA13, phaseB13, Pdeb, Pfin13, Np13, Va013, Ea013, Vb013, Eb013, segP13, segT13, segVa13, segEa13, segVb13, segEb13);
  
  // frontiere 23
  printf("PHASE 23 :\n");
  ligne_double(Nmax, epsilon, phaseA23, phaseB23, Pdeb, Pfin23, Np23, Va023, Ea023, Vb023, Eb023, segP23, segT23, segVa23, segEa23, segVb23, segEb23);
  
  
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


  /*
  fonction print
  */
  int res= print_double(Np12, segP12, segT12, segVa12, segEa12, segVb12, segEb12,  
                        Np13, segP13, segT13, segVa13, segEa13, segVb13, segEb13,    
                        Np23, segP23, segT23, segVa23, segEa23, segVb23, segEb23,
                        Nb, Nh, Ng, Nd, Tb, Th, Pg, Pd,
                        segPBAS, segPHAUT, segTGAUCHE, segTDROITE);  
   
  if(res){ 
    printf("Erreur dans l'impression"); return 1;}

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
  double Pfin;


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

/////////////////////
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


/////////////////////
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
  free(segP12);  free(segT12);
  free(segVa12);  free(segEa12);
  free(segVb12);  free(segEb12);

  free(segP13);  free(segT13);
  free(segVa13);  free(segEa13);
  free(segVb13);  free(segEb13);

  free(segP23);  free(segT23);
  free(segVa23);  free(segEa23);
  free(segVb23);  free(segEb23);

  free(segPHAUT);  free(segPBAS);
  free(segTGAUCHE);  free(segTDROITE);

}




// Calcul des fonction thermo pour chaque phase
int courb_thermo(void){
  int Nv=200;
  int Ne=200;

  double Vmin=90e-6;
  double Vmax=150e-6;
  double Emin=50e3;
  double Emax=200e3;
  
  
  int err=trace_phase(Nv, Ne, Vmin, Vmax, Emin, Emax);
  if (err){
    return err;
  }
  
  double he,hv;
  double E,V;
  
  he=(Emax-Emin)/(Ne-1);
  hv=(Vmax-Vmin)/(Nv-1);
  /*
  E=200e3;  
  for(int i=0; i<Nv; i++){
    V=Vmin+i*hv;
    printf("P= %.15lf\n",fP(V,E,1)*1e-9 );
  }
  int phase=3;
  V=130e-6;  
  double G;
  for(int i=0; i<Ne; i++){
    E=Emin+i*he;
    G=E+fP(V,E,phase)*V-fS(V,E,phase)*fT(V,E,phase);
    printf("phase= %d, E= %.15lf, T= %.15lf, P= %.15lf, G= %.15lf, S= %.15lf\n",phase, E ,fT(V,E,phase), fP(V,E,phase)*1e-9,G,fS(V,E,phase) );
  }
  */
}

  ///////////////////////////////////////////////////
  //           Calcul du Point TRIPLE
  ///////////////////////////////////////////////////
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





//  !!!!!!!!!!!!!!!!!!!!!    NE CONVERGE PAS   REVOIR LA METHODE  !!!!!!!!!!!!!!!!!!!!!

/*  pour une valeur de V renvoie E qui est sur la ligne frontière de la zone de mélange netre les phases A et B
  !!!!!!   Bien vérifié que ce E existe i.e. V>V_triple ou V<V_triple (voir en fonction des phases d'entrées)  !!!!!!!

  Valeur d'entree :
     - int phaseX, phaseY
     - double Vx

  Valeur sortie
     - double Ex

  Renvoie un int d'erreur si 
    1 : matrice non-inversible
    2 : Vx, phaseX, phaseY  mal choisi
*/
int fEfront(double Vx, int phaseX, int phaseY, double* EXres){
  int dim=3;
  double det_a;
  int n=1;
  int Nmax=1e4;
  double res;
  double epsilon=5e-6;

  // Triangle triple
  double Vta= 129.905671875844291e-6;
  double Eta= 72.168190127265447e3;
  double Vtb= 127.666112185063994e-6;
  double Etb= 99.480256119236515e3;
  double Vtc= 132.538478805442338e-6;
  double Etc= 136.154346228414312e3;

  double * tabVt=malloc(3*sizeof(double));
  double * tabEt=malloc(3*sizeof(double));
  
  tabVt[0]=129.905671875844291e-6; tabVt[1]=127.666112185063994e-6; tabVt[2]=132.538478805442338e-6;
  tabEt[0]=72.168190127265447e3;   tabEt[1]=99.480256119236515e3;   tabEt[2]=136.154346228414312e3;
 
  double* Delta=malloc(dim*sizeof(double));
  double* dX=malloc(dim*sizeof(double));

  double *A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));


  // Vérification que le Ea est bien défini
  if(phaseX==phaseY){printf("Erreur choix phaseX= %d, phaseY= %d (fEfront)\n",phaseX,phaseY ); return 2;}
  // phase B
  if(phaseX==2){
    if (Vx>Vtb){
      printf("Erreur dans le choix de Va > V_tripleA (fEfront)\n");
      return 2;
    }
  }
  // phase A
  if(phaseX==1){
    if (phaseY==2){
      if (Vx>Vta){
        printf("Erreur dans le choix de Vb > V_tripleB (fEfront)\n");
        return 2;
      }
    }
    else if (phaseY==3){
      if (Vx<Vta){
        printf("Erreur dans le choix de Vb < V_tripleB (fEfront)\n");
        return 2;
      }
    }
  }
  // phase C
  if(phaseX==3){
    if (phaseY==2){
      if (Vx>Vtc){
        printf("Erreur dans le choix de Vc > V_tripleC (fEfront)\n");
        return 2;
      }
    }
    else if (phaseY==1){
      if (Vx<Vtc){
        printf("Erreur dans le choix de Vc < V_tripleC (fEfront)\n");
        return 2;
      }
    }
  }
   

  double Ex, Px, Tx, Sx, Gx, dPex, dTex; 
  double Vy, Ey, Py, Ty, Sy, Gy, dPvy, dPey, dTvy, dTey;

  Ex=tabEt[phaseX-1]; 
  Vy=tabVt[phaseY-1]; 
  Ey=tabEt[phaseY-1]; 

  // phase X
  Px=fP(Vx,Ex, phaseX); Tx=fT(Vx,Ex, phaseX);
  Sx=fS(Vx,Ex, phaseX); Gx=Ex+Px*Vx-Tx*Sx;
  dPex=dPdE(phaseX);
  dTex=dTdE(phaseX);
  // phase Y
  Py=fP(Vy,Ey, phaseY); Ty=fT(Vy,Ey, phaseY);
  Sy=fS(Vy,Ey, phaseY); Gy=Ey+Py*Vy-Ty*Sy;
  dPvy=dPdV(Vy,phaseY); dPey=dPdE(phaseY);
  dTvy=dTdV(Vy,phaseY); dTey=dTdE(phaseY);


  double dGex, dGvy, dGey;
  dGex=Vx*dPex - Sx*dTex;

  dGvy=Vy*dPvy - Sy*dTvy; 
  dGey=Vy*dPey - Sy*dTey;

     
  double correc=1./2;

  double critere=1.0;


  while(critere>epsilon && n<Nmax){
    printf("  *** n= %d ***\n",n );


    // phase X
    Px=fP(Vx,Ex, phaseX); Tx=fT(Vx,Ex, phaseX);
    Sx=fS(Vx,Ex, phaseX); Gx=Ex+Px*Vx-Tx*Sx;
    dPex=dPdE(phaseX);
    dTex=dTdE(phaseX);
    // phase Y
    Py=fP(Vy,Ey, phaseY); Ty=fT(Vy,Ey, phaseY);
    Sy=fS(Vy,Ey, phaseY); Gy=Ey+Py*Vy-Ty*Sy;
    dPey=dPdE(phaseY);
    dTey=dTdE(phaseY);

    dGex=Vx*dPex - Sx*dTex;
    dGvy=Vy*dPvy - Sy*dTvy; 
    dGey=Vy*dPey - Sy*dTey;

    printf("  Iteration (fEfront) :\n");
    printf("    Ex= %.15lf\n",Ex );
    printf("    Ey= %.15lf\n",Ey );
    printf("\n");
    printf("    Vx= %.15lf\n",Vx*1e6 );
    printf("    Vy= %.15lf\n",Vy*1e6 );
    printf("\n");
    printf("    Tx= %.15lf\n",Tx );
    printf("    Ty= %.15lf\n",Ty );
    printf("\n");
    printf("    Px= %.15lf (GPa)\n",Px*1e-9 );
    printf("    Py= %.15lf (GPa)\n",Py*1e-9 );
    printf("\n");
    printf("    Sx= %.15lf\n",Sx );
    printf("    Sy= %.15lf\n",Sy );
    printf("\n");
    printf("    Gx= %.15lf\n",Gx );
    printf("    Gy= %.15lf\n",Gy );


    // Delt
    Delta[0]=Px-Py; Delta[1]=Tx-Ty;
    Delta[2]=Gx-Gy;   
     
    
    for(int j=0; j<dim; j++){
      printf("Delta[%d]= %.15lf  ",j, Delta[j] );
    }
    printf("\n");
    


    // on met le second menbre dans la dernière colonne
    //matrice A
    A[dim*0+0]=-dPex; A[dim*0+1]=dPvy; A[dim*0+2]=dPey;   

    A[dim*1+0]=-dTex; A[dim*1+1]=dTvy; A[dim*1+2]=dTey;       

    A[dim*2+0]=-dGex; A[dim*2+1]=dGvy; A[dim*2+2]=dGey;  
    
    
    printf("  * A = *\n");
    for(int i=0; i<dim; i++){
      for(int j=0; j<dim; j++){
        printf("A[%d][%d]= %.15lf  ",i,j, A[dim*i+j] );
      }
      printf("\n");
    }
    printf("\n");  
    


    det_a=inverse_matrice_pivot(A, dim, InvA);
    //printf("detA= %.15lf\n",det_a );
    if ( fabs(det_a)<1e-8 || isnan(det_a) ){printf("A non inversible : detA= %.15lf\n", det_a); return 1;}
    

    printf("  * InvA = *\n");
    for(int i=0; i<dim; i++){
      for(int j=0; j<dim; j++){
        printf("InvA[%d][%d]= %.15lf  ",i,j, InvA[dim*i+j] );
      }
      printf("\n");
    }
    printf("\n");
    
    
    
    // Mutiplcation matricielle de A et InvA
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
    


    // Mutiplcation matricielle
    for(int i=0; i<dim; i++){
      dX[i]=0;
      for(int j=0; j<dim; j++){
        dX[i]+=InvA[dim*i+j]*Delta[j];
      }
      //printf("  dX[%d]= %.15lf\n",i,dX[i] );
    }
    //printf("\n");
     

    Ex=Ex+correc*dX[1];
    Vy=Vy+correc*dX[2];
    Ey=Ey+correc*dX[3];

    critere=fabs(dX[1]);
    printf("critere= %.15lf\n",critere );

    n++;
  }
  
  *EXres=Ex;
  
  free(A); free(InvA);
  free(Delta);
  free(dX);
  return 0;
}



// fonction pour determiner la zone
/*  -int *pind : pointeur sur l'indice (valeur retour)  
*/
int fzone(double V, double E, int* pind){
  double epsilon=1e-6;

  // Points triple
  double VA= 129.905671875844291e-6;
  double EA= 72.168190127265447e3;
  double VB= 127.666112185063994e-6;
  double EB= 99.480256119236515e3;
  double VC= 132.538478805442338e-6;
  double EC= 136.154346228414312e3;

  // Calcul de G
  double Pa,Pb,Pc;
  double Ta,Tb,Tc;
  double Sa,Sb,Sc;
  double Ga,Gb,Gc;
  // phase A
  Pa=fP(V,E, 1); Ta=fT(V,E, 1);
  Sa=fS(V,E, 1); Ga=E+Pa*V-Ta*Sa;
  if(isnan(Ga)){printf("Ga = NaN\n"); Ga=1e30;}
  // phase B
  Pb=fP(V,E, 2); Tb=fT(V,E, 2);
  Sb=fS(V,E, 2); Gb=E+Pb*V-Tb*Sb;
  if(isnan(Gb)){printf("Gb = NaN\n"); Gb=1e30;}
  // phase C
  Pc=fP(V,E, 3); Tc=fT(V,E, 3);
  Sc=fS(V,E, 3); Gc=E+Pc*V-Tc*Sc;
  if(isnan(Gc)){printf("Gc = NaN\n"); Gc=1e30;}

  if(isnan(Ga) && isnan(Gb) && isnan(Gc)){ return 1;}


  // TRIANGLE
  double penteBA, penteBC, penteAC;

  penteBA=(EA-EB)/(VA-VB);
  penteBC=(EC-EB)/(VC-VB);
  penteAC=(EC-EA)/(VC-VA);

  double EBA=EB + penteBA*(V-VB);
  double EAC=EA + penteAC*(V-VA);
  double EBC=EB + penteBC*(V-VB);

  printf("Enthalpie de Gibbs : Ga= %.15lf, Gb= %.15lf, Gc= %.15lf\n",Ga,Gb,Gc);

  
  // triangle
  if(E>EBA && E<EBC && E>EAC){
    printf("triangle : |Ga-Gb|= %.15lf, |Ga-Gc|= %.15lf\n",fabs(Ga-Gb),fabs(Ga-Gc));
    *pind=0;
    return 0;
  } 
  // Zone de mélange
  else{
    if(fabs(Ga-Gb)<epsilon){
      *pind=12;
      printf("Zone de élange A/B : fabs(Ga-Gb)= %.15lf\n",fabs(Ga-Gb));
      return 0;
    }
    else if(fabs(Ga-Gc)<epsilon){
      *pind=13;
      printf("Zone de élange A/C : fabs(Ga-Gc)= %.15lf\n",fabs(Ga-Gc));
      return 0;
    }
    else if(fabs(Gb-Gc)<epsilon){
      *pind=23;
      printf("Zone de élange B/C : fabs(Gb-Gc)= %.15lf\n",fabs(Gb-Gc));
      return 0;
    }
    // PHASE PUR
    else{
      if(Ga<Gb && Ga<Gc){
        *pind=1;
        return 0;
      }
      else if(Gb<Ga && Gb<Gc){
        *pind=2;
        return 0;
      }
      else if(Gc<Ga && Gc<Gb){
        *pind=3;
        return 0;
      }
      else{
        printf("erreur (fzone)\n");
      }
    } 
  }

}



// ZONE 
// repere les zones (triangle pour le moment)
int diag_zone(){
  int err,ind;
  int Nv=200;
  int Ne=200;

  double Vmin=120e-6;
  double Vmax=135e-6;
  double Emin=40e3;
  double Emax=180e3;  

  double he,hv;
  double E,V;
  
  he=(Emax-Emin)/(Ne-1);
  hv=(Vmax-Vmin)/(Nv-1);

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
      E=Emin+j*he;
      printf("V= %.15lf, E=%.15lf\n",V*1e6,E );
      err=fzone(V, E, &ind);
      if(err){return err;}
      fprintf(fZONE, "%d ",ind );
      printf("  ind= %d\n",ind);
    }
    fprintf(fZONE, "\n");
  }
  
  fclose(fZONE);
  fclose(fV);
  fclose(fE);

}




  // Calcul de VE
int routine_VE(void){
  int err;
  int phase;
  double P;
  double T;
  double E,V;

  /*
  printf("PHASE 1 :\n");
  phase=1;
  P=0e9; T=300;
  err= fVE(phase, P, T, &V, &E);
  //err= fVE_NR(phase, P, T, &V, &E);
  printf("P= %.15lf, T= %.15lf\n",P*1e-9,T);
  printf("V= %.15lf, E= %.15lf\n",V*1e6,E);
  printf("fonction inverse :\n");
  printf("fP(V,E)= %.15lf, fT(V,E)= %.15lf\n",fP(V,E,phase),fT(V,E,phase));

  printf("PHASE 2 :\n");
  phase=2;
  P=20e9; T=300;
  err= fVE(phase, P, T, &V, &E);
  //err= fVE_NR(phase, P, T, &V, &E);
  printf("P= %.15lf, T= %.15lf\n",P*1e-9,T);
  printf("V= %.15lf, E= %.15lf\n",V*1e6,E);
  printf("fonction inverse :\n");
  printf("fP(V,E)= %.15lf, fT(V,E)= %.15lf\n",fP(V,E,phase)*1e-9,fT(V,E,phase));

  printf("PHASE 3 :\n");
  phase=3;
  P=0e9; T=1700;
  err= fVE(phase, P, T, &V, &E);
  //err= fVE_NR(phase, P, T, &V, &E);
  printf("P= %.15lf, T= %.15lf\n",P*1e-9,T);
  printf("V= %.15lf, E= %.15lf\n",V*1e6,E);
  printf("fonction inverse :\n");
  printf("fP(V,E)= %.15lf, fT(V,E)= %.15lf\n",fP(V,E,phase)*1e-9,fT(V,E,phase));
  */

  
  FILE* fres1;
  FILE* fres2;
  FILE* fres3;
  if((fres1 = fopen("VEfront1.txt", "w+")) == NULL){printf("erreur ouverture fichier (routine_VE)\n");return 1;}
  if((fres2 = fopen("VEfront2.txt", "w+")) == NULL){printf("erreur ouverture fichier (routine_VE)\n");return 1;}
  if((fres3 = fopen("VEfront3.txt", "w+")) == NULL){printf("erreur ouverture fichier (routine_VE)\n");return 1;}


  int N;
  double Tdeb, Tfin;
  double Pdeb, Pfin;


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

/////////////////////
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


/////////////////////
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

}


// fonction de maillage de la zone A ( qui a un nombre impair de sommets)
/* Arguments :
  - int N12, N13.
  - int NT1 : nombre de point sur la forntière P=0 T=[300, 505]
  - double* P12, T12 : tableaux de N12 éléments qui contients les points de la courbe C12
                       partant du point triple
  - double* P13, T13 : tableaux de N13 éléments qui contients les points de la courbe C13
                       partant du point triple

  Arguments de sorties: void
  - la fonction va écrire dans un fichier les sommets et dans l'autre les mailles 
*/
int maillageA(int N12, int N13, int NT1, double* P12, double* T12, double* P13, double* T13){
  
  /*
      AJOUT DE n_maille != nref
      dans l'écriture des  mailles de la zone interne
      ou pas (pas sur!!)
  */

  int N=2*(N12+N13)-3;
  int n=floor(N/2); //N=2*n+1
  
  int err;
  int phase=1;

  double P,T;
  double V,E;
  
  FILE* fptsPT; FILE* fmaill; FILE* fptsVE;
  if((fptsPT = fopen("ptsAPT.txt", "w+")) == NULL){printf("erreur ouverture fichier pts (maillageA)\n");return 1;}
  if((fptsVE = fopen("ptsAVE.txt", "w+")) == NULL){printf("erreur ouverture fichier pts (maillageA)\n");return 1;}
  if((fmaill = fopen("maillA.txt", "w+")) == NULL){printf("erreur ouverture fichier maill (maillgeA)\n");return 1;}


  double Pref, Tref12, Tref13;

  Pref=P13[N13-1];
  Tref12=T12[N12-1];
  Tref13=T13[N13-1];

  double ht=(Tref12-Tref13)/(NT1-1);
  
  // Liste où on mets le numéro des points pour la zone intermédiaire
  int NG=N13+NT1-1; int ND=N12-1;
  double* LG=malloc((N13+NT1-1)*sizeof(double));
  double* LD=malloc((N12-1)*sizeof(double));


  // **********   ECRITURE DES SOMMETS   *************
  // front en bas à gauche pression de C13 zone rectangle
  for(int j=0; j<NT1; j++){
    for(int i=N13-1; i>=0; i--){
      printf("  --i= %d\n",i );
      P=P13[i];
      T=Tref13+j*ht;
      fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
      err= fVE(phase, P, T, &V, &E);
      if(err){printf("Erreur de fVE (maillageA)\n");return err;}
      fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
    }
  }
  // zone haut gauche pression et temp 13
  for(int j=N13-2; j>=0; j--){
    for(int i=j; i>=0; i--){
      printf("  --i= %d\n",i );
      P=P13[i];
      T=T13[j];
      fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
      err= fVE(phase, P, T, &V, &E);
      if(err){printf("Erreur de fVE (maillageA)\n");return err;}
      fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
    }
  }

  // Zone droite , on commence d'en bas
  for(int j=N12-1; j>=0; j--){
    for(int i=1; i<=j; i++){
      printf("  --i= %d\n",i );
      P=P12[i];
      T=T12[j];
      fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
      err= fVE(phase, P, T, &V, &E);
      if(err){printf("Erreur de fVE (maillageA)\n");return err;}
      fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
    }
  }


  
  // **********   ECRITURE DES MAILLES   *************
  int nl=N13; // nombre d'éléments par ligne
  int MG;     // nombre de mailles à gauche
  int a,b,c;
/* maille triangle : (num des points coresp.)  (sens inverse des aiguilles d'une montre)
   c    b               c
             et           
   a               a    b
*/
  // Zone bas gauche
  for(int l=0; l<=NT1-2; l++){
    for(int k=1; k<=NT1-1; k++){
      a=k;
      b=nl+k+1;
      c=nl+k;
      fprintf(fmaill, "%d %d %d\n",a,b,c );
      a=k;
      b=k+1;
      c=nl+k+1;
      fprintf(fmaill, "%d %d %d\n",a,b,c );
    }

    LG[l]=N13*(l+1);
  }
  
  int nref=N13*NT1+1;  // point de debut de ligne basse (tout à gauche)
  
  // Zone haut gauche
  for(int l=0; l<=N13-2; l++){
    nl=N13-l; // nombre de points de la ligne
    // première maille de la ligne
    a=nref;
    b=nref+1;
    c=nref+nl;
    fprintf(fmaill, "%d %d %d\n",a,b,c );

    for(int k=1; k<=nl-2; k++){
      a=nref+k;
      b=nref+nl+k+1;
      c=nref+nl+k;
      fprintf(fmaill, "%d %d %d\n",a,b,c );
      a=nref+k;
      b=nref+k+1;
      c=nref+nl+k+1;
      fprintf(fmaill, "%d %d %d\n",a,b,c );
    }
    nref+=nl;
    LG[NT1-1+l]=nref-1;;
  }
  
  MG=nref;
  nref+=1;  // le point triple

  // Zone droite
  for(int l=0; l<=NT1-2; l++){
    nl=N12-l;;
    for(int k=1; k<=nl-2; k++){
      a=nref+k;
      b=nref+nl+k+1;
      c=nref+nl+k;
      fprintf(fmaill, "%d %d %d\n",a,b,c );
      a=nref+k;
      b=nref+k+1;
      c=nref+nl+k+1;
      fprintf(fmaill, "%d %d %d\n",a,b,c );
    }
    // dernier trianlge a gauche
    a=nref+nl-2;
    b=nref+nl-1;
    c=nref+nl+(nl-2);
    fprintf(fmaill, "%d %d %d\n",a,b,c );
    
    LD[l]=nref;
    nref+=nl;
  }

  // Zone interne entre les zones gauches et droite
  // Nb pts gauche  :  N13+NT-1     |   Nb pts droite : N12-1
  //  cas 1 :     nb mailles = 2*N12-3
  //  cas 2 :     nb mailles = 
  if(NG>=ND){   // cas 1
  
    for(int k=1; k<=NG-1; k++){
      a=LG[k-1];
      b=LD[k-1];
      c=LG[k];
      fprintf(fmaill, "%d %d %d\n",a,b,c );
      a=LG[k];
      b=LD[k-1];
      c=LD[k];
      fprintf(fmaill, "%d %d %d\n",a,b,c );
    }

    a=LG[N12-2];
    b=LD[N12-2];
    c=LG[N12-1];
    fprintf(fmaill, "%d %d %d\n",a,b,c );
  }
  else{   // cas 2
    for(int k=1; k<=ND-1; k++){
      a=LG[k-1];
      b=LD[k-1];
      c=LG[k];
      fprintf(fmaill, "%d %d %d\n",a,b,c );
      a=LG[k];
      b=LD[k-1];
      c=LD[k];
      fprintf(fmaill, "%d %d %d\n",a,b,c );
    }

    a=LG[N13+NT1-1];
    b=LD[N13+NT1-2];
    c=LD[N13+NT1-1];
    fprintf(fmaill, "%d %d %d\n",a,b,c );

  }
  

  fclose(fptsPT); 
  fclose(fptsVE); 
  fclose(fmaill);

  free(LG); free(LD);
}


// fonction de maillage de la zone B ( qui a un nombre impair de sommets)
/* Arguments :
  - int N12, N13.
  - double* P12, T12 : tableaux de N12 éléments qui contients les points de la courbe C12
                       partant du point triple
  - double* P13, T13 : tableaux de N13 éléments qui contients les points de la courbe C13
                       partant du point triple

  Arguments de sorties: void
  - la fonction va écrir dans un efichier les sommets et dans l'autre les mailles 
*/
int maillageB(int N12, int N23, int NP2, double* P12, double* T12, double* P23, double* T23){

  int N=2*(N23)-1;
  int n=floor(N/2); //N=2*n+

  int err;
  int phase=2;

  double P,T;
  double V,E;
  
  FILE* fptsPT; FILE* fmaill; FILE* fptsVE;
  if((fptsPT = fopen("ptsBPT.txt", "w+")) == NULL){printf("erreur ouverture fichier pts (maillageA)\n");return 1;}
  if((fptsVE = fopen("ptsBVE.txt", "w+")) == NULL){printf("erreur ouverture fichier pts (maillageA)\n");return 1;}
  if((fmaill = fopen("maillB.txt", "w+")) == NULL){printf("erreur ouverture fichier fmaill (maillgeA)\n");return 1;}


  double Pref12, Tref12;
  double Pref23, Tref23;

  Tref12=T12[N12-1];
  Tref23=T23[N23-1];

  Pref12=P12[N12-1];
  Pref23=P23[N23-1];

  int nj, nref,ind;

  double h=(Pref23-Pref12)/(NP2-1);


  int NH=N23-1;
  int NB=N23+NP2-1;

  double* LH=malloc(NH*sizeof(double));
  double* LB=malloc(NB*sizeof(double));


  // **********   ECRITURE DES SOMMETS   *************
  // Parcours de les lignes a T=cst de bas en haut et de les pressions de gauche à droit
  // jusqu'au point triple
  ind=0;
  nj=NP2;
  for(int j=0; j<=N12-1; j++){
    T=T12[N12-1-j];
    // zone avant Pref12
    for(int Iind=ind-1; Iind>=0; Iind--){
      P=P12[N12-1-Iind];
      fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
      err= fVE(phase, P, T, &V, &E);
      if(err){printf("Erreur de fVE (maillageB)\n");return err;}
      fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
    }
    // zone apres Pref12
    for(int i=0; i<NP2; i++){
      P=Pref12+h*i;
      fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
      err= fVE(phase, P, T, &V, &E);
      if(err){printf("Erreur de fVE (maillageB)\n");return err;}
      fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
    }
    ind++;
  }

  // Zone haute avec C23 a gauche
  nj=N23-1; // nombre de points sur la ligne j
  for(int j=1; j<=N23-1; j++){
    T=T23[j];
    for(int i=j; i<=N23-1; i++){
      P=P23[i];
      fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
      err= fVE(phase, P, T, &V, &E);
      if(err){printf("Erreur de fVE (maillageB)\n");return err;}
      fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
    }
  }



  

// **********   ECRITURE DES MAILLES   *************

  int a,b,c;
  int nl;

// MAILLE DU BAS
/* maille triangle : (num des points coresp.)
   c    b              c
             et         
   a              a    b
*/
  nref=1;
  nl=NP2;
  for(int l=0; l<=N12-2; l++){  // N12-1 lignes
    a=nref;
    b=nref+nl+1;
    c=nref+nl;
    fprintf(fmaill, "%d %d %d\n",a,b,c);
    for(int i=1; i<=nl-1; i++){
      a=nref+i;
      b=nref+nl+i+1;
      c=nref+nl+i;
      fprintf(fmaill, "%d %d %d\n",a,b,c);
      a=nref+i;
      b=nref+i+1;
      c=nref+nl+i+1;
      fprintf(fmaill, "%d %d %d\n",a,b,c);
    }
    nref+=nl;
    nl++;
  }


  for(int i=0; i<=N12+NP2-2; i++){
    LB[i]=nref+i;
  }

  
  nref+=N12+NP2-1;

  for(int i=0; i<=N23-2; i++){
    LH[i]=nref+i;
  }


  // MAILLE DU HAUT
/* maille triangle : (num des points coresp.)
   c    b              c
             et         
   a              a    b
*/
  nl=N23-1;
  for(int l=0; l<=N23-2; l++){  // N23-1 lignes
    a=nref;
    b=nref+1;
    c=nref+nl;
    fprintf(fmaill, "%d %d %d\n",a,b,c);
    for(int i=1; i<=nl; i++){
      a=nref+i;
      b=nref+nl+i+1;
      c=nref+nl+i;
      fprintf(fmaill, "%d %d %d\n",a,b,c);
      a=nref+i;
      b=nref+i+1;
      c=nref+nl+i+1;
      fprintf(fmaill, "%d %d %d\n",a,b,c);
    }
    nref+=nl;
    nl--;
  }
  

  // Zone interne entre les zones hautes et basses
  // Nb pts gauche  :  N13+NT-1     |   Nb pts droite : N12-1
  //  cas 1 :     nb mailles = 2*N12-3
  //  cas 2 :     nb mailles = 

  if(NB<=NH){   // cas 1
  
    for(int k=0; k<=NB-2; k++){
      a=LB[k];
      b=LH[k+1];
      c=LH[k];
      fprintf(fmaill, "%d %d %d\n",a,b,c );
      a=LB[k];
      b=LB[k+1];
      c=LH[k+1];
      fprintf(fmaill, "%d %d %d\n",a,b,c );
    }

    a=LB[NB-1];
    b=LH[NH-1];
    c=LH[NB-1];
    fprintf(fmaill, "%d %d %d\n",a,b,c );
  }
  else{   // cas 2
    for(int k=0; k<=NH-2; k++){
      a=LB[k];
      b=LH[k+1];
      c=LH[k];
      fprintf(fmaill, "%d %d %d\n",a,b,c );
      a=LB[k];
      b=LB[k+1];
      c=LH[k+1];
      fprintf(fmaill, "%d %d %d\n",a,b,c );
    }

    a=LB[NH-1];
    b=LB[NB-1];
    c=LH[NH-1];
    fprintf(fmaill, "%d %d %d\n",a,b,c );

  }


  fclose(fptsPT); fclose(fptsVE);
  fclose(fmaill);

  free(LH); free(LB);
}

// fonction de maillage de la zone B ( qui a un nombre pair de sommets)
/* Arguments :
  - int N12, N13.
  - double* P12, T12 : tableaux de N12 éléments qui contients les points de la courbe C12
                       partant du point triple
  - double* P13, T13 : tableaux de N13 éléments qui contients les points de la courbe C13
                       partant du point triple

  Arguments de sorties: void
  - la fonction va écrir dans un efichier les sommets et dans l'autre les mailles 
*/
int maillageC(int N13, int N23, int NT3, double* P13, double* T13, double* P23, double* T23){

  int N=2*(N13+N23-1);
  int n=floor(N/2); //N=2*n+

  int err;
  int phase=3;

  double P,T;
  double V,E;

  FILE* fptsPT; FILE* fmaill; FILE* fptsVE;
  if((fptsPT = fopen("ptsCPT.txt", "w+")) == NULL){printf("erreur ouverture fichier pts (maillageC)\n");return 1;}
  if((fptsVE = fopen("ptsCVE.txt", "w+")) == NULL){printf("erreur ouverture fichier pts (maillageC)\n");return 1;}
  if((fmaill = fopen("maillC.txt", "w+")) == NULL){printf("erreur ouverture fichier fmaill (maillgeC)\n");return 1;}


  double Tref23, Trefhaut;
  Trefhaut=2500;
  Tref23=T23[N23-1];

  double h=(Trefhaut-Tref23)/(NT3-1);


  // **********   ECRITURE DES SOMMETS   *************
  // Zone basse gauche , parcours de C13 jausuq'au point triple (dans le sens inverse)
  for(int j=1; j<=N13; j++){ // N13 points
    T=T13[N13-j];
    for(int i=1; i<=j; i++){
      P=P13[N13-i];
      fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
      err= fVE(phase, P, T, &V, &E);
      if(err){printf("Erreur de fVE (maillageC)\n");return err;}
      fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
    }
  }

  // Zone hauteur médiane , parcours de C23 jausuq'à la frontière droite (dans le bon sens)
  for(int j=1; j<=N23-1; j++){
    T=T23[j];
    // zone au dessus de C13
    for(int i=N13-1; i>=0; i--){
      P=P13[i];
      fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
      err= fVE(phase, P, T, &V, &E);
      if(err){printf("Erreur de fVE (maillageC)\n");return err;}
      fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
    }
    // zone au dessus de C23
    for(int i=1; i<=j; i++){
      P=P23[i];
      fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
      err= fVE(phase, P, T, &V, &E);
      if(err){printf("Erreur de fVE (maillageC)\n");return err;}
      fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
    }
  }


  // Zone haute 
  for(int j=1; j<=NT3-1; j++){
    T=Tref23+h*j;
    // zone au dessus de C13
    for(int i=N13-1; i>=0; i--){
      P=P13[i];
      fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
      err= fVE(phase, P, T, &V, &E);
      if(err){printf("Erreur de fVE (maillageC)\n");return err;}
      fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
    }
    // zone au dessus de C23
    for(int i=1; i<=N23-1; i++){
      P=P23[i];
      fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
      err= fVE(phase, P, T, &V, &E);
      if(err){printf("Erreur de fVE (maillageC)\n");return err;}
      fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
    }
  }
  

  // **********   ECRITURE DES MAILLES   *************
  int a,b,c;

  int nref=1; 

  for(int l=1; l<=N13+N23-2; l++){
    for(int i=1; i<=l-1; i++){
      a=nref+(i-1);
      b=nref+l+i;
      c=nref+l+(i-1);
      fprintf(fmaill, "%d %d %d\n",a,b,c);
      a=nref+(i-1);
      b=nref+i; 
      c=nref+l+i;
      fprintf(fmaill, "%d %d %d\n",a,b,c);
    }

    a=nref+l-1;
    b=nref+2*l; 
    c=nref+2*l-1;
    fprintf(fmaill, "%d %d %d\n",a,b,c);

    nref+=l;
  }

  


  int nl=N13+N23-1;
  for(int l=1; l<=NT3-1; l++){
    for(int i=1; i<=N13+N23-2; i++){
      a=nref+(i-1);
      b=nref+nl+i;
      c=nref+nl+(i-1);
      fprintf(fmaill, "%d %d %d\n",a,b,c);
      a=nref+(i-1);
      b=nref+i; 
      c=nref+nl+i;
      fprintf(fmaill, "%d %d %d\n",a,b,c);
    }
    nref+=nl;
  }



  fclose(fptsPT); fclose(fptsVE); 
  fclose(fmaill);
}




/*  Creation des fichiers de maillage

*/
int maillage(void){
  /*
  int err;
  int Nmax=1e2;
  double epsilon=1e-7;
  
  int N12=10;
  int N13=10;
  int N23=10;
  
  double* P12; double* T12; 
  double* P13; double* T13; 
  double* P23; double* T23;
  
  P12=malloc(N12*sizeof(double));
  T12=malloc(N12*sizeof(double));

  P13=malloc(N13*sizeof(double));
  P13=malloc(N13*sizeof(double));
  
  P23=malloc(N23*sizeof(double));
  P23=malloc(N23*sizeof(double));


  // pas utile ici
  double* segVx, *segEx, *segVy, *segEy;
  segVx=malloc((N12+N13+N23)*sizeof(double));
  segEx=malloc((N12+N13+N23)*sizeof(double));
  segVy=malloc((N12+N13+N23)*sizeof(double));
  segEy=malloc((N12+N13+N23)*sizeof(double));

  // Points triple
  double VA= 129.905671875844291e-6;
  double EA= 72.168190127265447e3;
  double VB= 127.666112185063994e-6;
  double EB= 99.480256119236515e3;
  double VC= 132.538478805442338e-6;
  double EC= 136.154346228414312e3;

  double Pdeb=4.5104217816091e9;  // pression du point triple
  double Pfin12=9.4e9;
  double Pfin13=0.;
  double Pfin23=70e9;

  printf("ligne double A :\n");
  ligne_double(Nmax, epsilon, 1, 2, Pdeb, Pfin12, N12, VA, EA, VB, EB, P12, T12, segVx, segEx, segVy, segEy);
  printf("ligne double B :\n");
  ligne_double(Nmax, epsilon, 1, 3, Pdeb, Pfin13, N13, VA, EA, VC, EC, P13, T13, segVx, segEx, segVy, segEy);
  printf("ligne double C :\n");
  ligne_double(Nmax, epsilon, 2, 3, Pdeb, Pfin23, N23, VB, EB, VC, EC, P23, T23, segVx, segEx, segVy, segEy);
  */

  int err;
  int Nmax=1e4;
  double epsilon=5e-5;
  double nb_points=4;
  int NT1=nb_points, NP2=3*nb_points, NT3=floor(nb_points/2);
  int N12=nb_points; int N13=nb_points; int N23=5*nb_points;
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
  double* segVb13=malloc(N13*sizeof(double));
  double* segEb13=malloc(N13*sizeof(double));

  double* segP23=malloc(N23*sizeof(double));
  double* segT23=malloc(N23*sizeof(double));
  double* segVa23=malloc(N23*sizeof(double));
  double* segEa23=malloc(N23*sizeof(double));
  double* segVb23=malloc(N23*sizeof(double));
  double* segEb23=malloc(N23*sizeof(double));

  double* segPHAUT=malloc(Nh*sizeof(double));
  double* segPBAS=malloc((Nb1+Nb2)*sizeof(double));
  double* segTGAUCHE=malloc((Ng1+Ng2)*sizeof(double));
  double* segTDROITE=malloc((Nd1+Nd2)*sizeof(double));

  // frontiere 12
  printf("PHASE 12 :\n");
  ligne_double(Nmax, epsilon, phaseA12, phaseB12, Pdeb, Pfin12, N12, Va012, Ea012, Vb012, Eb012, segP12, segT12, segVa12, segEa12, segVb12, segEb12);

  // frontiere 13
  printf("PHASE 13 :\n");
  ligne_double(Nmax, epsilon, phaseA13, phaseB13, Pdeb, Pfin13, N13, Va013, Ea013, Vb013, Eb013, segP13, segT13, segVa13, segEa13, segVb13, segEb13);
  
  // frontiere 23
  printf("PHASE 23 :\n");
  ligne_double(Nmax, epsilon, phaseA23, phaseB23, Pdeb, Pfin23, N23, Va023, Ea023, Vb023, Eb023, segP23, segT23, segVa23, segEa23, segVb23, segEb23);
  
  
  

  printf("maillage A\n");
  err=maillageA(N12, N13, NT1, segP12, segT12, segP13, segT13);
  if(err){return err;}
  
  printf("maillage B\n");
  err=maillageB(N12, N23, NP2, segP12, segT12, segP23, segT23);
  if(err){return err;}

  printf("maillage C\n");
  err=maillageC(N13, N23, NT3, segP13, segT13, segP23, segT23);
  if(err){return err;}



    // Liberation
  free(segP12);  free(segT12);
  free(segVa12);  free(segEa12);
  free(segVb12);  free(segEb12);

  free(segP13);  free(segT13);
  free(segVa13);  free(segEa13);
  free(segVb13);  free(segEb13);

  free(segP23);  free(segT23);
  free(segVa23);  free(segEa23);
  free(segVb23);  free(segEb23);

  free(segPHAUT);  free(segPBAS);
  free(segTGAUCHE);  free(segTDROITE);


  return 0;
}




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
  err=maillage();



}
