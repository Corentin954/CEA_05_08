#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "Matrice.h"
#include "EOS_etain.h"
#include "head_etain.h"

double inverse_matrice_pivot(double *mat, int dim, double *mat_inv); /* détermine la matrice inverse d'une matrice carrée non nulle par la méthode du pivot*/


//  Calcul d'une isentrope 
/* 2 isentropes : 1 et 2
*/
int isentrope(){

  int err;
  int N_dis=2e4;

  int sch_cin_phase=0;
  
  double V1,E1,P1,T1;
  double V2,E2,P2,T2;
  double P01,T01;
  double P02,T02;
  P01=0; T01=50;
  P02=1e5; T02=300;
  double V01,E01;
  double V02,E02;

  err=fVE(1, P01, T01, &V01, &E01);
  if(err){return 1;}
  err=fVE(1, P02, T02, &V02, &E02);
  if(err){return 1;}

  //printf("P0= %g, T0= %g, V0= %g, E0= %g\n",P0,T0,V0,E0 );

  double Vfin=90e-6;
  double hv1=(Vfin-V01)/(N_dis-1);
  double hv2=(Vfin-V02)/(N_dis-1);

  /*
  double tabV=malloc(N_dis*sizeof(double));  double tabE=malloc(N_dis*sizeof(double));
  double tabP=malloc(N_dis*sizeof(double));  double tabT=malloc(N_dis*sizeof(double));
  */

  double C1,C2, dX, dt, Ttau; 
  double* lambda_test=malloc(3*sizeof(double));
  int  Ntau, zone; 
  

  
  // ******* Calcul des frontières ***********


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


  double hv_dis=(Vmax_dis-Vmin_dis)/(Nv_dis-1);
  double he_dis=(Emax_dis-Emin_dis)/(Ne_dis-1);

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





  //  ****************************************
  // Pour plot les isentropes
  FILE *fP, *fT, *fXA, *fXB, *fXC, *fC;
  if((fP = fopen("fichiers/Pisen.txt", "w+")) == NULL){printf("erreur ouverture fichier P test\n");return 1;}
  if((fT = fopen("fichiers/Tisen.txt", "w+")) == NULL){printf("erreur ouverture fichier T test\n");return 1;}
  if((fC = fopen("fichiers/Cisen.txt", "w+")) == NULL){printf("erreur ouverture fichier C test\n");return 1;}
  if((fXA = fopen("fichiers/XAisen.txt", "w+")) == NULL){printf("erreur ouverture fichier XA test\n");return 1;}
  if((fXB = fopen("fichiers/XBisen.txt", "w+")) == NULL){printf("erreur ouverture fichier XB test\n");return 1;}
  if((fXC = fopen("fichiers/XCisen.txt", "w+")) == NULL){printf("erreur ouverture fichier XC test\n");return 1;}
  FILE *fV, *fE;
  if((fV = fopen("fichiers/Visen.txt", "w+")) == NULL){printf("erreur ouverture fichier V test \n");return 1;}
  if((fE = fopen("fichiers/Eisen.txt", "w+")) == NULL){printf("erreur ouverture fichier E test \n");return 1;}



  // *****************************************

  

  V1=V01; E1=E01;
  P1=P01; T1=T01;

  V2=V02; E2=E02;
  P2=P02; T2=T02;

  fprintf(fV, "%g %g\n",V1,V2 );  fprintf(fE, "%g %g\n",E1,E2 );

  for(int i=1; i<N_dis; i++){
    printf("--i= %d\n",i );
    V1=V01+i*hv1;
    E1-=P1*hv1;

    V2=V02+i*hv2;
    E2-=P2*hv2;
    //printf(" V=%g; E=%g;\n",V,E );
    
    fprintf(fP, "%g %g\n",P1,P2 );  fprintf(fT, "%g %g\n",T1,T2 );  fprintf(fC, "%g %g\n",C1,C2 );
    fprintf(fXA, "%g\n",lambda_test[0] );   fprintf(fXB, "%g\n",lambda_test[1] );   fprintf(fXC, "%g\n",lambda_test[2] ); 

    err=fPTC_cin_phase(0, 0, sch_cin_phase, V1, E1, &P1, &C1, &T1, lambda_test, lambda_test, Ntau, dX, dt, Ttau, &zone,
                       nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis,
                       SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
    if(err){return err;}

    err=fPTC_cin_phase(0, 0, sch_cin_phase, V2, E2, &P2, &C2, &T2, lambda_test, lambda_test, Ntau, dX, dt, Ttau, &zone,
                       nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis,
                       SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
    if(err){return err;}

    fprintf(fV, "%g %g\n",V1,V2 );  fprintf(fE, "%g %g\n",E1,E2 );
    


  }
  fprintf(fP, "%g %g\n",P1,P2 );  fprintf(fT, "%g %g\n",T1,T2 );  fprintf(fC, "%g %g\n",C1,C2 );
  fprintf(fXA, "%g\n",lambda_test[0] );   fprintf(fXB, "%g\n",lambda_test[1] );   fprintf(fXC, "%g\n",lambda_test[2] ); 




  return 0;
}


//  HUGONIOT


/* Valeur d'entree :
    _ phaseX, phaseY
    - Vx,Ex, Vy,Ey
    - Vpts, Epts : plot d'Hogoniot
    - x : la fraction massique de phase Y est imposé (ici 0 ou 1)

    - A,InvA,Delta,dX : pour limiter les allocations 


  Valeur sortie
    - Vx,Ex, Vy,Ey, 


*/
int  Newton_hugniot_pt12(int cas, double* A, double* InvA, double* Delta, int phaseX, int phaseY, double Vpt, double Ept, double Ppt, double Vx, double Ex, double Vy, double Ey, double* pVX, double* pEX, double* pVY, double* pEY, double* dX){
  int dim=4;

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
  if (cas==1){      Delta[3]=2*(Ept-Ex) + (Ppt+Px)*(Vpt-Vx) ; }
  else if(cas==2){  Delta[3]=2*(Ept-Ey) + (Ppt+Py)*(Vpt-Vy) ; }
   
  /*
  for(int j=0; j<dim; j++){
    printf("Delta[%d]= %g  ",j, Delta[j] );
  }
  printf("\n");
  */
  
  // Dérive des equation Hugoniot
  double dHvx, dHex;
  double dHvy, dHey;

  if(cas==1){
    dHvx= Ppt - Vpt*dPvx + Px + Vx*dPvx;
    dHex= 2.0 - Vpt*dPex + Vx*dPex;
    dHvy= 0;
    dHey= 0;
  }
  else if(cas==2){   
    dHvx= 0;
    dHex= 0;
    dHvy= Ppt - Vpt*dPvy + Py + Vy*dPvy;
    dHey= 2.0 - Vpt*dPey + Vy*dPey;
  }
  
  
  // on met le second menbre dans la dernière colonne
  //matrice A
  A[dim*0+0]=dPvx; A[dim*0+1]=dPex; A[dim*0+2]=-dPvy; A[dim*0+3]=-dPey;   
  
  A[dim*1+0]=dTvx; A[dim*1+1]=dTex; A[dim*1+2]=-dTvy; A[dim*1+3]=-dTey; 
  
  A[dim*2+0]=dGvx; A[dim*2+1]=dGex; A[dim*2+2]=-dGvy; A[dim*2+3]=-dGey; 
  
  A[dim*3+0]=dHvx; A[dim*3+1]=dHex; A[dim*3+2]=dHvy;  A[dim*3+3]=dHey;     
  
  /*
  printf("  * A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){  printf("A[%d][%d]= %g  ",i,j, A[dim*i+j] );  }  printf("\n");
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
      printf("InvA[%d][%d]= %g  ",i,j, InvA[dim*i+j] );
    }
    printf("\n");
  }
  printf("\n");
  */
  
  
  // Print de la mutiplcation matricielle de A et InvA
  /*
  double res;
  printf("  * InvA*A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      res=0;
      for(int k=0; k<dim; k++){
        res+=InvA[dim*i+k]*A[dim*k+j];
      }
      printf("A*InvA[%d][%d]= %g ",i,j,res );
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
    //printf("  dX[%d]= %g\n",i,dX[i] );
  }
  //printf("\n");
   
  double correc=1./50;

  *pVX=Vx+correc*dX[0];
  *pEX=Ex+correc*dX[1];
  *pVY=Vy+correc*dX[2];
  *pEY=Ey+correc*dX[3];

  
  //printf("  Vx= %g, Ex= %g ,  Vy= %g, Ey= %g \n",*pVX,*pEX,*pVY,*pEY );

  return 0;
}


// CAS=
//     1 : lance la calcul du point 1 avec Newton_hugniot_pt1
//     2 : lance la calcul du point 2 avec Newton_hugniot_pt2
int hugoniot_pts12(int cas, int Nmax, double epsilon, int phaseX, int phaseY, double Vpt, double Ept, double Ppt, double V0x, double E0x, double V0y, double E0y, double* pVx, double* pEx, double* pVy, double* pEy){
  int n=1;
  int dim=4;
  int err;
  
  double* A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));
  double* dX=malloc(dim*sizeof(double));
  double* Delta=malloc(dim*sizeof(double));
  
  double Px, Tx, Sx, Gx;
  double Py, Ty, Sy, Gy;
  
  double Vx, Ex, Vy, Ey;
  double VX, EX, VY, EY;
  
  Vx=V0x; Ex=E0x; 
  Vy=V0y; Ey=E0y; 
  
  
  // phase X
  Px=fP(Vx,Ex, phaseX); Tx=fT(Vx,Ex, phaseX);
  Sx=fS(Vx,Ex, phaseX); Gx=Ex+Px*Vx-Tx*Sx;
  // phase Y
  Py=fP(Vy,Ey, phaseY); Ty=fT(Vy,Ey, phaseY);
  Sy=fS(Vy,Ey, phaseY); Gy=Ey+Py*Vy-Ty*Sy;
  
  
  double critere=1.;
  
  //printf("  critere= %.14lf\n",critere );
    
  n=0;
  while(critere>=epsilon && n<Nmax){
    
    err= Newton_hugniot_pt12(cas, A, InvA, Delta, phaseX, phaseY, Vpt, Ept, Ppt, Vx, Ex, Vy, Ey, &VX, &EX, &VY, &EY, dX); 
    if(err){return 1;}    

    Vx=VX; Ex=EX;
    Vy=VY; Ey=EY;

    //critere=fabs(Gx-Gy)+fabs(Px-Py)+fabs(Ty-Tx);
    critere=fabs(dX[0]/Vx) +fabs(dX[1]/Ex) +fabs(dX[2]/Vy) +fabs(dX[3]/Ey);
    //critere=fabs(dX[0]) +fabs(dX[1]) +fabs(dX[2]) +fabs(dX[3])+fabs(dX[4]);
    //printf(" --n= %d  critere= %g\n",n,critere );

    
    n++;
  }
  printf(" --n= %d, critere= %g\n",n,critere );
  //for(int i=0; i<dim; i++){  printf(" dX[%d]= %g  ",i,dX[i] );  } printf("\n");
  
  if(n==Nmax){printf("Newton_mixte non convergent : -n= %d , critere= %g\n",n,critere );return 1;}

  *pVx=Vx; *pEx=Ex;
  *pVy=Vy; *pEy=Ey;


  free(A); 
  free(InvA);
  free(dX);
  free(Delta);
  
  return 0;
}


int  Newton_hugniot_pt3(double* A, double* InvA, double* Delta, int phase, double* Vpt, double* Ept, int* PHASEpt, double V, double E, double* pV, double* pE, double* dX){
  int dim=2;

  // Grandeurs thermodynamiques
  double P, dPv, dPe;
  // phase X
  P=fP(V,E, phase); 
  dPv=dPdV(V, phase); dPe=dPdE(phase);

  double* Ppt=malloc(2*sizeof(double));
  for(int i=0; i<2; i++){ Ppt[i]=fP(Vpt[i],Ept[i], PHASEpt[i]); } 
  
  double D1=Vpt[0]*sqrt((Ppt[1]-Ppt[0])/(Vpt[0]-Vpt[1]));
  //double D1=Vpt[0]*sqrt((Ppt[1]-Ppt[0])/(Vpt[0]-Vpt[1]));
  
  /*
  printf("Px= %.10lf, Tx= %.10lf, Sx= %.10lf, Gx= %.10lf\n",Px,Tx,Sx,Gx );
  printf("Py= %.10lf, Ty= %.10lf, Sy= %.10lf, Gy= %.10lf\n",Py,Ty,Sy,Gy );
  */
  
  double u1=Vpt[0]*(Ppt[1]-Ppt[0])/D1;
  double DuV2=((D1-u1)/Vpt[1])*((D1-u1)/Vpt[1]);
  double DPDV=((Ppt[1]-Ppt[0])/(Vpt[1]-Vpt[0]));
  printf("u1= %g\n",u1 );
  
  Delta[0]=2*(Ept[0]-E) + (Ppt[0]+P)*(Vpt[0]-V) ;
  //Delta[1]=Ppt[1]-P + DuV2*(V-Vpt[1]) ;
  Delta[1]=Ppt[1]-P + DPDV*(V-Vpt[1]) ;


  /*
  for(int j=0; j<dim; j++){
    printf("Delta[%d]= %g  ",j, Delta[j] );
  }
  printf("\n");
  */
  
  // on met le second menbre dans la dernière colonne
  //matrice A
  A[dim*0+0]= (V-Vpt[0])*dPv + Ppt[0] + P ;  A[dim*0+1]= 2.0 + (V- Vpt[0])*dPe ; 
  
  //A[dim*1+0]= dPv-DuV2 ;  A[dim*1+1]= dPe ; 
  A[dim*1+0]= dPv-DPDV ;  A[dim*1+1]= dPe ; 
  
  /*
  printf("  * A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){  printf("A[%d][%d]= %g  ",i,j, A[dim*i+j] );  }  printf("\n");
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
      printf("InvA[%d][%d]= %g  ",i,j, InvA[dim*i+j] );
    }
    printf("\n");
  }
  printf("\n");
  */
  
  
  // Print de la mutiplcation matricielle de A et InvA
  /*
  double res;
  printf("  * InvA*A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      res=0;
      for(int k=0; k<dim; k++){
        res+=InvA[dim*i+k]*A[dim*k+j];
      }
      printf("A*InvA[%d][%d]= %g ",i,j,res );
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
    //printf("  dX[%d]= %g\n",i,dX[i] );
  }
  //printf("\n");
   
  double correc=1.;

  *pV=V+correc*dX[0];
  *pE=E+correc*dX[1];
  
  //printf("  V= %g, E= %g \n",*pV,*pE );
  
  free(Ppt);
  return 0;
}


int hugniot_pt3(int Nmax, double epsilon, int phase, double* Vpt, double* Ept, int* PHASEpt, double V0, double E0, double* pV, double* pE){
  int n=0;
  int dim=2;
  int err;

  double* A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));
  double* dX=malloc(dim*sizeof(double));
  double* Delta=malloc(dim*sizeof(double));
  
  double P;

  double V, E;
  double Vin, Ein;

  Vin=V0; Ein=E0; 

  // phase X
  P=fP(V,E, phase); 


  double critere;
  
  critere=1.;
  //printf("  critere= %.14lf\n",critere );
    
  n=1;
  while(critere>=epsilon && n<Nmax){
    //printf("  --n= %d\n",n );

    err= Newton_hugniot_pt3(A, InvA, Delta, phase, Vpt, Ept, PHASEpt, Vin, Ein, &V, &E, dX);
    if(err){return 1;}    

    Vin=V; Ein=E;

    //critere=fabs(Gx-Gy)+fabs(Px-Py)+fabs(Ty-Tx);
    critere=fabs(dX[0]/V) +fabs(dX[1]/E);
    //critere=fabs(dX[0]) +fabs(dX[1]) +fabs(dX[2]) +fabs(dX[3])+fabs(dX[4]);
    printf(" --n= %d  critere= %g, dV/V= %g, dE/E= %g\n",n,critere,dX[0]/V,dX[1]/E );

    
    n++;
  }
  //printf(" --n= %d, critere= %g\n",n,critere );
  //for(int i=0; i<dim; i++){  printf(" dX[%d]= %g  ",i,dX[i] );  } printf("\n");
  
  if(n==Nmax){printf("Newton_mixte non convergent : -n= %d , critere= %g\n",n,critere );return 1;}

  *pV=V; *pE=E;
  printf(" P3= %g, T3= %g\n", fP(V,E,2), fT(V,E,2) );


  free(A); 
  free(InvA);
  free(dX);
  free(Delta);
  
  return 0;
}


// Discretisation de courbes sur phase pur en cherchant E
int courbe_hugoniot_pure_E(int Nmax, double epsilon, int N, int phase, double V0, double E0, double Vpt, double Ept, double Ppt, double Vfin, double *tabV, double* tabE, double* tabP, double* tabT){
  // portion de courbe 0-1
  int n=0;
  double f,f_prime,dE;
  double V,E;
  V=V0; E=E0;
  double hv=(Vfin-V0)/N;
  double critere=1.;
  tabV[0]=V0;  tabE[0]=E0;
  tabP[0]=fP(V0,E0,phase);  tabT[0]=fT(V0,E0,phase);

  for(int i=1; i<N; i++){
    n=0;  critere=1.;
    V=tabV[i-1];  E=tabE[i-1];
    V+=hv;
    tabV[i]=V;
    while(critere>epsilon && n<Nmax){
      f=E-Ept - fP(V,E,phase)*(Vpt-V)/2;
      f_prime=1-dPdE(phase)*(Vpt-V)/2;

      dE=-f/f_prime;
      critere=fabs(dE/E);
      E+=dE;
      //printf(" --n= %d , dE= %g\n",n,dE );
      n++;
    }
    if(n==Nmax){printf("Newton 1D non convergent : critere= %g (courbe_hugoniot_pure)\n",critere ); return 1;}
    tabE[i]=E;
    tabP[i]=fP(V,E,phase);
    tabT[i]=fT(V,E,phase);
    printf(" **i= %d, critere= %g, V= %g, E= %g, P= %g\n",i,critere,V,E,tabP[i] );
  }

  
  return 0;
}


// Discretisation de courbes sur phase pur en cherchant P
int courbe_hugoniot_pure_P(int Nmax, double epsilon, int N, int phase, double V0, double P0, double Vpt, double Ept, double Ppt, double Vfin, double *tabV, double* tabE, double* tabP, double* tabT){
  // portion de courbe 0-1
  int n=0;
  double f,f_prime,dP;
  double V,P;
  V=V0; P=P0;
  double hv=(Vfin-V0)/N;
  double critere=1.;
  tabV[0]=V0;  tabP[0]=P0;
  tabE[0]=fE_VP(V0,P0,phase);  tabT[0]=fT(V0,tabE[0],phase);  

  double correc=1./10;

  for(int i=1; i<N; i++){
    n=0;  critere=1.;
    V=tabV[i-1];  
    P=tabP[i-1];  
    V+=hv;
    tabV[i]=V;
    while(critere>epsilon && n<Nmax){
      f=fE_VP(V,P,phase)-Ept - P*(Vpt-V)/2;
      f_prime=(1./dPdE(phase)) -(Vpt-V)/2;

      dP=correc*(-f/f_prime);
      critere=fabs(dP/P);
      P+=dP;
      //printf(" --n= %d , dE= %g\n",n,dE );
      n++;
    }
    if(n==Nmax){printf("Newton 1D non convergent : critere= %g (courbe_hugoniot_pure)\n",critere ); return 1;}
    tabP[i]=P;
    tabE[i]=fE_VP(V,P,phase); 
    tabT[i]=fT(V,tabE[i],phase);
    printf(" **i= %d, n= %d, critere= %g, V= %g, P= %g, E= %g\n",i,n,critere,V,P,tabE[i] );
  }

  
  return 0;
}


int Newton_hugoniot_mixte(double* A, double* InvA, double* Delta, double* dX, double Pcst, int phaseX, int phaseY, double Vpt, double Ept, double Ppt, double Vx, double Ex, double Vy, double Ey, double* VX, double* EX, double* VY, double* EY, double* px){
  int dim=5;

  // Grandeurs thermodynamiques
  double Px, Tx, Sx, Gx, dPvx, dPex, dTvx, dTex;
  double Py, Ty, Sy, Gy, dPvy, dPey, dTvy, dTey;
  // phase A
  Px=fP(Vx,Ex, phaseX); Tx=fT(Vx,Ex, phaseX);
  Sx=fS(Vx,Ex, phaseX); Gx=Ex+Px*Vx-Tx*Sx;
  dPvx=dPdV(Vx,phaseX); dPex=dPdE(phaseX);
  dTvx=dTdV(Vx,phaseX); dTex=dTdE(phaseX);
  // phase B
  Py=fP(Vy,Ey, phaseY); Ty=fT(Vy,Ey, phaseY);
  Sy=fS(Vy,Ey, phaseY); Gy=Ey+Py*Vy-Ty*Sy;
  dPvy=dPdV(Vy,phaseY); dPey=dPdE(phaseY);
  dTvy=dTdV(Vy,phaseY); dTey=dTdE(phaseY);


  double dGvx, dGex, dGvy, dGey;
  dGvx=Vx*dPvx - Sx*dTvx;
  dGex=Vx*dPex - Sx*dTex;

  dGvy=Vy*dPvy - Sy*dTvy; 
  dGey=Vy*dPey - Sy*dTey;
  
  double P2=(Pcst+Ppt)/2;

  double x=*px;

  double Vmel=(1-x)*Vx + x*Vy;
  double Emel=(1-x)*Ex + x*Ey;

  Delta[0]=Px-Pcst; Delta[1]=Py-Pcst;
  Delta[2]=Ty-Tx;   Delta[3]=Gy-Gx;
  Delta[4]=Ept - Emel + P2*(Vpt - Vmel);
   
  
  //for(int j=0; j<dim; j++){
  //  printf("Delta[%d]= %.15lf  ",j, Delta[j] );
  //}
  //printf("\n");
  


  // on met le second menbre dans la dernière colonne
  //matrice A
  A[dim*0+0]=-dPvx;     A[dim*0+1]=-dPex; A[dim*0+2]=0.;    A[dim*0+3]=0.;    A[dim*0+4]=0.;
 
  A[dim*1+0]=0.;        A[dim*1+1]=0.;    A[dim*1+2]=-dPvy; A[dim*1+3]=-dPey; A[dim*1+4]=0.;

  A[dim*2+0]=dTvx;      A[dim*2+1]=dTex;  A[dim*2+2]=-dTvy; A[dim*2+3]=-dTey; A[dim*2+4]=0.;

  A[dim*3+0]=dGvx;      A[dim*3+1]=dGex;  A[dim*3+2]=-dGvy; A[dim*3+3]=-dGey; A[dim*3+4]=0.;

  A[dim*4+0]=P2*(1-x); A[dim*4+1]=(1-x); A[dim*4+2]=P2*x; A[dim*4+3]=x;     A[dim*4+4]=Ey - Ex + P2*(Vy-Vx);
  
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
    //printf("  dX[%d]= %g\n",i,dX[i] );
  }
  //printf("\n");
   
  double correc=1.;

  *VX=Vx+correc*dX[0];
  *EX=Ex+correc*dX[1];
  *VY=Vy+correc*dX[2];
  *EY=Ey+correc*dX[3];
  *px=x +correc*dX[4];

  return 0;
}

// Discretisation de courbes sur phase pur
int courbe_hugoniot_mixte(int Nmax, double epsilon, int N, int phaseX, int phaseY, double Pdeb, double Pfin, double Vpt, double Ept, double Ppt, double V0x, double E0x, double V0y, double E0y, double *tabV, double* tabE, double* tabP, double* tabT, double* tabX){
  // portion de courbe 0-1
  int n=0;
  
  int dim=5;
  int err;

  double* A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));
  double* dX=malloc(dim*sizeof(double));
  double* Delta=malloc(dim*sizeof(double));

  double Vx,Ex,Vy,Ey;
  double VX,EX,VY,EY;
  double x;

  double Vmel, Emel;
  double hp=(Pfin-Pdeb)/N;
  double critere=1.;
  tabV[0]=V0x;  tabE[0]=E0x;
  tabP[0]=Pdeb;  tabT[0]=fT(V0x,E0x,phaseX);
  tabX[0]=0;
  double Pcst=Pdeb;
  
  Vx=V0x;   Vy=V0y; 
  Ex=E0x;   Ey=E0y; 

  for(int i=1; i<N; i++){
    n=0;  critere=1.;
    Pcst+=hp;
    tabP[i]=Pcst;
    x=1./2;
    while(critere>epsilon && n<Nmax){
      // Ici dim=4   // MODIFER Newton_double pour liliter les allocations
      //Newton_double(A, InvA, Delta, dX, Pcst, phaseX, phaseY, Vx, Ex, Vy, Ey, &VX, &EX, &VY, &EY);
      err=Newton_hugoniot_mixte(A, InvA, Delta, dX, Pcst, phaseX, phaseY, Vpt, Ept, Ppt, Vx, Ex, Vy, Ey, &VX, &EX, &VY, &EY, &x);
      if(err){return 1;}
      //printf(" --n= %d , critere= %g\n",n,critere );
      Vx=VX;   Vy=VY; 
      Ex=EX;   Ey=EY; 

      critere=(fabs(dX[0]/Vx) + fabs(dX[1]/Ex) + fabs(dX[2]/Vy) + fabs(dX[3]/Ey) + fabs(dX[4]/x))/dim;
      //printf("*i= %d -n= %d, critere= %g, Vx= %g, Ex= %g, Vy= %g, Ey= %g, x= %g\n",i,n,critere, Vx,Ex,Vy,Ey,x );
      n++;
    }
    if(n==Nmax){printf("Newton non convergent : critere= %g (courbe_hugoniot_pure)\n",critere ); return 1;}     
    
    Vmel=(1.-x)*VX + x*VY;
    Emel=(1.-x)*EX + x*EY;
    
    tabV[i]=Vmel;
    tabE[i]=Emel;
    tabT[i]=fT(Vx,Ex,phaseX);
    tabX[i]=x;
    printf(" **i= %d, critere= %g, Vmel= %g, Emel= %g, x= %g,  P= %g\n",i,critere,Vmel,Emel,x,tabP[i] );
  }
  

  free(A); free(InvA);
  free(Delta); free(dX);
  
  return 0;
}

int Newton_hugoniot_mixte_Xfixe(double* A, double* InvA, double* Delta, double* dX, double Xcst, int phaseX, int phaseY, double Vpt, double Ept, double Ppt, double Vx, double Ex, double Vy, double Ey, double* VX, double* EX, double* VY, double* EY){
  int dim=4;

  // Grandeurs thermodynamiques
  double Px, Tx, Sx, Gx, dPvx, dPex, dTvx, dTex;
  double Py, Ty, Sy, Gy, dPvy, dPey, dTvy, dTey;
  // phase A
  Px=fP(Vx,Ex, phaseX); Tx=fT(Vx,Ex, phaseX);
  Sx=fS(Vx,Ex, phaseX); Gx=Ex+Px*Vx-Tx*Sx;
  dPvx=dPdV(Vx,phaseX); dPex=dPdE(phaseX);
  dTvx=dTdV(Vx,phaseX); dTex=dTdE(phaseX);
  // phase B
  Py=fP(Vy,Ey, phaseY); Ty=fT(Vy,Ey, phaseY);
  Sy=fS(Vy,Ey, phaseY); Gy=Ey+Py*Vy-Ty*Sy;
  dPvy=dPdV(Vy,phaseY); dPey=dPdE(phaseY);
  dTvy=dTdV(Vy,phaseY); dTey=dTdE(phaseY);


  double dGvx, dGex, dGvy, dGey;
  dGvx=Vx*dPvx - Sx*dTvx;
  dGex=Vx*dPex - Sx*dTex;

  dGvy=Vy*dPvy - Sy*dTvy; 
  dGey=Vy*dPey - Sy*dTey;
  
  double P=(Px+Py)/2;
  P=Px;
  double P2=(P+Ppt)/2;

  double x=Xcst;

  double Vmel=(1-x)*Vx + x*Vy;
  double Emel=(1-x)*Ex + x*Ey;

  Delta[0]=Py-Px;  Delta[1]=Ty-Tx;   Delta[2]=Gy-Gx;
  Delta[3]=Ept - Emel + P2*(Vpt - Vmel);
   
  
  //for(int j=0; j<dim; j++){
  //  printf("Delta[%d]= %.15lf  ",j, Delta[j] );
  //}
  //printf("\n");
  


  // on met le second menbre dans la dernière colonne
  //matrice A
  A[dim*0+0]=dPvx;     A[dim*0+1]=dPex;  A[dim*0+2]=-dPvy; A[dim*0+3]=-dPey;

  A[dim*1+0]=dTvx;     A[dim*1+1]=dTex;  A[dim*1+2]=-dTvy; A[dim*1+3]=-dTey;

  A[dim*2+0]=dGvx;     A[dim*2+1]=dGex;  A[dim*2+2]=-dGvy; A[dim*2+3]=-dGey; 

  A[dim*3+0]=P2*(1-x) + ((Vmel-Vpt)/2)*dPvx; A[dim*3+1]=(1-x) + ((Vmel-Vpt)/2)*dPex; A[dim*3+2]=P2*x;  A[dim*3+3]=x;     
  
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
  if ( fabs(det_a)<1e-8){printf("A non inversible : detA= %g\n", det_a); return 1;}
  
  /*
  printf("  * InvA = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      printf("InvA[%d][%d]= %g  ",i,j, InvA[dim*i+j] );
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
      printf("A*InvA[%d][%d]= %g ",i,j,res );
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
    //printf("  dX[%d]= %g ",i,dX[i] );
  }
  //printf("\n");
   
  double correc=1.;

  *VX=Vx+correc*dX[0];
  *EX=Ex+correc*dX[1];
  *VY=Vy+correc*dX[2];
  *EY=Ey+correc*dX[3];

  return 0;
}


// Discretisation de courbes sur phase pur
int Pt_hugoniot_mixte_Xfixe(int Nmax, double epsilon, int N, int phaseX, int phaseY, double* tabX, double Vpt, double Ept, double Ppt, double V0x, double E0x, double V0y, double E0y, double* tabV, double* tabE, double* tabP, double* tabT){
  // portion de courbe 0-1
  int n=0;
  
  int dim=4;
  int err;

  double* A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));
  double* dX=malloc(dim*sizeof(double));
  double* Delta=malloc(dim*sizeof(double));

  double Vx,Ex,Vy,Ey;
  double VX,EX,VY,EY;

  double Vmel, Emel, x;
  double critere;
  
  Vx=V0x;   Vy=V0y; 
  Ex=E0x;   Ey=E0y; 
  
  for (int i = 0; i <N; i++){
    x=tabX[i];
    Vx=V0x;   Vy=V0y; 
    Ex=E0x;   Ey=E0y; 
    critere=1.; n=0;
    while(critere>epsilon && n<Nmax){
      // Ici dim=4   // MODIFER Newton_double pour liliter les allocations
      err=Newton_hugoniot_mixte_Xfixe(A, InvA, Delta, dX, x, phaseX, phaseY, Vpt, Ept, Ppt, Vx, Ex, Vy, Ey, &VX, &EX, &VY, &EY);
      if(err){return 1;}
      //printf(" --n= %d , critere= %g\n",n,critere );
      Vx=VX;   Vy=VY; 
      Ex=EX;   Ey=EY; 

      critere=(fabs(dX[0]/Vx) + fabs(dX[1]/Ex) + fabs(dX[2]/Vy) + fabs(dX[3]/Ey) )/dim;
      //printf("    -n= %d, x= %g, critere= %g, Vx= %g, Ex= %g, Vy= %g, Ey= %g\n",n,x,critere, Vx,Ex,Vy,Ey );
      n++;
    }
    if(n==Nmax){printf("Newton non convergent : critere= %g (courbe_hugoniot_pure)\n",critere ); return 1;}     
    
    Vmel=(1.-x)*VX + x*VY;
    Emel=(1.-x)*EX + x*EY;
    
    tabV[i]=Vmel;
    tabE[i]=Emel;
    tabP[i]=fP(Vx,Ex,phaseX);
    tabT[i]=fT(Vx,Ex,phaseX);

    printf(" --i= %d, x= %g : n_iter= %d, critere= %g, Vmel= %g, Emel= %g, P= %g , T= %g\n",i,x,n,critere,Vmel,Emel,tabP[i],tabT[i] );
  }

  free(A); free(InvA);
  free(Delta); free(dX);
  
  return 0;
}

//  Calcul d'ue courbe d'Hugoniot
int hugoniot(){

  int err;
  int N, Nx, Nmax=1e2;
  double epsilon=1e-7;
  int phase, phaseX, phaseY;
  double Vfin;
  double Pdeb, Pfin;
  double Vdeb,Edeb; 
  double Vinit, Einit;
  double D,U,u;

  // Point 0
  double P0,T0; 
  P0=1.01325e5; T0=300;
  P0=0.; T0=300;
  int zone0=1;

  double V0,E0;
  err=fVE(1, P0, T0, &V0, &E0);
  if(err){return 1;}
  
  double V1,E1, V2,E2, V3,E3, V4,E4, V5,E5;
  double P1,T1, P2,T2, P3,T3, P4,T4, P5,T5;

  double V0x, E0x, V0y, E0y;
  double Vx, Ex, Vy, Ey;

  double E,V;

  // Valeur d'initialisation
  double VaT=129.905671875844291e-6;  double EaT=72.168190127265447e3;
  double VbT=127.666112185063994e-6;  double EbT=99.480256119236515e3;
  double VcT=132.538478805442338e-6;  double EcT=136.154346228414312e3;


  // CALCUL DU POINT 1
  /*
  V0x=123.8e-6;  E0x=51.4e3;
  V0y=121.9e-6;  E0y=79e3;
  
  Nmax=1e4;  epsilon=4e-14;
  err= hugoniot_pts12(1, Nmax, epsilon, 1, 2, V0, E0, P0, V0x, E0x, V0y, E0y, &Vx, &Ex, &Vy, &Ey);
  if(err){return 1;}
  
  V1=Vx; E1=Ex;
  P1=fP(V1,E1,1);  T1=fT(V1,E1,1);
  printf("  V1= %g,E1= %g , P1= %g, T1= %g\n",V1,E1,P1,T1 );


  // CALCUL DU POINT 2
  V0x=122.5e-6;  E0x=49.6e3;
  V0y=120.6e-6;  E0y=77e3;

  Nmax=1e4;  epsilon=2e-14;
  err= hugoniot_pts12(2, Nmax, epsilon, 1, 2, V1, E1, P1, V0x, E0x, V0y, E0y, &Vx, &Ex, &Vy, &Ey);
  if(err){return 1;}
  
  V2=Vy; E2=Ey;
  P2=fP(V2,E2,2);  T2=fT(V2,E2,2);
  printf("  V2= %g, E2= %g , P2= %g, T2= %g\n",V2,E2,P2,T2 );
  
  // CALCUL DU POINT 3
  double* Vpt, *Ept; int* PHASEpt;
  Vpt=malloc(2*sizeof(double));  Ept=malloc(2*sizeof(double));
  PHASEpt=malloc(2*sizeof(int));  
  
  Vpt[0]=V0; Vpt[1]=V1;  Ept[0]=E0; Ept[1]=E1;
  PHASEpt[0]=1;  PHASEpt[1]=1;

  phase=2;
  Vinit=126e-6;  Einit=100e3;
  
  Nmax=1e3; epsilon=5e-15;
  err= hugniot_pt3(Nmax, epsilon, phase, Vpt, Ept, PHASEpt, Vinit, Einit, &V3, &E3);
  P3=fP(V3,E3,2);  T3=fT(V3,E3,2);
  printf("  V0= %g, E0= %g , P0= %g, T0= %g\n",V0,E0,P0,T0 );
  printf("  V1= %g, E1= %g , P1= %g, T1= %g\n",V1,E1,P1,T1 );
  printf("  V2= %g, E2= %g , P2= %g, T2= %g\n",V2,E2,P2,T2 );
  printf("  V3= %g, E3= %g , P3= %g, T3= %g\n",V3,E3,P3,T3 );
  
  free(Vpt);  free(Ept);  free(PHASEpt);
  */

  // CALCUL DU POINT 4  debut du mélange liquide
  V0x=96e-6;  E0x=8.85e5;
  V0y=95e-6;  E0y=1e6;
  
  Nmax=1e4;  epsilon=4e-14;
  err= hugoniot_pts12(1, Nmax, epsilon, 2, 3, V0, E0, P0, V0x, E0x, V0y, E0y, &Vx, &Ex, &Vy, &Ey);
  if(err){return 1;}  

  V4=Vx; E4=Ex;
  P4=fP(V4,E4,2);  T4=fT(V4,E4,2);
  printf("  V4= %g,E4= %g , P4= %g, T4= %g\n",V4,E4,P4,T4 );

  double* lambda_test=malloc(3*sizeof(double)); // tableau des fractions massiques
  lambda_test[0]=0.0; 
  lambda_test[1]=1.0; 
  lambda_test[2]=0.0; 

  double Va0, Ea0, Vb0, Eb0, Vc0, Ec0; 
  Va0=125e-6;  Ea0=2.5e4;
  Vb0=108e-6;  Eb0=3.4e5;
  Vc0=112e-6;  Ec0=6.8e5;

  int zone=2;
  double dV=1e-8;
  double C4,P;

  int c_celerite_son_cin_3_phase_all_dom(int Nmax, double epsilon, int zone, double P, double V, double E, double dV, 
                             double Va0, double Ea0, double Vb0, double Eb0, double Vc0, double Ec0, double* lambda, double *pc);


  err=c_celerite_son_cin_3_phase_all_dom(Nmax, epsilon, zone, P, V4, E4, dV, Va0, Ea0, Vb0, Eb0, Vc0, Ec0, lambda_test, &C4);  
  if(err){return 1;}
  printf("\n**********\n CELERIT2 C= %g, V=%g, E= %g, P= %g \n**********\n",C4,V4,E4,P );
  free(lambda_test);
  

  
  /*
  // CALCUL DU POINT 5  debut du mélange liquide
  V0x=95e-6;  E0x=1.02e6;
  V0y=94.2e-6;  E0y=1.22e6;
  
  Nmax=1e4;  epsilon=4e-14;
  err= hugoniot_pts12(2, Nmax, epsilon, 2, 3, V0, E0, P0, V0x, E0x, V0y, E0y, &Vx, &Ex, &Vy, &Ey);
  if(err){return 1;}
  
  V5=Vy; E5=Ey;
  P5=fP(V5,E5,3);  T5=fT(V5,E5,3);
  printf("  V5= %g,E5= %g , P5= %g, T5= %g\n",V5,E5,P5,T5 );


  // Calcul de point dans les zones mixtes à fraction massique fixee

  // Point sur la zone mixte 1/2 avec x fixee
  Nx=11;
  double* tabXm12=malloc(Nx*sizeof(double));  
  double* tabVm12=malloc(Nx*sizeof(double));  double* tabEm12=malloc(Nx*sizeof(double));
  double* tabPm12=malloc(Nx*sizeof(double));  double* tabTm12=malloc(Nx*sizeof(double));

  for(int i=0; i<Nx; ++i){
    tabXm12[i]=i*(1./(Nx-1));
  }

  V0x=123.8e-6;  E0x=51.4e3;
  V0y=121.9e-6;  E0y=79e3;
  
  Nmax=1e4;  epsilon=4e-14;
  phaseX=1; phaseY=2;
  printf("\n\n mixte 12\n");
  err=Pt_hugoniot_mixte_Xfixe(Nmax, epsilon, Nx, phaseX, phaseY, tabXm12, V1, E1, P1, V0x, E0x, V0y, E0y, tabVm12, tabEm12, tabPm12, tabTm12);
  

  // Point sur la zone mixte 2/3 avec x fixee
  Nx=11;
  double* tabXm23=malloc(Nx*sizeof(double));  
  double* tabVm23=malloc(Nx*sizeof(double));  double* tabEm23=malloc(Nx*sizeof(double));
  double* tabPm23=malloc(Nx*sizeof(double));  double* tabTm23=malloc(Nx*sizeof(double));

  for(int i=0; i<Nx; ++i){
    tabXm23[i]=i*(1./(Nx-1));
  }

  V0x=V4;  E0x=E4;
  V0y=94.2e-6;  E0y=1.22e6;
  Nmax=2e2;  epsilon=1e-13; 

  phaseX=2; phaseY=3;
  printf("\n\n mixte 23\n");
  err=Pt_hugoniot_mixte_Xfixe(Nmax, epsilon, Nx, phaseX, phaseY, tabXm23, V0, E0, P0, V0x, E0x, V0y, E0y, tabVm23, tabEm23, tabPm23, tabTm23);
  

  double D1=V0*sqrt((P1-P0)/(V0-V1));
  double D2=D1;
  double D3=V0*sqrt((P3-P0)/(V0-V3));
  double D4=V0*sqrt((P4-P0)/(V0-V4));
  double D5=V0*sqrt((P5-P0)/(V0-V5));

  double u1=V0*(P1-P0)/D1;
  double u2=u1+sqrt((P2-P1)*(V1-V2));
  double u3=V0*(P3-P0)/D3;
  double u4=V0*(P4-P0)/D4;
  double u5=V0*(P5-P0)/D5;


  printf("\n\n\n");
  printf("  rho0= %.6lf, V0= %.6lf, E0= %.6lf , P0= %.6lf, T0= %.6lf \n",1./V0, V0*1e6,E0,P0*1e-9,T0 );
  printf("  rho1= %.6lf, V1= %.6lf, E1= %.6lf , P1= %.6lf, T1= %.6lf , D1= %.6lf, u1= %.6lf\n",1./V1,V1*1e6,E1,P1*1e-9,T1,D1,u1 );
  printf("  rho2= %.6lf, V2= %.6lf, E2= %.6lf , P2= %.6lf, T2= %.6lf , D2= %.6lf, u2= %.6lf\n",1./V2,V2*1e6,E2,P2*1e-9,T2,D2,u2 );
  printf("  rho3= %.6lf, V3= %.6lf, E3= %.6lf , P3= %.6lf, T3= %.6lf , D3= %.6lf, u3= %.6lf\n",1./V3,V3*1e6,E3,P3*1e-9,T3,D3,u3 );
  printf("  rho4= %.6lf, V4= %.6lf, E4= %.6lf , P4= %.6lf, T4= %.6lf , D4= %.6lf, u4= %.6lf\n",1./V4,V4*1e6,E4,P4*1e-9,T4,D4,u4 );
  printf("  rho4= %.6lf, V5= %.6lf, E5= %.6lf , P5= %.6lf, T5= %.6lf , D5= %.6lf, u5= %.6lf\n",1./V5,V5*1e6,E5,P5*1e-9,T5,D5,u5 );
  
  printf("\n\n ZONE MIXTE beta/gamma (1/2)\n");
  for(int i=0; i<Nx; i++){
    u=u1 + sqrt((tabPm12[i]-P1)*(V1-tabVm12[i]));
    printf("  x= %g, P= %g, u= %g, rho= %g, V= %g, E= %g, T= %g\n",tabXm12[i], tabPm12[i], u, 1./tabVm12[i], tabVm12[i], tabEm12[i], tabTm12[i] );
  }

  printf("\n\n ZONE MIXTE gamma/liquide (2/3)\n");
  for(int i=0; i<Nx; i++){
    D=V0*sqrt((tabPm23[i]-P0)/(V0-tabVm23[i]));
    u=V0*(tabPm23[i]-P0)/D;
    printf(" x= %g, P= %g, u= %g, rho= %g, V= %g, E= %g, T= %g\n",tabXm23[i],tabPm23[i],u, 1./tabVm23[i], tabVm23[i], tabEm23[i], tabTm23[i] );
  }
  */

  /*
  //****** Déclaration des fichiers pour plot octave  **********
  // Pour tester les valeurs
  FILE *fPh, *fTh, *fXA, *fXB, *fXC, *fC;
  if((fPh = fopen("fichiers/Phug.txt", "w+")) == NULL){printf("erreur ouverture fichier P test\n");return 1;}
  if((fTh = fopen("fichiers/Thug.txt", "w+")) == NULL){printf("erreur ouverture fichier T test\n");return 1;}
  if((fXA = fopen("fichiers/XAhug.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  if((fXB = fopen("fichiers/XBhug.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  if((fXC = fopen("fichiers/XChug.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  if((fC = fopen("fichiers/Chug.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  FILE *fV, *fE, *fU, *fD;
  if((fV = fopen("fichiers/Vhug.txt", "w+")) == NULL){printf("erreur ouverture fichier V test \n");return 1;}
  if((fE = fopen("fichiers/Ehug.txt", "w+")) == NULL){printf("erreur ouverture fichier E test \n");return 1;}
  if((fU = fopen("fichiers/Uhug.txt", "w+")) == NULL){printf("erreur ouverture fichier U test \n");return 1;}
  if((fD = fopen("fichiers/Dhug.txt", "w+")) == NULL){printf("erreur ouverture fichier D test \n");return 1;}

  // ********************************************************


  // DISCRETISATION DES PORTIONS DE COURBES ENTRE LES =! POINTS
  N=2e2;
  double* tabV=malloc(N*sizeof(double));  double* tabE=malloc(N*sizeof(double));
  double* tabP=malloc(N*sizeof(double));  double* tabT=malloc(N*sizeof(double));
  double* tabX=malloc(N*sizeof(double));  
  double* tabU=malloc(N*sizeof(double));  double* tabD=malloc(N*sizeof(double));  
  
  
  // ZONE ENTRE 0 ET 1  (phase 1)
  
  Nmax=2e3; epsilon=1e-13; 
  phase=1;
  Vfin=V1;
  Vdeb=V0; Edeb=E0;
  err=courbe_hugoniot_pure_E(Nmax, epsilon, N, phase, Vdeb, Edeb, V0, E0, P0, Vfin, tabV, tabE, tabP, tabT);
  if(err){return 1;}

  for(int i=0; i<N; i++){
    D=V0*sqrt((tabP[i]-P0)/(V0-tabV[i]));
    U=V0*(tabP[i]-P0)/D;

    fprintf(fV, "%g \n",tabV[i] );
    fprintf(fE, "%g \n",tabE[i] );
    fprintf(fPh, "%g \n",tabP[i] );
    fprintf(fTh, "%g \n",tabT[i] );
    
    fprintf(fU, "%g \n",U );
    fprintf(fD, "%g \n",D );

    fprintf(fXA, "%g \n",1.0 );
    fprintf(fXB, "%g \n",0. );
    fprintf(fXC, "%g \n",0. ); 
  }
  

  // ZONE BIPHASIQUE ENTRE 1 ET 2
  V0x=V1;  E0x=E1;
  V0y=121.9e-6;  E0y=79e3;
  Nmax=2e3; epsilon=5e-13; 
  phaseX=1;  phaseY=2;
  Pdeb=P1; Pfin=P2;
  err= courbe_hugoniot_mixte(Nmax, epsilon, N, phaseX, phaseY, Pdeb, Pfin, V1, E1, P1, V0x, E0x, V0y, E0y,  tabV, tabE, tabP, tabT, tabX);
  if(err){return 1;}
  
  for(int i=0; i<N; i++){
    D=0;
    U=0;

    fprintf(fV, "%g \n",tabV[i] );
    fprintf(fE, "%g \n",tabE[i] );
    fprintf(fPh, "%g \n",tabP[i] );
    fprintf(fTh, "%g \n",tabT[i] );

    fprintf(fU, "%g \n",U );
    fprintf(fD, "%g \n",D );

    fprintf(fXA, "%g \n",1-tabX[i] );
    fprintf(fXB, "%g \n",tabX[i] );
    fprintf(fXC, "%g \n",0. );
  }
  

  
   
  // ZONE ENTRE 2 ET 3 (phase 2)
  Nmax=2e3; epsilon=1e-13; 
  phase=2;
  Vfin=V3;
  Vdeb=V2; Edeb=E2; Pdeb=P2;
  err=courbe_hugoniot_pure_E(Nmax, epsilon, N, phase, Vdeb, Edeb, V2, E2, P2, Vfin, tabV, tabE, tabP, tabT);
  //err=courbe_hugoniot_pure_P(Nmax, epsilon, N, phase, Vdeb, Pdeb, V2, E2, P2, Vfin, tabV, tabE, tabP, tabT);
  if(err){return 1;}

  for(int i=0; i<N; i++){
    D=u1+V1*sqrt((tabP[i]-P1)/(V1-tabV[i]));
    u=u1+V1*sqrt((tabP[i]-P1)*(V1-tabV[i]));

    fprintf(fV, "%g \n",tabV[i] );
    fprintf(fE, "%g \n",tabE[i] );
    fprintf(fPh, "%g \n",tabP[i] );
    fprintf(fTh, "%g \n",tabT[i] );

    fprintf(fU, "%g \n",U );
    fprintf(fD, "%g \n",D );

    fprintf(fXA, "%g \n",0. );
    fprintf(fXB, "%g \n",1.0 );
    fprintf(fXC, "%g \n",0.); 
  }
  


  
  // ZONE ENTRE 3 et 4  (phase 2)
  Nmax=2e3; epsilon=1e-13; 
  phase=2;
  Vfin=V4;
  Vdeb=V3; Edeb=E3;
  err=courbe_hugoniot_pure_E(Nmax, epsilon, N, phase, Vdeb, Edeb, V0, E0, P0, Vfin, tabV, tabE, tabP, tabT);
  if(err){return 1;}

  for(int i=0; i<N; i++){
    D=V0*sqrt((tabP[i]-P0)/(V0-tabV[i]));
    U=V0*(tabP[i]-P0)/D;


    fprintf(fV, "%g \n",tabV[i] );
    fprintf(fE, "%g \n",tabE[i] );
    fprintf(fPh, "%g \n",tabP[i] );
    fprintf(fTh, "%g \n",tabT[i] );
    
    fprintf(fU, "%g \n",U );
    fprintf(fD, "%g \n",D );

    fprintf(fXA, "%g \n",0. );
    fprintf(fXB, "%g \n",1.0 );
    fprintf(fXC, "%g \n",0.); 
  }

  // ZONE BIPHASIQUE ENTRE 4 ET 5
  V0x=V4;  E0x=E4;
  V0y=94.2e-6;  E0y=1.22e6;
  Nmax=2e3;  epsilon=1e-13; 
  phaseX=2;  phaseY=3;
  Pdeb=P4;   Pfin=P5;
  err= courbe_hugoniot_mixte(Nmax, epsilon, N, phaseX, phaseY, Pdeb, Pfin, V0, E0, P0, V0x, E0x, V0y, E0y,  tabV, tabE, tabP, tabT, tabX);
  if(err){return 1;}
  
  for(int i=0; i<N; i++){
    D=V0*sqrt((tabP[i]-P0)/(V0-tabV[i]));
    U=V0*(tabP[i]-P0)/D;

    fprintf(fV, "%g \n",tabV[i] );
    fprintf(fE, "%g \n",tabE[i] );
    fprintf(fPh, "%g \n",tabP[i] );
    fprintf(fTh, "%g \n",tabT[i] );   

    fprintf(fU, "%g \n",U );
    fprintf(fD, "%g \n",D );

    fprintf(fXA, "%g \n",0. );
    fprintf(fXB, "%g \n",1-tabX[i] );
    fprintf(fXC, "%g \n",tabX[i] );
  }
  


  // ZONE APRES 5  (phase 3)
  Nmax=2e3; epsilon=1e-13; 
  phase=3;
  Vfin=90e-6;
  Vdeb=V5;  Edeb=E5;
  err=courbe_hugoniot_pure_E(Nmax, epsilon, N, phase, Vdeb, Edeb, V0, E0, P0, Vfin, tabV, tabE, tabP, tabT);
  if(err){return 1;}

  for(int i=0; i<N; i++){
    D=V0*sqrt((tabP[i]-P0)/(V0-tabV[i]));
    U=V0*(tabP[i]-P0)/D;

    fprintf(fV, "%g \n",tabV[i] );
    fprintf(fE, "%g \n",tabE[i] );
    fprintf(fPh, "%g \n",tabP[i] );
    fprintf(fTh, "%g \n",tabT[i] );

    fprintf(fU, "%g \n",U );
    fprintf(fD, "%g \n",D );

    fprintf(fXA, "%g \n",0.0 );
    fprintf(fXB, "%g \n",0.0 );
    fprintf(fXC, "%g \n",1.0 );
  }
  
  
  fclose(fV);  fclose(fE);
  fclose(fPh);  fclose(fTh);
  fclose(fXA); fclose(fXB);  fclose(fXC);
  fclose(fU);  fclose(fD);  


  // PLOT PREMIERE HUGONIOT

  //****** Déclaration des fichiers pour plot octave  **********
  // Pour tester les valeurs
  //FILE *fPh, *fTh, *fXA, *fXB, *fXC, *fC;
  if((fPh = fopen("fichiers/Phug1.txt", "w+")) == NULL){printf("erreur ouverture fichier P test\n");return 1;}
  if((fTh = fopen("fichiers/Thug1.txt", "w+")) == NULL){printf("erreur ouverture fichier T test\n");return 1;}
  if((fXA = fopen("fichiers/XAhug1.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  if((fXB = fopen("fichiers/XBhug1.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  if((fXC = fopen("fichiers/XChug1.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  if((fC = fopen("fichiers/Chug1.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  //FILE *fV, *fE, *fU, *fD;
  if((fV = fopen("fichiers/Vhug1.txt", "w+")) == NULL){printf("erreur ouverture fichier V test \n");return 1;}
  if((fE = fopen("fichiers/Ehug1.txt", "w+")) == NULL){printf("erreur ouverture fichier E test \n");return 1;}
  if((fU = fopen("fichiers/Uhug1.txt", "w+")) == NULL){printf("erreur ouverture fichier U test \n");return 1;}
  if((fD = fopen("fichiers/Dhug1.txt", "w+")) == NULL){printf("erreur ouverture fichier D test \n");return 1;}

  // ********************************************************


  // DISCRETISATION DES PORTIONS DE COURBES ENTRE LES =! POINTS
  N=2e2;
  /*
  double* tabV=malloc(N*sizeof(double));  double* tabE=malloc(N*sizeof(double));
  double* tabP=malloc(N*sizeof(double));  double* tabT=malloc(N*sizeof(double));
  double* tabX=malloc(N*sizeof(double));  
  double* tabU=malloc(N*sizeof(double));  double* tabD=malloc(N*sizeof(double));  
  */
  /*
  // ZONE ENTRE 0 ET 1  (phase 1)
  
  Nmax=2e3; epsilon=1e-13; 
  phase=1;
  Vfin=90e-6;
  Vdeb=V0; Edeb=E0;
  err=courbe_hugoniot_pure_E(Nmax, epsilon, N, phase, Vdeb, Edeb, V0, E0, P0, Vfin, tabV, tabE, tabP, tabT);
  if(err){return 1;}

  for(int i=0; i<N; i++){
    D=V0*sqrt((tabP[i]-P0)/(V0-tabV[i]));
    U=V0*(tabP[i]-P0)/D;

    fprintf(fV, "%g \n",tabV[i] );
    fprintf(fE, "%g \n",tabE[i] );
    fprintf(fPh, "%g \n",tabP[i] );
    fprintf(fTh, "%g \n",tabT[i] );
    
    fprintf(fU, "%g \n",U );
    fprintf(fD, "%g \n",D );

    fprintf(fXA, "%g \n",1.0 );
    fprintf(fXB, "%g \n",0. );
    fprintf(fXC, "%g \n",0. ); 
  }
  
  
  fclose(fV);  fclose(fE);
  fclose(fPh);  fclose(fTh);
  fclose(fXA); fclose(fXB);  fclose(fXC);
  fclose(fU);  fclose(fD);  
  
  */
}