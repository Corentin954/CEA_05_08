#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "head_etain.h"


  

  // ZONE 
  // repere les zones (triangle pour le moment)
  //err=diag_zone();
  //if(err){return err;}
  
  // Frontière VE
  //err=routine_VE();
  //if(err){return err;}
  
  
  //**  MAILLAGE  **
  //err=maillage();
  
  
  // Zone mixte
  //mapsPT();

  // Ecriture de la tabulation
  
int main(){
  int err;
  
  // ZONE 
  // repere les zones (triangle pour le moment)
  //
  
  int diag_zone();
  err=diag_zone();
  if(err){return err;}  
  

  // Ecriture de la tabulation  
  /*
  int nb_points=80;

  int Nv_dis=5e2;
  int Ne_dis=5e2;
  double Vmin_dis=81.5e-6, Vmax_dis=215e-6;
  double Emin_dis=-7e4,  Emax_dis=6.7e6;
  
  int ecr_benchmark_PT(int nb_points, int Nv ,int Ne, double Vmin, double Vmax, double Emin, double Emax);
  err=ecr_benchmark_PT(nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis);
  if(err){return err;}
  */

  // USE BENCHMARK PT
  // Discrétisation du test
  /*
  int Nv_test=3e2;
  int Ne_test=3e2;

  double Vmin_test, Vmax_test;
  double Emin_test, Emax_test;

  Vmin_test=90e-6, Vmax_test=160e-6;
  Emin_test=-5e3,   Emax_test=10e5;
  
  int use_benchmark_PT(int Nv_test, int Ne_test, double Vmin_test, double Emin_test, double Vmax_test, double Emax_test);

  err=use_benchmark_PT(Nv_test, Ne_test, Vmin_test, Emin_test, Vmax_test, Emax_test);
  if(err){return err;}
  */

  // TEST POUR LA CIN2TIQUE DES PHASES
  // Discrétisation du test
  /*
  int Nv_test=200;
  int Ne_test=200;

  double Vmin_test, Vmax_test;
  double Emin_test, Emax_test;

  Vmin_test=110e-6, Vmax_test=145e-6;
  Emin_test=-5e3,   Emax_test=3e5;

  int use_benchmark_PT_cin_phase(int Nv_test, int Ne_test, double Vmin_test, double Emin_test, double Vmax_test, double Emax_test);

  err=use_benchmark_PT_cin_phase(Nv_test, Ne_test, Vmin_test, Emin_test, Vmax_test, Emax_test);
  if(err){return err;}
  */

  
  //int debug_Newton();
  //err=debug_Newton();


  // TEST FRACTION MASSIQUE DANS LE TRIANGLE TRIPLE
  /*
  int coorbary(double *xy, double **sommets, double **pcoorbar);
  // Valeur d'initialisation
  double VaT=129.905671875844291e-6;  double EaT=72.168190127265447e3;
  double VbT=127.666112185063994e-6;  double EbT=99.480256119236515e3;
  double VcT=132.538478805442338e-6;  double EcT=136.154346228414312e3;
  double* Xfm=malloc(3*sizeof(double));

  double* VE=malloc(2*sizeof(double));
  double V= 0.000129443; double E= 88522.7;
  VE[0]=V;   VE[1]=E;

  double** TRI=alloctabd(3,2);
  TRI[0][0]=VaT; TRI[0][1]=EaT;
  TRI[1][0]=VbT; TRI[1][1]=EbT;
  TRI[2][0]=VcT; TRI[2][1]=EcT;  
  err=coorbary(VE, TRI, &Xfm); if(err){return 1;} 
  
  printf(" V= %g, E= %g, xA= %f, xB= %g, xC= %g\n",VE[0],VE[1],Xfm[0],Xfm[1],Xfm[2] );
  */
  
  // Calcul d'une isentrope
  //int isentrope();
  //err=isentrope();


  //********   Calcul d'une courbe d'Hugoniot  *********
  int hugoniot();
  err=hugoniot();


  // **********  TEST DES DIFFERENTES METHODES  ******************** 
  // Discrétisation du test
  /*
  int Nv_test=3e2;
  int Ne_test=3e2;

  double Vmin_test, Vmax_test;
  double Emin_test, Emax_test;

  Vmin_test=98e-6, Vmax_test=150e-6;
  Emin_test=-5e3,   Emax_test=4e5;

  int methode=1;
  
  int use_fPTC_METH_cin_phase(int methode, int Nv_test, int Ne_test, double Vmin_test, double Emin_test, double Vmax_test, double Emax_test);

  err=use_fPTC_METH_cin_phase(methode, Nv_test, Ne_test, Vmin_test, Emin_test, Vmax_test, Emax_test);
  if(err){return err;}
  */

  
  


}