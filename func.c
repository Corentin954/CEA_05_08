#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "head.h"
#include "head_etain.h"

/* Retourne la valeur du pas de temps

*/
double pastemps(double T, double t, double dt, double dt_old){
  double res;
  if(dt>1.1*dt_old){
    dt=1.1*dt_old;
  }
  if ((T-t)<dt){ dt=T-t; }
  return dt;
}

// ind_u = 0 : vitesse U défini sur maille primal
// ind_u=1   : vitesse U défini sur maille duale
int print_sol(int ind_u, int ideb, int ifin, double* X, double* Xc, double* TAU, double* U, double* E, double* P, double* EPS){
  // ecriture dans un fichier des points (x,u(x)) de la solution
  FILE* Res;
  /*
  L1 : X  maillage
  L2 : tau  volume spécifique
  L3 : u  vitesse
  L4 : e   energie interne
  L5 : p   pression
  L6 : epsilon   energie interne spécifique
  */

  if((Res = fopen("sol.txt", "w+")) != NULL){
    // x
    for(int i=ideb; i<=ifin+1; i++){ 
      fprintf(Res,"%.15lf ",X[i]);
    }
    fprintf(Res, "\n\n");
    // X duale
    for (int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15lf ", (X[i+1]+X[i])/2);
    }
    fprintf(Res, "1e20\n");
    // tau
    for(int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15lf ",TAU[i]);
    }
    fprintf(Res, "1e20\n");
    fprintf(Res, "\n");
    // u
    for(int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15lf ",U[i]);
    }
    if(ind_u==0){
      fprintf(Res,"%.15lf ",U[ifin+1]);
      fprintf(Res, "\n\n");
    }
    else if(ind_u==1){
      fprintf(Res, "0\n\n");
    }
    // e
    for(int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15lf ",E[i]);
    }
    fprintf(Res, "1e20\n");
    fprintf(Res, "\n");
    // p
    for(int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15lf ",P[i]);
    }
    fprintf(Res, "1e20\n");
    fprintf(Res, "\n");
    // epsilon
    for(int i=ideb; i<=ifin; i++){ 
      fprintf(Res,"%.15lf ",EPS[i]);
      // fprintf(Res,"%.15lf ",REPS[i]/RHO0cAV[i]);
    }
    fprintf(Res, "1e20\n");
    fclose(Res);
    return 0;
  } 
  else{
    printf("Erreur d'ouverture de sol.txt (affiche_sol)\n");
    return 1;
  }
}

int print_nrj(int n, double* Etot, double* IMPUL, double* Mtot, double* VOL){
  FILE* Res;
  
  // Conservation au cours du temps
  if((Res = fopen("NRJ.txt", "w+")) != NULL){
    // energie totale
    for(int i=0; i<=n-1; i++){
      fprintf(Res,"%.15lf ",Etot[i]);
    }
    fprintf(Res, "\n");
    // impulsion
    for(int i=0; i<=n-1; i++){ 
      fprintf(Res,"%.15lf ",IMPUL[i]);
    }
    fprintf(Res, "\n");
    // masse totale
    for(int i=0; i<=n-1; i++){ 
      fprintf(Res,"%.15lf ",Mtot[i]);
    }
    fprintf(Res, "\n");
    // volume total
    for(int i=0; i<=n-1; i++){ 
      fprintf(Res,"%.15lf ",VOL[i]);
    }
    fprintf(Res, "\n");
    fclose(Res);
    return 0;
  } 
  else{
    printf("Erreur d'ouverture de NRJ.txt (affiche_nrj)\n");
    return 2;
  }
}

void print_err(int err, int Rsch, int sch){
  if(Rsch==0){
    printf("Godunov energie totale :\n");
  }
  else if(Rsch==1){
    printf("RUNGE-KUTTA energie interne :\n");
  }
  else if(Rsch==2){
    if(sch==0){
      printf("BBC JCP 2009 :\n");
    }
    else if(sch==1){
      printf("BBC Predictor-Corrector :\n");
    }
    else if(sch==2){
      printf("BBC RK2 average :\n");
    }
    else{printf("sch inconnu\n");}
    }
  else if(Rsch==3){
    printf("von Neumann-Ritchmyer :\n");
  }
  else if(Rsch==4){
    printf("Cauchy-Kovalevskaya :\n");
  }
  else{printf("sch inconnu\n");}
    

  if(err==1){
    printf("Erreur dans affichage_sol\n");
  }
  else if(err==2){
    printf("Erreur dans affichage_nrj\n");
  }
  else if(err==3){
    printf("Erreur dans u0\n");
  }
  else if(err==4){
    printf("Erreur dans initButcher\n");
  }
  else if(err==5){
    printf("Erreur dans condlim\n");
  }
  else{
    printf("err inconnu\n");
  }
}

// Met a jour les condition aux limites pour la variables Z
/* Si Z est centré aux maille  jfin= ifin
      Z est frontières aux mailles jfin=ifin+1

  int iu = 0 si cond sym
      iu = 1 si cond anti-sym
*/
int condlim(int ind_cond, int jdeb,int jfin, int nbghosts, double* Z, int ind){
  if(ind_cond==0){
    if (ind==0){
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=Z[jdeb+j];
        // à droite
        Z[jfin+1+j]=Z[jfin-j];
      }
    }
    else if (ind==1){  // U est frontière aux mailles, la valeur à la frontière n'intervient pas
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=-Z[jdeb+1+j];
        // à droite
        Z[jfin+1+j]=-Z[jfin-1-j];
      }
    }
    return 0;
  }
  // Condition périodique 
  else if(ind_cond==1){
    if(ind==0){
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=Z[jfin-j];
        // à droite
        Z[jfin+1+j]=Z[jdeb+j];
      }
    }
    else if(ind==1){ // pour l'energie (variable au carré)
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=Z[jfin+1+j];
        // à droite
        Z[jfin+1+j]=Z[jdeb-1-j];
      }
    }
    return 0;
  }
  // Condition aux limites pour l'etain : Mur à gauche et libre à droite
  else if(ind_cond==2){
    if (ind==0){
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=Z[jdeb+j];
        // à droite
        Z[jfin+1+j]=Z[jfin-j];
      }
    }
    else if (ind==1){  // U est frontière aux mailles, la valeur à la frontière n'intervient pas
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=-Z[jdeb+1+j];
        // à droite
        Z[jfin+1+j]=Z[jfin-1-j];
      }
    }
    return 0;
  }
  else{
    printf("Erreur choix de ind_cond (condlim)\n");
    return 5;
  }
}


// The Minmod limiter
double MinMod(double r){
  return fmax(0.,fmin(1.,r));
}

// The Superbee limiter
double Superbee(double r){
  return fmax(0., fmax(fmin(1.,2.*r),fmin(2.,r)) );
}

// The Christensen limiter
double Christensen(double r, double arg){
  return fmax(0., fmin(2.,fmin(2*r,arg)) );
}


double qvis(int q, double Cq, double Cl, double* U, double tau, double c, double* DM){
  int third=q%10;               // pseudo-viscosity
  int second=(q%100-third)/10;  // limiteur pour l'ordre 2
  int first=q/100;              // ordre


  double du=U[1]-U[0];
  // ordre 1
  //   first==0
  // Ordre 2
  if (first==1){
    // delta vitesse
    double dUm=U[0]-U[-1];
    double dUp=U[2]-U[1];
    // delta masse
    double dM=DM[0];
    double dMm=DM[-1];
    double dMp=DM[1];

    double rplus=(dUp/dMp)*(dM/du);
    double rmoins=(dUm/dMm)*(dM/du);
    
    double phiplus;
    double phimoins;

    if(second==0){
      phiplus=MinMod(rplus);
      phimoins=MinMod(rmoins);
    }
    else if(second==1){
     phiplus=Superbee(rplus);
     phimoins=Superbee(rmoins);
    }
    else if(second==2){
      double argplus=(dM*rplus+dMp)/(dM+dMp);
      double argmoins=(dM*rmoins+dMm)/(dM+dMm);
      phiplus=Christensen(rplus,argplus);
      phimoins=Christensen(rmoins, argmoins);
    }
    else{ printf("Erreur q= %d (qvis)\n",q );  }

    du*=(1-0.5*(phiplus + phimoins));
  }

  
  // von Neumann-Ritchmyer
  if (third == 0){
    return  -Cq*(du)*fabs(du)/tau;
  }
  // Rosenbluth
  else if(third == 1){
  	if (du<0){
  	  return -Cq*(du)*fabs(du)/tau;
  	}
  	else{ return 0.;}
  }
  // Landschoff
  else if(third == 2){
  	if (du<0){
  	  return -Cq*(du)*fabs(du)/tau -Cl*c*du/tau;
  	}
  	else{ return 0;}
  }
  // "Magical" q combination
  else if(third == 3){
  	if (du<0){
  	  return -Cq*(du)*fabs(du)/tau -Cl*c*du/tau;
  	}
  	else{ return -Cl*c*du/tau;}
  }
  // Landschoff
  else if(third == 4){
    return -Cq*(du)*fabs(du)/tau -Cl*c*du/tau;
  }
  else if(third == 5){
    return 0.;
  }
  else{
    printf("Erreur dans le choix de pseudo q= %d (qvis)\n",q);
  }
}

// Evaluation de la moyenne (Average) ou ponctuelle (valable jusqu'a l'ordre 4) de la Z
// Le type d'evaluation depend des arguments de la focntion
// k : spatial order
// si . Coef=Ck retourne Point-Wise
//    . Coef=Ckbar retourne Average
// Z = &Z[i]
double phi(int so, double* Z, double** Coef){
  double res=0.0;
  int k=so-2;
  for(int j=-k; j<=k; j++){
    res+=Coef[k][abs(j)]*Z[j];
  }
  return res;
}


// Même principe que la fonction au dessus (renvoie la moyenne ou le valeur ponctuelle)
// sauf qu'elle prend en arguments les tableaux : rho.Z et rho 
// il faut faire leur division pour obtenir Z
double phidiv(int so, double* RZ, double* RHO, double** Coef){
  double res=0.0;
  int k=so-2;
  for(int j=-k; j<=k; j++){
    res+=Coef[k][abs(j)]*(RZ[j]/RHO[j]);
  }
  return res;
}

// Calcul la maille dual (centrée) à partir de la maille primal (frontière) 
double fxc(int so,double* X, double** rk){
  // voir CRAS 2016 p214
  double res=0.0;
  int k=so-2;
  for(int l=0; l<=k; l++){
    res+=rk[k][l]*(X[l+1] + X[-l]);
  }
  return res;
}

// Pour les variables centrées aux mailles :   ind=0
//   delta i+1/2 = Z[i]-Z[i-1]
// Pour les variables frontières aux mailles :   ind=1
//   delta i = Z[i+1]-Z[i]
double delta(int so,int ind, double* Z, double** dk){
  double dz=0.0;
  int k=so-2;
  if (ind==0){
    for(int l=0; l<=k; l++){
      //dz+=dk[k][l]*(Z[l+1] - Z[l]);
      dz+=dk[k][l]*(Z[l] - Z[-l-1]);
    }
  }
  else if (ind==1){
    for(int l=0; l<=k; l++){
      //dz+=dk[k][l]*(Z[l+1] - Z[l]);
      dz+=dk[k][l]*(Z[l+1] - Z[-l]);
    }
  }
  return dz;
}


// Pour rho.eps  , retourne AV(PI.\delta u)
// Retourne la moyenne de YdZ  
// soit : . Y=PI et Z=u  pour EPS
//        . Y=u  et Z=PI pour EKIN
// int    cas==0 : EPS    cas==1 : EKIN
double YdZ(int so, int cas, double* P, double* Q, double* u, double** dk, double** Ckbar){
  double du, res=0.0;
  int k=so-2;
  if(cas==0){ // EPS+=(dt/dx).PIdu
    for(int j=-k; j<=k; j++){
      // calcul de du
      du=0.0;
      for(int l=0; l<=k; l++){
        du+=dk[k][l]*(u[j+(l+1)]-u[j-l]);
      }
      // PI.du
      res+=Ckbar[k][abs(j)]*(P[j]+Q[j])*du;
    }
  }
  else if(cas==1){ // EKIN+=(dt/dx).udPI
    for(int j=-k; j<=k; j++){
      // calcul de dPI
      du=0.0;
      for(int l=0; l<=k; l++){
        du+=dk[k][l]*( (P[j+l]+Q[j+l])-(P[j-l-1]+Q[j-l-1]) );
      }
      // u.dPI
      res+=Ckbar[k][abs(j)]*u[j]*du;
    }
  }
  return res;
}

// Fonction uniquement pour le kinetic energy fix  
// sum Qbar.Delta K
// k : spatial order
// Z = &Z[i]
double phiQK(int so, double* DK, double** Q){
  double res=0.0;
  int k=so-2;
  for(int j=0; j<=k; j++){
    res+=Q[k][j]*(DK[j+1]+DK[-j]);
  }
  return res;
}


void initCoef(double** Ck, double** Ckbar, double** dk, double** Qbar, double** rk){
  // Ck
  Ck[0][0]=1.0;               Ck[0][1]=0.0;               Ck[0][2]=0.0;             Ck[0][3]=0.0;           Ck[0][4]=0.0; 
  Ck[1][0]=13./12;            Ck[1][1]=-1./24;            Ck[1][2]=0.0;             Ck[1][3]=0.0;           Ck[1][4]=0.0; 
  Ck[2][0]=1067./960;         Ck[2][1]=-29./480;          Ck[2][2]=3./640.;         Ck[2][3]=0.0;           Ck[2][4]=0.0; 
  Ck[3][0]=30251./26880;      Ck[3][1]=-7621./107520;     Ck[3][2]=159./17920;      Ck[3][3]=-5./7168;      Ck[3][4]=0.0; 
  Ck[4][0]=5851067./5160960;  Ck[4][1]=-100027./1290240;  Ck[4][2]=31471./2580480;  Ck[4][3]=-425./258048;  Ck[4][4]=35./294912; 
  // Ckbar
  Ckbar[0][0]=1.0;                 Ckbar[0][1]=0.0;                Ckbar[0][2]=0.0;                 Ckbar[0][3]=0.0;              Ckbar[0][4]=0.0; 
  Ckbar[1][0]=11./12;              Ckbar[1][1]=1./24;              Ckbar[1][2]=0.0;                 Ckbar[1][3]=0.0;              Ckbar[1][4]=0.0; 
  Ckbar[2][0]=863./960;            Ckbar[2][1]=77./1440;           Ckbar[2][2]=-17./5760;           Ckbar[2][3]=0.0;              Ckbar[2][4]=0.0;
  Ckbar[3][0]=215641./241920;      Ckbar[3][1]=6361./107520;       Ckbar[3][2]=-281./53760;         Ckbar[3][3]=367./967680;      Ckbar[3][4]=0.0; 
  Ckbar[4][0]=41208059./46448640;  Ckbar[4][1]=3629953./58060800;  Ckbar[4][2]=-801973./116121600;  Ckbar[4][3]=49879./58060800;  Ckbar[4][4]=-27859./464486400;  
  // dk
  dk[0][0]=1.0;           dk[0][1]=0.0;         dk[0][2]=0.0;         dk[0][3]=0.0;           dk[0][4]=0.0; 
  dk[1][0]=9./8;          dk[1][1]=-1./24;      dk[1][2]=0.0;         dk[1][3]=0.0;           dk[1][4]=0.0; 
  dk[2][0]=75./64;        dk[2][1]=-25./384;    dk[2][2]=3./640;      dk[2][3]=0.0;           dk[2][4]=0.0; 
  dk[3][0]=1225./1024;    dk[3][1]=-245./3072;  dk[3][2]=49./5120;    dk[3][3]=-5./7168;      dk[3][4]=0.0; 
  dk[4][0]=19845./16384;  dk[4][1]=-735./8192;  dk[4][2]=567./40960;  dk[4][3]=-405./229376;  dk[4][4]=35./294912; 
  // Qbar
  Qbar[0][0]=1./2;              Qbar[0][1]=0.0;             Qbar[0][2]=0.0;           Qbar[0][3]=0.0;              Qbar[0][4]=0.0; 
  Qbar[1][0]=13./24;            Qbar[1][1]=-1./24;          Qbar[1][2]=0.0;           Qbar[1][3]=0.0;              Qbar[1][4]=0.0; 
  Qbar[2][0]=401./720;          Qbar[2][1]=-31./480;        Qbar[2][2]=11./1440;      Qbar[2][3]=0.0;              Qbar[2][4]=0.0; 
  Qbar[3][0]=68323./120960;     Qbar[3][1]=-353./4480;      Qbar[3][2]=1879./120960;  Qbar[3][3]=-191./120960;     Qbar[3][4]=0.0; 
  Qbar[4][0]=2067169./3628800;  Qbar[4][1]=-40111./453600;  Qbar[4][2]=581./25920;    Qbar[4][3]=-28939./7257600;  Qbar[4][4]=2497./7257600; 
  // rk
  rk[0][0]=1./2;          rk[0][1]=0.0;           rk[0][2]=0.0;         rk[0][3]=0.0;          rk[0][4]=0.0; 
  rk[1][0]=9./16;         rk[1][1]=-1./16;        rk[1][2]=0.0;         rk[1][3]=0.0;          rk[1][4]=0.0; 
  rk[2][0]=75./128;       rk[2][1]=-25./256;      rk[2][2]=3./256.;     rk[2][3]=0.0;          rk[2][4]=0.0; 
  rk[3][0]=1225./2048;    rk[3][1]=-245./2048;    rk[3][2]=49./2048;    rk[3][3]=-5./2048;     rk[3][4]=0.0; 
  rk[4][0]=19845./32768;  rk[4][1]=-2205./16384;  rk[4][2]=567./16384;  rk[4][3]=-405./65536;  rk[4][4]=35./65536; 
}


int u0(int tst, double a, double b,double x, double* W){ 
  double eps=1e-8; 
  double k=2*M_PI; //////////////////////////
  double rho0=1.0; double p0=5./7;
  double Gamma;
  if (tst==0){ // Sod
    if(x<=a+(b-a)/2.0){
      W[0]=1.0; // rho
      W[1]=0.0; // u
      W[2]=1.0; // p
    }
    else{
      W[0]=0.125; // rho
      W[1]=0.0; // u
      W[2]=0.1; // p
      //W[2]=1.0;
    }
    return 0;
  }
  else if (tst==2){ // Bizarrium
    if(x<=a+(b-a)/2.0){
      W[0]=1./0.7e-4; // rho
      W[1]=0.0;       // u
      W[2]=1e11;      // p
    }
    else{
      W[0]=1e4; // rho
      W[1]=250.0; // u
      W[2]=0.0; // p
    }
    return 0;
  }
  else if (tst==3){ // Onde acoustique
    W[0]=rho0 + eps*sin(k*x);  // rho
    W[1]=eps*sin(k*x);  // u
    W[2]=p0 + eps*sin(k*x);  // p

    return 0;
  }
  else if (tst==4){ // Sod symétrisé
    if(x<=a+2.0*(b-a)/5){
      W[0]=0.125; // rho
      W[1]=0.0;       // u
      W[2]=0.1;      // p
    }
    else if(a+3.0*(b-a)/5<=x){
      W[0]=0.125; // rho
      W[1]=0.0;  // u
      W[2]=0.1; // p
    }
    else{  // if(2.0*(a+b)/5<x || x<3.0*(a+b)/5){
      W[0]=1.0; // rho
      W[1]=0.0; // u
      W[2]=1.0; // p
    }
    return 0;
  }
  else if (tst==1){ // LeBlanc
    Gamma=5./3;
    if(x<=a+(b-a)/3.0){
      W[0]=1.0; // rho
      W[1]=0.0; // u
      W[2]=(Gamma-1)/10; // p
    }
    else{
      W[0]=1e-3; // rho
      W[1]=0.0; // u
      W[2]=(Gamma-1)*1e-3*1e-9; // p
      //W[2]=1.0;
    }
    return 0;
  }
  else if (tst==5){ // Woodward 3 etats 
    if(x<=a+(b-a)/10){
      W[0]=1.0; // rho
      W[1]=0.0;       // u
      W[2]=1e3;      // p
    }
    else if(a+9.0*(b-a)/10<=x){
      W[0]=1.0; // rho
      W[1]=0.0;  // u
      W[2]=1e2; // p
    }
    else{  // if(2.0*(a+b)/5<x || x<3.0*(a+b)/5){
      W[0]=1.0; // rho
      W[1]=0.0; // u
      W[2]=1e-2; // p
    }
    return 0;
  }
  else if (tst==6){ // Shu_Oscher 
    if(x<=a+(b-a)/10){
      W[0]=27./7;        // rho
      W[1]=4*sqrt(35)/9;        // u
      W[2]=31./3; // p
    }
    else{  
      W[0]=1.0+sin(5.*x)/5; // rho
      W[1]=0; // u
      W[2]=1.0; // p
    }
    return 0;
  }
  else{ 
    printf("Erreur dans le choix de tst (u0)\n");
    return 3;
  }
}


int u0_multi(int tst, int ind, int ind_air, int ind_etain, double a, double b,double x, double* W){ 
  double eps=1e-8; 
  double k=2*M_PI; //////////////////////////
  double rho0=1.0; double p0=5./7;
  double Gamma;
  if (tst==0){ // Sod
    if(x<=a+(b-a)/2.0){
      W[0]=1.0; // rho
      W[1]=0.0; // u
      W[2]=1.0; // p
    }
    else{
      W[0]=0.125; // rho
      W[1]=0.0; // u
      W[2]=0.1; // p
      //W[2]=1.0;
    }
    return 0;
  }
  else if (tst==2){ // Bizarrium
    if(x<=a+(b-a)/2.0){
      W[0]=1./0.7e-4; // rho
      W[1]=0.0;       // u
      W[2]=1e11;      // p
    }
    else{
      W[0]=1e4; // rho
      W[1]=250.0; // u
      W[2]=0.0; // p
    }
    return 0;
  }
  else if (tst==3){ // Onde acoustique
    W[0]=rho0 + eps*sin(k*x);  // rho
    W[1]=eps*sin(k*x);  // u
    W[2]=p0 + eps*sin(k*x);  // p

    return 0;
  }
  else if (tst==4){ // Sod symétrisé
    if(x<=a+2.0*(b-a)/5){
      W[0]=0.125; // rho
      W[1]=0.0;       // u
      W[2]=0.1;      // p
    }
    else if(a+3.0*(b-a)/5<=x){
      W[0]=0.125; // rho
      W[1]=0.0;  // u
      W[2]=0.1; // p
    }
    else{  // if(2.0*(a+b)/5<x || x<3.0*(a+b)/5){
      W[0]=1.0; // rho
      W[1]=0.0; // u
      W[2]=1.0; // p
    }
    return 0;
  }
  else if (tst==1){ // LeBlanc
    Gamma=5./3;
    if(x<=a+(b-a)/3.0){
      W[0]=1.0; // rho
      W[1]=0.0; // u
      W[2]=(Gamma-1)/10; // p
    }
    else{
      W[0]=1e-3; // rho
      W[1]=0.0; // u
      W[2]=(Gamma-1)*1e-3*1e-9; // p
      //W[2]=1.0;
    }
    return 0;
  }
  else if (tst==5){ // Woodward 3 etats 
    if(x<=a+(b-a)/10){
      W[0]=1.0; // rho
      W[1]=0.0;       // u
      W[2]=1e3;      // p
    }
    else if(a+9.0*(b-a)/10<=x){
      W[0]=1.0; // rho
      W[1]=0.0;  // u
      W[2]=1e2; // p
    }
    else{  // if(2.0*(a+b)/5<x || x<3.0*(a+b)/5){
      W[0]=1.0; // rho
      W[1]=0.0; // u
      W[2]=1e-2; // p
    }
    return 0;
  }
  else if (tst==6){ // Shu_Oscher 
    if(x<=a+(b-a)/10){
      W[0]=27./7;        // rho
      W[1]=4*sqrt(35)/9;        // u
      W[2]=31./3; // p
    }
    else{  
      W[0]=1.0+sin(5.*x)/5; // rho
      W[1]=0; // u
      W[2]=1.0; // p
    }
    return 0;
  }
  double u;
  if(tst==100 ){ u=320.756; }
  else if(tst==101){ u=320.756; } // x=0    beta/gamma      POLE 1
  else if(tst==102){ u=331.046; } // x=0.2  beta/gamma
  else if(tst==103){ u=346.306; } // x=0.5  beta/gamma
  else if(tst==104){ u=361.359; } // x=0.8  beta/gamma
  else if(tst==105){ u=371.266; } // x=1    beta/gamma      POLE 2
  else if(tst==106){ u=474.570; } // phase gamma            POLE 3
  else if(tst==107){ u=1333.77; } // x=0    gamma/liquide   POLE 4
  else if(tst==108){ u=1383.13; } // x=0.2  gamma/liquide
  else if(tst==109){ u=1453.67; } // x=0.5  gamma/liquide
  else if(tst==110){ u=1520.7;  } // x=0.8  gamma/liquide
  else if(tst==111){ u=1563.72; } // x=1    gamma/liquide   POLE 5
  else if(tst==112){ u=900; }     // x=1    gamma/liquide   POLE 5
  else if(tst==113){ u=1700; }    // x=1    gamma/liquide   POLE 5

  if (tst==100){ // CHOC SIMPLE
    if(ind==ind_air){  W[0]=1.27;} // rho
    else if(ind==ind_etain){  W[0]=7287.;} // rho
    W[1]=-u;  //-1100; //-400.761; // u
    W[2]=1e5; // P
    W[3]=300.0; // T
    W[4]=0.;  // C
    W[5]=1.;  // xA
    W[6]=0.0; // xB
    W[7]=0.0; // xC
    return 0;
  } 
  else if(200>tst && tst>100) { // CHOC SYMETRISE
    //printf(" choc sym\n");
    if(x<=a+(b-a)/2.0){
      if(ind==ind_air){  W[0]=1.27;} // rho
      else if(ind==ind_etain){ W[0]=7287.;} // rho
      W[1]=+u;  //-1100; //-400.761; // u
      W[2]=1e5; // P
      W[3]=300.0; // T
      W[4]=0.;  // C
      W[5]=1.;  // xA
      W[6]=0.0; // xB
      W[7]=0.0; // xC
    }
    else{
      if(ind==ind_air){  W[0]=1.27;} // rho
      else if(ind==ind_etain){  W[0]=7287.;} // rho
      W[1]=-u;  //-1100; //-400.761; // u
      W[2]=1e5; // P
      W[3]=300.0; // T
      W[4]=0.;  // C
      W[5]=1.;  // xA
      W[6]=0.0; // xB
      W[7]=0.0; // xC
    }
    return 0;
  }    
  else if(tst==200) { // CHOC SYMETRISE    ETAIN DANS L'ETAT P4
    //printf(" choc sym\n");
    if(ind==ind_air){  
    W[0]=1.27;
    W[1]=0;  //-1100; //-400.761; // u
    W[2]=1e5; // P
    W[3]=1965.674770; // T
    W[4]= 3897.63;  // C
    W[5]=0.;  // xA
    W[6]=1.0; // xB
    W[7]=0.0; // xC

    } // rho
    else if(ind==ind_etain){  
      W[0]=10247.863117; // rho
      W[1]=0;  //-1100; //-400.761; // u
      W[2]=44.866639e9; // P
      W[3]=1965.674770; // T
      W[4]= 3897.63;  // C
      W[5]=0.;  // xA
      W[6]=1.0; // xB
      W[7]=0.0; // xC
    }
    /*
    else{
      if(ind==ind_air){  W[0]=1.27;} // rho
      else if(ind==ind_etain){  W[0]=1./10247.863117;} // rho
      W[1]=0;  //-1100; //-400.761; // u
      W[2]=44.866639e9; // P
      W[3]=1965.674770; // T
    }
    */
    return 0;
  }
  

  printf("Erreur dans le choix de tst (u0)\n");
  return 3;

}


int u0_etain(int tst, double a, double b,double x, double* W){ 
  double u;
  if(tst==100 ){ u=320.756; }
  else if(tst==101){ u=320.756; } // x=0    beta/gamma      POLE 1
  else if(tst==102){ u=331.046; } // x=0.2  beta/gamma
  else if(tst==103){ u=346.306; } // x=0.5  beta/gamma
  else if(tst==104){ u=361.359; } // x=0.8  beta/gamma
  else if(tst==105){ u=371.266; } // x=1    beta/gamma      POLE 2
  else if(tst==106){ u=474.570; } // phase gamma            POLE 3
  else if(tst==107){ u=1333.77; } // x=0    gamma/liquide   POLE 4
  else if(tst==108){ u=1383.13; } // x=0.2  gamma/liquide
  else if(tst==109){ u=1453.67; } // x=0.5  gamma/liquide
  else if(tst==110){ u=1520.7;  } // x=0.8  gamma/liquide
  else if(tst==111){ u=1563.72; } // x=1    gamma/liquide   POLE 5
  else if(tst==112){ u=900; } // x=1    gamma/liquide   POLE 5
  else if(tst==113){ u=1700; } // x=1    gamma/liquide   POLE 5
  
  if (tst==100){ // CHOC SIMPLE
    W[0]=7287.; // rho
    W[1]=-u;  //-1100; //-400.761; // u
    W[2]=101325; // P
    W[3]=300.0; // T
    return 0;
  } 
  else if(tst>100) { // CHOC SYMETRISE
    //printf(" choc sym\n");
    if(x<=a+(b-a)/2.0){
      W[0]=7287.; // rho
      W[1]=+u;  //-1100; //-400.761; // u
      W[2]=101325; // P
      W[3]=300.0; // T
    }
    else{
      W[0]=7287.; // rho
      W[1]=-u;  //-1100; //-400.761; // u
      W[2]=101325; // P
      W[3]=300.0; // T
    }
    return 0;
  }
}
  /*
  POLE
  rho0= 7287.000021, V0= 137.230684, E0= -0.009371 , P0= 0.000000, T0= 300.000000
  rho1= 8076.006891, V1= 123.823569, E1= 51443.000837 , P1= 7.673987, T1= 404.147829 , D1= 3283.175368, u1= 320.758508
  rho2= 8288.370443, V2= 120.650978, E2= 77064.926756 , P2= 8.478062, T2= 354.880625 , D2= 3283.175368, u2= 371.265969
  rho3= 8518.284342, V3= 117.394532, E3= 112608.341615 , P3= 11.353850, T3= 388.435466 , D3= 3283.175368, u3= 474.570018
  rho4= 10247.863117, V4= 97.581319, E4= 889466.874467 , P4= 44.866639, T4= 1965.674770 , D4= 4616.309074, u4= 1333.766759
  rho4= 10609.008403, V5= 94.259516, E5= 1222604.764353 , P5= 56.903493, T5= 2231.310798 , D5= 4993.811106, u5= 1563.716582


  ZONE MIXTE beta/gamma (1/2)
  x= 0, P= 7.67399e+09, u= 320.759, rho= 8076.01, V= 0.000123824, E= 51443, T= 404.148
  x= 0.1, P= 7.75565e+09, u= 325.913, rho= 8097.28, V= 0.000123498, E= 53952.9, T= 399.069
  x= 0.2, P= 7.83711e+09, u= 331.046, rho= 8118.55, V= 0.000123175, E= 56474.8, T= 394.022
  x= 0.3, P= 7.91833e+09, u= 336.157, rho= 8139.8, V= 0.000122853, E= 59008.6, T= 389.007
  x= 0.4, P= 7.99928e+09, u= 341.245, rho= 8161.04, V= 0.000122533, E= 61554.1, T= 384.026
  x= 0.5, P= 8.07994e+09, u= 346.31, rho= 8182.28, V= 0.000122215, E= 64111.2, T= 379.079
  x= 0.6, P= 8.16029e+09, u= 351.351, rho= 8203.51, V= 0.000121899, E= 66679.8, T= 374.167
  x= 0.7, P= 8.2403e+09, u= 356.368, rho= 8224.73, V= 0.000121584, E= 69259.6, T= 369.29
  x= 0.8, P= 8.31994e+09, u= 361.359, rho= 8245.95, V= 0.000121272, E= 71850.5, T= 364.45
  x= 0.9, P= 8.3992e+09, u= 366.325, rho= 8267.16, V= 0.00012096, E= 74452.3, T= 359.646
  x= 1, P= 8.47806e+09, u= 371.266, rho= 8288.37, V= 0.000120651, E= 77064.9, T= 354.881


 ZONE MIXTE gamma/liquide (2/3)
 x= 0, P= 4.48666e+10, u= 1333.77, rho= 10247.9, V= 9.75813e-05, E= 889467, T= 1965.67
 x= 0.1, P= 4.61233e+10, u= 1358.7, rho= 10287.4, V= 9.72059e-05, E= 923037, T= 1995.08
 x= 0.2, P= 4.73658e+10, u= 1383.13, rho= 10326.1, V= 9.6842e-05, E= 956520, T= 2023.74
 x= 0.3, P= 4.85952e+10, u= 1407.07, rho= 10363.9, V= 9.64888e-05, E= 989930, T= 2051.71
 x= 0.4, P= 4.98125e+10, u= 1430.58, rho= 10400.9, V= 9.61455e-05, E= 1.02328e+06, T= 2079.02
 x= 0.5, P= 5.10185e+10, u= 1453.67, rho= 10437.2, V= 9.58113e-05, E= 1.05658e+06, T= 2105.73
 x= 0.6, P= 5.22141e+10, u= 1476.37, rho= 10472.8, V= 9.54858e-05, E= 1.08983e+06, T= 2131.87
 x= 0.7, P= 5.33997e+10, u= 1498.71, rho= 10507.7, V= 9.51683e-05, E= 1.12306e+06, T= 2157.47
 x= 0.8, P= 5.45762e+10, u= 1520.7, rho= 10542, V= 9.48584e-05, E= 1.15626e+06, T= 2182.56
 x= 0.9, P= 5.57439e+10, u= 1542.36, rho= 10575.8, V= 9.45556e-05, E= 1.18944e+06, T= 2207.16
 x= 1, P= 5.69035e+10, u= 1563.72, rho= 10609, V= 9.42595e-05, E= 1.2226e+06, T= 2231.31
 */




int initButcher(int sch, double** A, double* THETA, double* ALPHA){
  // RK1  (Matsuno forward-backward)
  if (sch==0){
    A[0][0]=1.0;

    THETA[0]=0.0; THETA[1]=1.0;

    ALPHA[0]=1.;
    return 0;
  }
  // RK2 (SSP) HEUN
  else if (sch==1){
    A[0][0]=1.0;

    THETA[0]=1.0/2.0; THETA[1]=1.0/2.0;

    ALPHA[0]=1.;
    return 0;
  }
  // RK3 (SSP)
  else if (sch==2){
    A[0][0]=1.0; 
    A[1][0]=1.0/4.0; A[1][1]=1.0/4.0;

    THETA[0]=1.0/6.0; THETA[1]=1.0/6.0; THETA[2]=2.0/3.0;

    ALPHA[0]=1.; ALPHA[1]=1./2;
    return 0;
  }
  // KUTTA ORDRE 4 (Classic RK4)
  else if (sch==3){
    A[0][0]=1.0/2.0; 
    A[1][0]=0.0;     A[1][1]=1.0/2.0;
    A[2][0]=0.0;     A[2][1]=0.0;      A[2][2]=1.0;

    THETA[0]=1.0/6.0; THETA[1]=1.0/3.0; THETA[2]=1.0/3.0; THETA[3]=1.0/6.0;

    ALPHA[0]=1./2; ALPHA[1]=1./2; ALPHA[2]=1.;
    return 0;
  }
  // Cash-Karp (ordre 5)
  else if (sch==4){
    A[0][0]=1.0/5.0; 
    A[1][0]=3./40;        A[1][1]=9./40;
    A[2][0]=3./10;        A[2][1]=-9./10;    A[2][2]=6./5;
    A[3][0]=-11./54;      A[3][1]=5./2;      A[3][2]=-70./27;     A[3][3]=-35./27;
    A[4][0]=1631./55296;  A[4][1]=175./512;  A[4][2]=575./13824;  A[4][3]=44275./110592;  A[4][4]=253./4096;

    THETA[0]=37./378; THETA[1]=0.; THETA[2]=250./621; THETA[3]=125./594; THETA[4]=0.; THETA[5]=512./1771;  

    ALPHA[0]=1./5; ALPHA[1]=3./10; ALPHA[2]=3./5; ALPHA[3]=1.; ALPHA[4]=7./8;
    return 0;
  }
  // Dormand-Prince
  else if (sch==5){
    A[0][0]=1.0/5.0; 
    A[1][0]=3./40;         A[1][1]=9./40;
    A[2][0]=44./45;        A[2][1]=-56./15;        A[2][2]=32./9;
    A[3][0]=19372./6561;   A[3][1]=-25360./2187;   A[3][2]=64448./6561;   A[3][3]=-212./729;
    A[4][0]=9017./3168;    A[4][1]=-355./33;       A[4][2]=46732./5247;  A[4][3]=49./176;    A[4][4]=-5103./18656;
    A[5][0]=35./384;       A[5][1]=0.;             A[5][2]=500./1113;    A[5][3]=125./192;   A[5][4]=-2187./6784;  A[5][5]=11./84;

    THETA[0]=35./384; THETA[1]=0.; THETA[2]=500./1113; THETA[3]=125./192; THETA[4]=-2187./6784; THETA[5]=11./84; THETA[6]=0.;

    ALPHA[0]=1./5; ALPHA[1]=3./10; ALPHA[2]=4./5; ALPHA[3]=8./9; ALPHA[4]=1.; ALPHA[5]=1.;
    return 0;
  }
  // HEUN  ordre 3
  else if (sch==6){
    A[0][0]=1.0/3; 
    A[1][0]=0; A[1][1]=2.0/3;

    THETA[0]=1.0/4; THETA[1]=0; THETA[2]=3./4;

    ALPHA[0]=1./3; ALPHA[1]=2./3;
    return 0;
  }
  // RALSTON  ordre 3
  else if (sch==7){
    A[0][0]=1.0/2; 
    A[1][0]=0; A[1][1]=3.0/4;

    THETA[0]=2.0/9; THETA[1]=1./3; THETA[2]=4./9;

    ALPHA[0]=1./2; ALPHA[1]=3./4;
    return 0;
  }
  // Bogacki-Shampine ordre 3 
  else if (sch==8){
    A[0][0]=1.0/2.0; 
    A[1][0]=0.0;     A[1][1]=3.0/4;
    A[2][0]=2./9;    A[2][1]=1./3;      A[2][2]=4./9;

    THETA[0]=2./9; THETA[1]=1.0/3; THETA[2]=4./9; THETA[3]=0;

    ALPHA[0]=1./2; ALPHA[1]=3./4; ALPHA[2]=1.; 
    return 0;
  }
  // Bogacki-Shampine ordre 2 
  else if (sch==9){
    A[0][0]=1.0/2.0; 
    A[1][0]=0.0;     A[1][1]=3.0/4;
    A[2][0]=2./9;    A[2][1]=1./3;      A[2][2]=4./9;

    THETA[0]=7./24; THETA[1]=1.0/4; THETA[2]=1./3; THETA[3]=1./8;

    ALPHA[0]=1./2; ALPHA[1]=3./4; ALPHA[2]=1.; 
    return 0;
  }
  // KUTTA odre 3
  else if (sch==10){
    A[0][0]=1./2; 
    A[1][0]=-1.; A[1][1]=2.0;

    THETA[0]=1./6; THETA[1]=2./3; THETA[2]=1./6;

    ALPHA[0]=1./2; ALPHA[1]=1.;
    return 0;
  }
  // KUTTA ordre 2 (Explicit Midpoint Method) 
  else if (sch==11){
    A[0][0]=1.0/2;

    THETA[0]=0.0; THETA[1]=1.0;

    ALPHA[0]=1./2;
    return 0;
  }
  // KUTTA ordre 2 Ralston
  else if (sch==12){
    A[0][0]=2.0/3;

    THETA[0]=1.0/4; THETA[1]=3./4;

    ALPHA[0]=2./3;
    return 0;
  }
  // EULER forward
  else if (sch==13){
    A[0][0]=0;

    THETA[0]=0; THETA[1]=1.0;

    ALPHA[0]=0.0;
    return 0;
  }
  // SSPRK4(5,4) Spiteri-Ruuth
  else if (sch==14){
    A[0][0]=0.39175222700392; 
    A[1][0]=0.21766909633821;  A[1][1]=0.36841059262959;
    A[2][0]=0.08269208670950;  A[2][1]=0.13995850206999;  A[2][2]=0.25189177424738;
    A[3][0]=0.06796628370320;  A[3][1]=0.11503469844438;  A[3][2]=0.20703489864929;  A[3][3]=0.54497475021237;

    THETA[0]=0.14681187618661; THETA[1]=0.24848290924556; THETA[2]=0.10425883036650; THETA[3]=0.27443890091960; THETA[4]=0.22600748319395; 

    ALPHA[0]=0.39175222700392; ALPHA[1]=0.58607968896779; ALPHA[2]=0.47454236302687; ALPHA[3]=0.93501063100924; 
    return 0;
  }
  // SSPRK3(4,3) Spiteri-Ruuth
  else if (sch==15){
    A[0][0]=0.5; 
    A[1][0]=0.5;   A[1][1]=0.5;
    A[2][0]=1./6;  A[2][1]=1./6;  A[2][2]=1./6;

    THETA[0]=1./6; THETA[1]=1./6; THETA[2]=1./6; THETA[3]=0.5; 

    ALPHA[0]=0.5; ALPHA[1]=1.0; ALPHA[2]=0.5; 
    return 0;
  }
  // SSPRK3(5,3) Spiteri-Ruuth
  else if (sch==16){
    A[0][0]=0.37726891511710; 
    A[1][0]=0.75453783023419;  A[1][1]=0.37726891511710;
    A[2][0]=0.49056882269314;  A[2][1]=0.16352294089771;  A[2][2]=0.16352294089771;
    A[3][0]=0.78784303014311;  A[3][1]=0.14831273384724;  A[3][2]=0.14831273384724;  A[3][3]=0.34217696850008;

    THETA[0]=0.19707596384481; THETA[1]=0.11780316509765; THETA[2]=0.11709725193772; THETA[3]=0.27015874934251; THETA[4]=0.29786487010104; 

    ALPHA[0]=0.37726891511710; ALPHA[1]=0.75453783023419; ALPHA[2]=0.49056882269314; ALPHA[3]=0.78784303014311; 
    return 0;
  }
  // SPPRK1(3,1) Spiteri-Ruuth
  else if (sch==17){
    A[0][0]=1./3; 
    A[1][0]=1./3; A[1][1]=1./3;

    THETA[0]=1./3; THETA[1]=1./3; THETA[2]=1./3;

    ALPHA[0]=1./3; ALPHA[1]=2./3;
    return 0;
  }
  // SPPRK2(3,2) Spiteri-Ruuth
  else if (sch==18){
    A[0][0]=1./2; 
    A[1][0]=1./2; A[1][1]=1./2;

    THETA[0]=1./3; THETA[1]=1./3; THETA[2]=1./3;

    ALPHA[0]=1./2; ALPHA[1]=1.;
    return 0;
  }
  // SPPRK2(4,2) Spiteri-Ruuth
  else if (sch==19){
    A[0][0]=1./3; 
    A[1][0]=1./3; A[1][1]=1./3;
    A[2][0]=1./3; A[2][1]=1./3;  A[2][2]=1./3;

    THETA[0]=1./4; THETA[1]=1./4; THETA[2]=1./4; THETA[3]=1./4;

    ALPHA[0]=1./3; ALPHA[1]=2./3; ALPHA[2]=1.;
    return 0;
  }
  // SPPRK1(2,1) Spiteri-Ruuth
  else if (sch==20){
    A[0][0]=1./2; 

    THETA[0]=1./2; THETA[1]=1./2; 

    ALPHA[0]=1./2;
    return 0; 
  }
  // SPPRK3(4,3) Spiteri-Ruuth
  else if (sch==21){
    A[0][0]=1./2; 
    A[1][0]=1./2; A[1][1]=1./2;
    A[2][0]=1./6; A[2][1]=1./6;  A[2][2]=1./6;

    THETA[0]=1./6; THETA[1]=1./6; THETA[2]=1./6; THETA[3]=1./2;

    ALPHA[0]=1./2; ALPHA[1]=1.; ALPHA[2]=1./2;
    return 0;
  }
  else{ 
    printf("Erreur choix sch (initButcher)\n");
    return 4;
  }
}


int init_fluides(int tst_multi, int* indFLUIDE, int N, int ind_air, int ind_etain){
  double Ndtemp=N/6;
  int Ntemp=Ndtemp;
  printf("ind_etain= %d, Ndtemp= %g, Ntemp= %d\n",ind_etain,Ndtemp,Ntemp);

  if(tst_multi==0){ // Gaz parfait air
    for (int i = 0; i<N; i++){
      indFLUIDE[i]=ind_air;
    }
  }
  else if(tst_multi==1){  // Bizaruim
    for (int i = 0; i<N; i++){
      indFLUIDE[i]=ind_air;
    }
  }
  else if(tst_multi==2){
    for (int i = 0; i<N; i++){
      indFLUIDE[i]=ind_etain;
    }
  }
  else if(tst_multi==3){
    for(int i=0; i<Ntemp; i++){
      indFLUIDE[i]=ind_air;
      indFLUIDE[N-1-i]=ind_air;
    }
    for(int i=Ntemp; i<N-Ntemp; i++){
      indFLUIDE[i]=ind_etain ;
    }
  }
  else if(tst_multi==4){
    for(int i=0; i<Ntemp; i++){
      indFLUIDE[N-1-i]=ind_air;
    }
    for(int i=0; i<N-Ntemp; i++){
      indFLUIDE[i]=ind_etain ;
    }
  }
  else{printf("Erreur choix tst_multi= %d\n",tst_multi ); return 1;}
  return 0;
}


/*
int maillage(int tst, int N, double fac, int ind_etain, int ind_air, int* indFLUIDE, double* X){
  
  // Maillage
  for(int i=0; i<=N; i++){
    if(indFLUIDE[i]==ind_etain){  X[i]=(a-nbghosts*dx)+dx*i;  }
    if(indFLUIDE[i]==ind_etain){  X[i]=(a-nbghosts*dx)+ dx*i;  }
    //printf("i= %d, X=%lf\n",i,X[i]);
  }

  err=init_fluides(tst_multi, indFLUIDE, N, ind_air, ind_etain);
  if(err){printf("Erreur fluide\n"); return 1;}
*/