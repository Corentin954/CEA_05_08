#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// Isotherme Inversible

// Phase 1
double K01=5.473e10;
double N01=5.75;
double gamma01=2.27;
double Cv1=210.;
double theta01=250.;
double T01=300.;
double P01=0.;

double rho01=7287.;   // 1./rho01;
double v01=1./7287.;   // 1./rho01;
double E01=0.; 
double Sr1=-38.2875;
//double F03= 300*38.2875;          // E0 - T0*S0;
 

// Phase 2
double K02=5.9e10;
double N02=4.8;
double gamma02=1.96;
double Cv2=210.;
double theta02=300.;
double T02=300.;
double P02=9.4e9;
double Deltav12=-1.9e-6;
double dPdT12=-1.7e7;

double rho02=7401.45;   // 1./rho01;
double v02=1./7401.45;   // 1./rho01;
double E02=39033.; 
double Sr2=24.7813;
//double F03=39033. - 300*24.7813;          // E0 - T0*S0;


// Phase 3
double K03=4.8e10;
double N03=5.;
double gamma03=2.25;
double Cv3=200.;
double theta03=350.;
double T03=505.;
double P03=0.;
double Deltav13=4e-6;
double dPdT13=3.125e7;

double rho03=7026.17;   // 1./rho01;
double v03=1./7026.17;   // 1./rho01;
double E03=90660.5; 
double Sr3=165.11;
//double F03=90660.5 - 505*165.11;          // E0 - T0*S0;

double ur=1.2;

// Coeff MABIRE
/*
K02=94e9;
N02=4.88;

K03=42e9;
*/




void coeff(int phase, double* K0, double* N0, double* gamma0, double* Cvr, double* theta0, double* T0, double* P0, double* rho0, double* v0, double* E0, double* Sr){
  if(phase==1){
    *gamma0=gamma01; *rho0=rho01; *Cvr=Cv1; *theta0=theta01; *E0=E01; *K0=K01; *N0=N01;
  }
  else if(phase==2){
    *gamma0=gamma02; *rho0=rho02; *Cvr=Cv2; *theta0=theta02; *E0=E02; *K0=K02; *N0=N02; 
  }
  else if(phase==3){
    *gamma0=gamma03; *rho0=rho03; *Cvr=Cv3; *theta0=theta03; *E0=E03; *K0=K03; *N0=N03;
  }
  else{printf("Erreur choix phase (coeff)\n");}
}


// x
double fx(double V, int phase){
  double rho0;
  if(phase==1){
    rho0=rho01; 
  }
  else if(phase==2){
    rho0=rho02;
  }
  else if(phase==3){
    rho0=rho03;
  }
  else{printf("Erreur choix phase (fx)\n");}

  return 1-rho0*V;
}

// THETA
double THETA(double V, int phase){
  double gamma0, theta0;
  if(phase==1){
    gamma0=gamma01; theta0=theta01;
  }
  else if(phase==2){
    gamma0=gamma02; theta0=theta02;
  }
  else if(phase==3){
    gamma0=gamma03; theta0=theta03;
  }
  else{printf("Erreur choix phase (THETA)\n");}

  
  return theta0*exp( gamma0*fx(V,phase) );
}

// Es
double fEs(double V, int phase){
  double E0, K0, rho0, N0;
  if(phase==1){
    E0=E01; K0=K01; rho0=rho01; N0=N01;
  }
  else if(phase==2){
    E0=E02; K0=K02; rho0=rho02; N0=N02;
  }
  else if(phase==3){
    E0=E03; K0=K03; rho0=rho03; N0=N03;
  }
  else{printf("Erreur choix phase (fEs)\n");}

  double x=fx(V,phase);

  double res=E0;
  double fac=K0/(rho0*(N0+1.0));
  res+= fac*( exp( (N0+1.0)*x ) -1.0 )/(N0+1.0);
  res-= fac*x;

  return res;
}


// Ps
double fPs(double V, int phase){
  double E0, K0, N0;
  if(phase==1){
    E0=E01; K0=K01; N0=N01;
  }
  else if(phase==2){
    E0=E02; K0=K02; N0=N02;
  }
  else if(phase==3){
    E0=E03; K0=K03; N0=N03;
  }
  else{printf("Erreur choix phase (fPs)\n");}


  double fac=K0/(N0+1);
  return fac*( exp( (N0+1)*fx(V,phase) ) -1.0 );
}

// u
double u(double V, double E, int phase){
  double Cvr;
  if(phase==1){
  	Cvr=Cv1;
  }
  else if(phase==2){
  	Cvr=Cv2;
  }
  else if(phase==3){
  	Cvr=Cv3;
  }
  else{printf("Erreur choix phase (u)\n");}

  double res=ur;
  res+=( E-fEs(V,phase) )/( Cvr*THETA(V,phase) ); 

  return res;
}

// S
double fS(double V, double E, int phase){
  double Sr, Cvr;
  if(phase==1){
    Sr=Sr1; Cvr=Cv1;
  }
  else if(phase==2){
    Sr=Sr2; Cvr=Cv2;
  }
  else if(phase==3){
    Sr=Sr3; Cvr=Cv3;
  }
  else{printf("Erreur choix phase (fS)\n");}

  double res=Sr;
  res+= Cvr*log(u(V,E,phase));

  return res;
}


// P
double fP(double V, double E, int phase){
  double gamma0, rho0;
  if(phase==1){
    gamma0=gamma01; rho0=rho01;
  }
  else if(phase==2){
    gamma0=gamma02; rho0=rho02;
  }
  else if(phase==3){
    gamma0=gamma03; rho0=rho03;
  }
  else{printf("Erreur choix phase (fP)\n");}

  double res=fPs(V,phase);
  res+= gamma0*rho0*( E - fEs(V,phase) );

  return res;
}


// T
double fT(double V, double E, int phase){
  double Cvr;
  if(phase==1){
    Cvr=Cv1;
  }
  else if(phase==2){
    Cvr=Cv2;
  }
  else if(phase==3){
    Cvr=Cv3;
  }
  else{printf("Erreur choix phase (fT)\n");}

  double res=(E-fEs(V,phase))/Cvr;
  res+= ur*THETA(V,phase);

  return res;
}


// Es'
double Es_prime(double V, int phase){
  double K0, N0;
  if(phase==1){
    K0=K01; N0=N01;
  }
  else if(phase==2){
    K0=K02; N0=N02;
  }
  else if(phase==3){
    K0=K03; N0=N03;
  }
  else{printf("Erreur choix phase (Es_prime)\n");}

  double x=fx(V,phase);

  double fac=K0/(N0+1);
  double res= fac*( 1 - exp( (N0+1)*x ) );

  return res;
}

// Ps'
double Ps_prime(double V, int phase){
  double K0, N0, rho0;
  if(phase==1){
    K0=K01; N0=N01; rho0=rho01;
  }
  else if(phase==2){
    K0=K02; N0=N02; rho0=rho02;
  }
  else if(phase==3){
    K0=K03; N0=N03; rho0=rho03;
  }
  else{printf("Erreur choix phase (Ps_prime)\n");}

  double x=fx(V,phase);

  return -rho0*K0*exp( (N0+1)*x );
}


// dP/dE
double dPdE(int phase){
  double gamma0, rho0;
  if(phase==1){
    gamma0=gamma01; rho0=rho01;
  }
  else if(phase==2){
    gamma0=gamma02; rho0=rho02;
  }
  else if(phase==3){
    gamma0=gamma03; rho0=rho03;
  }
  else{printf("Erreur choix phase (dPdE)\n");}

  return gamma0*rho0;
}


// dPdV
double dPdV(double V, int phase){
  double gamma0, rho0;
  if(phase==1){
    gamma0=gamma01; rho0=rho01;
  }
  else if(phase==2){
    gamma0=gamma02; rho0=rho02;
  }
  else if(phase==3){
    gamma0=gamma03; rho0=rho03;
  }
  else{printf("Erreur choix phase (dPdV)\n");}

  return Ps_prime(V,phase) - gamma0*rho0*Es_prime(V,phase);
}

// dT/dE
double dTdE(int phase){
  double Cvr;
  if(phase==1){
    Cvr=Cv1; 
  }
  else if(phase==2){
    Cvr=Cv2; 
  }
  else if(phase==3){
    Cvr=Cv3; 
  }
  else{printf("Erreur choix phase (dTdE)\n");}

  return 1./Cvr;
}


// dTdV
double dTdV(double V, int phase){
  double gamma0, rho0, Cvr, theta0;
  if(phase==1){
    gamma0=gamma01; rho0=rho01; Cvr=Cv1; theta0=theta01;
  }
  else if(phase==2){
    gamma0=gamma02; rho0=rho02; Cvr=Cv2; theta0=theta02;
  }
  else if(phase==3){
    gamma0=gamma03; rho0=rho03; Cvr=Cv3; theta0=theta03;
  }
  else{printf("Erreur choix phase (dTdV)\n");}

  double x=fx(V,phase);
  
  double res=(-1./Cvr)*Es_prime(V,phase);
  res-=ur*theta0*gamma0*rho0*exp(gamma0*x);
  return res;
}



// EPSILON fnction de V et P  (pour l'initailisation dans les schemas)
double fE_VP(double V, double P, int phase){
  double gamma0, rho0, Cvr, theta0;
  if(phase==1){
    gamma0=gamma01; rho0=rho01; Cvr=Cv1; theta0=theta01;
  }
  else if(phase==2){
    gamma0=gamma02; rho0=rho02; Cvr=Cv2; theta0=theta02;
  }
  else if(phase==3){
    gamma0=gamma03; rho0=rho03; Cvr=Cv3; theta0=theta03;
  }
  else{printf("Erreur choix phase (fE_VP)\n");}


  double res=fEs(V,phase);
  res+= (P-fPs(V,phase))/(gamma0*rho0);
  return res;
}




/*
// fonction de maillage de la zone A ( qui a un nombre impair de sommets)
/* Arguments :
  - int N12, N13.
  - double* P12, T12 : tableaux de N12 éléments qui contients les points de la courbe C12
                       partant du point triple
  - double* P13, T13 : tableaux de N13 éléments qui contients les points de la courbe C13
                       partant du point triple

  Arguments de soities: void
  - la fonction va écrire dans un fichier les sommets et dans l'autre les mailles 
*/
/*
int maillageA(int N12, int N13, double* P12, double* T12, double* P13, double* T13){

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


  double Pref, Tref;

  Pref=P13[N13-1];
  Tref=T12[N12-1];


  // **********   ECRITURE DES SOMMETS   *************
  // front en bas à gauche pression de C13
  for(int i=N13-1; i>=0; i--){
    printf("  --i= %d\n",i );
    P=P13[i];
    T=Tref;
    fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (maillageA)\n");return err;}
    fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
  }
  // front en bas à gauche mais pression de C12
  for(int i=1; i<N12; i++){
    printf("  --i= %d\n",i );
    P=P12[i]; 
    T=Tref;
    fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (maillageA)\n");return err;}
    fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
  }
  // On remonte jusqu'a point triple par C12
  for(int i=N12-1; i>=0; i--){
    printf("  --i= %d\n",i );
    P=P12[i]; 
    T=T12[i];
    fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (maillageA)\n");return err;}
    fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
  }
  // On suit C13 jusqu'à la limite gauche
  for(int i=1; i<N13; i++){
    printf("  --i= %d\n",i );
    P=P13[i];
    T=T13[i];
    fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (maillageA)\n");return err;}
    fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
  }
  


  int a,b,c,d;
  */
/* maille : (num des points coresp.)
  d    c

  a    b
*/
  // **********   ECRITURE DES MAILLES   *************
  /*
  for(int i=1; i <=n; i++){
    a=i; 
    b=i+1; 
    c=N-i;
    d=N-i+1;
    fprintf(fmaill, "%d %d %d %d\n",a,b,c,d );
  }
  */
/* maille triangle : (num des points coresp.)
  c

  a    b
*/
  /*
  a=n; 
  b=n+1; 
  c=n+2;
  fprintf(fmaill, "%d %d %d\n",a,b,c );


  fclose(fptsPT); 
  fclose(fptsVE); 
  fclose(fmaill);
}
*/

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
/*
int maillageB(int N12, int N23, double* P12, double* T12, double* P23, double* T23){

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


  double Pref, Tref;

  Tref=T12[N12-1];


  // **********   ECRITURE DES SOMMETS   *************
  // Parcours de C12
  for(int i=0; i<N12; i++){
    P=P12[i];
    T=T12[i];
    fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (maillageB)\n");return err;}
    fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
  }
  // front en bas à droite mais pression de C23
  for(int i=N12; i<N23; i++){
    P=P23[i]; // +1 pour eviter le points triple
    T=Tref;
    fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (maillageB)\n");return err;}
    fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
  }
  // On remonte jusqu'a point triple par C12
  for(int i=0; i<N12-1; i++){
    P=P12[N12-2-i]; // +1 pour eviter le points triple
    T=T12[N12-2-i];
    fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (maillageB)\n");return err;}
    fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
  }
  // On suit C23 jusqu'au point triple (sens inverse)
  for(int i=N23-1; i>=2; i--){
    P=P23[i]; // +1 pour eviter le points triple
    T=T23[i];
    fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (maillageB)\n");return err;}
    fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
  }
  */

// **********   ECRITURE DES MAILLES   *************

  //int a,b,c,d;
/* maille triangle : (num des points coresp.)
  c

  a    b
*/
  /*
  a=1; 
  b=2; 
  c=N;
  fprintf(fmaill, "%d %d %d\n",a,b,c );
*/
/* maille carre : (num des points coresp.)
  d    c

  a    b
*/
  /*
  for(int i=2; i<=n; i++){
    a=i;
    b=i+1; 
    c=N-i+1;
    d=N-i+2;
    fprintf(fmaill, "%d %d %d %d\n",a,b,c,d );
  }
  */

/* maille triangle : (num des points coresp.)
   c    b              c
             et           
   a              a    b
*/
  /*
  for(int i=1; i<=n-1; i++){
    a=i;
    b=N-i;
    c=N-i+1;
    fprintf(fmaill, "%d %d %d\n",a,b,c);
    a=i;
    b=i+1;
    c=N-i;
    fprintf(fmaill, "%d %d %d\n",a,b,c);
  }


  fclose(fptsPT); fclose(fptsVE);
  fclose(fmaill);
}
*/
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
/*
int maillageC(int N13, int N23, double* P13, double* T13, double* P23, double* T23){

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


  double Pref, Trefhaut;
  Trefhaut=2500;


  // **********   ECRITURE DES SOMMETS   *************
  // Parcours de C13 jausuq'au point triple (dans le sens inverse)
  for(int i=N13-1; i>=0; i--){
    P=P13[i];
    T=T13[i];
    fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (maillageC)\n");return err;}
    fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
  }
  // Parcours de C23 jausuq'à la frontière droite (dans le bon sens)
  for(int i=1; i<N23; i++){
    P=P23[i];
    T=T23[i];
    fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (maillageC)\n");return err;}
    fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
  }
  // Parcours la fontière haute de la droite vers la presion du point triple
  for(int i=N23-1; i>=0; i--){
    P=P23[i];
    T=Trefhaut;
    fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (maillageC)\n");return err;}
    fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
  }
  // Parcours la fontière haute du point triple vers la gauche
  for(int i=0; i<N13; i++){
    P=P13[i];
    T=Trefhaut;
    fprintf(fptsPT, "%.15lf %.15lf\n",P,T );
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (maillageC)\n");return err;}
    fprintf(fptsVE, "%.15lf %.15lf\n",V,E );
  }
  

  int a,b,c,d;
  */
/* maille carre : (num des points coresp.)
  d    c

  a    b
*/
  // **********   ECRITURE DES MAILLES   *************
  /*
  for(int i=1; i<=n-1; i++){
    a=i;
    b=i+1; 
    c=N-i;
    d=N-i+1;
    fprintf(fmaill, "%d %d %d %d\n",a,b,c,d );
  }
  */


/* maille triangle : (num des points coresp.)
   c    b              c
             et    
   a              a    b
*/
  /*
  // **********   ECRITURE DES MAILLES   *************
  for(int i=1; i<=n-1; i++){
    a=i;
    b=N-i; 
    c=N-i+1;
    fprintf(fmaill, "%d %d %d\n",a,b,c);
    a=i;
    b=i+1; 
    c=N-i+1;
    fprintf(fmaill, "%d %d %d\n",a,b,c);
  }

  fclose(fptsPT); fclose(fptsVE); 
  fclose(fmaill);
}
*/