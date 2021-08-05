#include <stdio.h>
#include <stdlib.h>
#include <math.h>



// Phase 1
double rho01=7.287e3;
double P01=0;
double T01=300;
double K01=54.73e9;
double N01=5.75;
double sigma01=2.27;
double Cv1=210;

double v01=1./7.287e3;   // 1./rho01;
double E01=0; 
double S01=0;
double F01=0;          // E01 - T01*S01;
 

// Phase 2
double P02=9.4e9;
double T02=300;
double Deltav12=-0.0031e-3;
double dPdT12=-0.017e9;
double K02=94e9;
double N02=4.88;
double sigma02=1.96;
double Cv2=210;

double v02=1.184e-4;   
double E02=123226.473; 
double S02=-65.224;
double F02=142793.615;          // E01 - T01*S01;


// Phase 3
double P03=0;
double T03=505;
double Deltav13=0.004e-3;
double dPdT13=0.03125e9;
double K03=42e9;
double N03=5;
double sigma03=2.25;
double Cv3=200;

double v03=1.431e4;   // 1./rho01;
double E03=118541.806; 
double S03=129.750;
double F03=53018.048;          // E01 - T01*S01;



// Isotherme de Birch 
double Pkb(double v, int phase){
  double v0, K0, N0, P0;
  if(phase==1){
  	v0=v01; K0=K01; N0=N01; P0=P01;
  }
  else if(phase==2){
  	v0=v02; K0=K02; N0=N02; P0=P02;
  }
  else if(phase==3){
  	v0=v03; K0=K03; N0=N03; P0=P03;
  }
  else{printf("Erreur choix phase (Pkb)\n");}

  double res=P0; 
  double fac=(3./2.)*K0*( pow((v/v0),-7./3.) - pow((v/v0),-5./3.) );
  res+=fac*(1. - (3./4.)*(4.-N0)*( pow((v/v0),-2./3.) - 1. ));
  return res;
}

// dérivée de l'Isotherme de Birch 
double Pkb_prime(double v, int phase){
  double v0, K0, N0, P0;
  if(phase==1){
  	v0=v01; K0=K01; N0=N01; P0=P01;
  }
  else if(phase==2){
  	v0=v02; K0=K02; N0=N02; P0=P02;
  }
  else if(phase==3){
  	v0=v03; K0=K03; N0=N03; P0=P03;
  }
  else{printf("Erreur choix phase (Pkb_prime)\n");}

  double fac2=(1. - (3./4.)*(4.-N0)*( pow((v/v0),-2./3.)-1 ));
  double fac1=(3./2)*K0*( pow((v/v0),-7./3) - pow((v/v0),-5./3) );
  double res=fac2*(3./2)*K0*( -(7./3)*(1./v0)*pow((v/v0),-10./3) + (5./3)*(1./v0)*pow((v/v0),-8./3) );
  res+=fac1*(3./4)*(4-N0)*(-2./3)*(1./v0)*pow((v/v0),-5./3);
  return res;
}


// Volume spécifique
double fv(double P, double T, int phase){
  double Cv, sigma0, T0, v0, K0, N0, P0;
  if(phase==1){
  	v0=v01; K0=K01; N0=N01; P0=P01, Cv=Cv1, sigma0=sigma01, T0=T01;
  }
  else if(phase==2){
  	v0=v02; K0=K02; N0=N02; P0=P02, Cv=Cv2, sigma0=sigma02, T0=T02;
  }
  else if(phase==3){
  	v0=v03; K0=K03; N0=N03; P0=P03, Cv=Cv3, sigma0=sigma03, T0=T03;
  }
  else{printf("Erreur choix phase (fv)\n");}

  int n=0;
  double epsilon=1e-8;
  int Nmax=5e5;
  double num, den;

  double v=v0; // initialisation


  // racine de f
  // Algo de Newton
  double delta=1.0;
  while(( fabs(delta)>=epsilon )&&(n<Nmax)){
    num=Pkb(v, phase) + Cv*sigma0*(T-T0)/v0 -P;
    den=Pkb_prime(v, phase);
    delta = num/den;
    v-=delta;
    n++;
  }

  return v; 
}

double Ek(double v, int phase){
  double v0, K0, N0, P0, E0;
  if(phase==1){
  	v0=v01; K0=K01; N0=N01; P0=P01, E0=E01;
  }
  else if(phase==2){
  	v0=v02; K0=K02; N0=N02; P0=P02, E0=E02;
  }
  else if(phase==3){
  	v0=v03; K0=K03; N0=N03; P0=P03, E0=E03;
  }
  else{printf("Erreur choix phase (Ek)\n");}

  double res=E0 - P0*(v-v0); 
  res-=(3./16)*K0*(48-9*N0)*pow(v0,5./3)*(pow(v,-2./3) - pow(v0,-2./3));
  res+=(9./16)*K0*(14-3*N0)*pow(v0,7./3)*(pow(v,-4./3) - pow(v0,-4./3));
  res-=(9./16)*K0*(4-N0)*pow(v0,3)*(pow(v,-2) - pow(v0,-2));
  return res;
}

// Potentiel energie libre
double fF(double v, double T, int phase){
  double Cv, sigma0, T0, v0, K0, N0, P0, E0, S0;
  if(phase==1){
  	v0=v01; K0=K01; N0=N01; P0=P01, Cv=Cv1, sigma0=sigma01, T0=T01, E0=E01, S0=S01;
  }
  else if(phase==2){
  	v0=v02; K0=K02; N0=N02; P0=P02, Cv=Cv2, sigma0=sigma02, T0=T02, E0=E01, S0=S01;
  }
  else if(phase==3){
  	v0=v03; K0=K03; N0=N03; P0=P03, Cv=Cv3, sigma0=sigma03, T0=T03, E0=E01, S0=S01;
  }
  else{printf("Erreur choix phase (fF)\n");}

  double res=E0 + Ek(v, phase);
  res+=Cv*T0*sigma0*(v-v0)/v0;
  res+=Cv*(T-T0);
  res-=T*(S0+Cv*log(T/T0)+Cv*sigma0*(v-v0)/v0);
  return res;
}

// Entropie
double fS(double v, double T, int phase){
  double Cv, sigma0, T0, v0, K0, N0, P0, E0, S0;
  if(phase==1){
  	v0=v01; K0=K01; N0=N01; P0=P01, Cv=Cv1, sigma0=sigma01, T0=T01, E0=E01, S0=S01;
  }
  else if(phase==2){
  	v0=v02; K0=K02; N0=N02; P0=P02, Cv=Cv2, sigma0=sigma02, T0=T02, E0=E01, S0=S01;
  }
  else if(phase==3){
  	v0=v03; K0=K03; N0=N03; P0=P03, Cv=Cv3, sigma0=sigma03, T0=T03, E0=E01, S0=S01;
  }
  else{printf("Erreur choix phase (fS)\n");}

  return S0+Cv*log(T/T0)+Cv*sigma0*(v-v0)/v0;
}

// Pression
double fP(double v, double T, int phase){
  double Cv, sigma0, T0, v0, K0, N0, P0, E0, S0;
  if(phase==1){
  	v0=v01; K0=K01; N0=N01; P0=P01, Cv=Cv1, sigma0=sigma01, T0=T01, E0=E01, S0=S01;
  }
  else if(phase==2){
  	v0=v02; K0=K02; N0=N02; P0=P02, Cv=Cv2, sigma0=sigma02, T0=T02, E0=E01, S0=S01;
  }
  else if(phase==3){
  	v0=v03; K0=K03; N0=N03; P0=P03, Cv=Cv3, sigma0=sigma03, T0=T03, E0=E01, S0=S01;
  }
  else{printf("Erreur choix phase (fP)\n");}

  return Pkb(v, phase) + Cv*sigma0*(v-v0)/v0;
}

// Energie interne
double fE(double v, double T, int phase){
  double Cv, sigma0, T0, v0, K0, N0, P0, E0, S0;
  if(phase==1){
  	v0=v01; K0=K01; N0=N01; P0=P01, Cv=Cv1, sigma0=sigma01, T0=T01, E0=E01, S0=S01;
  }
  else if(phase==2){
  	v0=v02; K0=K02; N0=N02; P0=P02, Cv=Cv2, sigma0=sigma02, T0=T02, E0=E01, S0=S01;
  }
  else if(phase==3){
  	v0=v03; K0=K03; N0=N03; P0=P03, Cv=Cv3, sigma0=sigma03, T0=T03, E0=E01, S0=S01;
  }
  else{printf("Erreur choix phase (fE)\n");}

  double res=E0 + Ek(v, phase);
  res+=Cv*T0*sigma0*(v-v0)/v0;
  res+=Cv*(T-T0);
  return res;
}



// Fonction pour la phase 1
void EOS(double v, double T, int phase, double* P, double* S, double* E, double* F, double* G){

  *P=fP(v, T, phase);
  *S=fS(v, T, phase); 
  *E=fE(v, T, phase);
  *F=(*E)-T*(*S);
  *G=(*F)+(*P)*v;
}




// Matrice jacobienne de G
//   DeltaG[0] : dG/dv
//   DeltaG[1] : dG/dT
//
void JacG(double v, double T, int phase, double* DeltaG){
  double Cv, sigma0, T0, v0, K0, N0, P0, E0, S0;
  if(phase==1){
  	v0=v01; K0=K01; N0=N01; P0=P01, Cv=Cv1, sigma0=sigma01, T0=T01, E0=E01, S0=S01;
  }
  else if(phase==2){
  	v0=v02; K0=K02; N0=N02; P0=P02, Cv=Cv2, sigma0=sigma02, T0=T02, E0=E01, S0=S01;
  }
  else if(phase==3){
  	v0=v03; K0=K03; N0=N03; P0=P03, Cv=Cv3, sigma0=sigma03, T0=T03, E0=E01, S0=S01;
  }
  else{printf("Erreur choix phase (fE)\n");}

  DeltaG[0]=Pkb(v, phase);
  DeltaG[0]-=Cv*sigma0*(T0-T)/v0;
  DeltaG[0]+=v*Pkb_prime(v, phase);
  DeltaG[0]+=fP(v, T, phase);
  
  DeltaG[1]=-fS(v,T, phase);
  DeltaG[1]+=v*Cv*sigma0/v0;
}


// Fonction qui trouve le point triple
/*
  pas : pas constant de la méthode du gradient
  Nmax : nombre max d'itértion
  epsilon : tolérence à zéro
  P0, T0 : point intiale 
  ( ATTENTION : pour calculer G, il me faut v , j'obtient donc 3 valeurs de v pour les 3 phases
                puis je choisi celle qui minimise G du coup pour les 2 autres phases
                le point initial n'est plus le même )

  Valeur retour :
  f : valeur de la fonction
  T,v : racine de f
*/
void triple(double pas, int Nmax, double epsilon, double P0, double T0, double* fres, double* vres, double* Tres){

  printf("ok 1 (triple)\n");

  int n=1;
  int phase;
  double v0;
  double f,v,T;
  /*
  double P0=4.0;
  double T0=600;
  */

  // On calcule v en fct de P et T
  // on a besion de savoir dans quelle phase on est 
  // on regarde celle sui minimise l'entalpie libre de Gibbs
  double P1, P2, P3;
  double E1, E2, E3;
  double S1, S2, S3; 
  double F1, F2, F3;
  double G1, G2, G3;

  double v1=fv(P0, T0, 1);
  double v2=fv(P0, T0, 2);
  double v3=fv(P0, T0, 3);
 
  EOS(v1, T0, 1, &P1, &E1, &S1, &F1, &G1);
  EOS(v2, T0, 2, &P2, &E2, &S2, &F2, &G2);
  EOS(v3, T0, 3, &P3, &E3, &S3, &F3, &G3);

  if(G1<=fmin(G2,G3)){
  	phase=1;
  	v0=v1;
  }
  if(G2<=fmin(G1,G3)){
  	phase=2;
  	v0=v2;
  }
  if(G3<=fmin(G1,G2)){
  	phase=3;
  	v0=v3;
  }

  printf("Phase intiale : %d\n",phase );
  printf("P0 = %.15lf\n",P0*1e-9 );
  printf("P1 = %.15lf, P2 = %.15lf, P3 = %.15lf (GPa)\n",P1*1e-9, P2*1e-9, P3*1e-9 );
  printf("v1 = %.15lf, v2 = %.15lf, v3 = %.15lf\n",v1, v2, v3 );


  // IDEE CALCULE la bonne 
  
  // Cette fonctin calcul la phase à partir de v et T 
  // on ne peut donc pas l'appliquer au dessus car on ne connait que P et T 
  // il nous faut la phase pour connaitre v
  ////////////////////////////////////////////////////////////////////////////
  /*
  int cal_phase(v,T){
  	  double *P1, *P2, *P3;
	  double *E1, *E2, *E3;
	  double *S1, *S2, *S3; 
	  double *F1, *F2, *F3;
	  double *G1, *G2, *G3;
	 
	  EOS(v, T, P1, 1, E1, S1, F1, G1);
	  EOS(v, T, P2, 2, E2, S2, F2, G2);
	  EOS(v, T, P3, 3, E3, S3, F3, G3);

	  if(G1<=fmin(G2,G3)){
	  	phase=1;
	  }
	  else if(G2<=fmin(G1,G3)){
	  	phase=2;
	  }
	  else if(G3<=fmin(G1,G2)){
	  	phase=3;
	  }
	  else{printf("Erreur calcul minimum (cal_phase)\n");}

	  return phase;
  } 
  */
  ///////////////////////////////////////////////////////////////

  // Initialisation
  double *DeltaG1=malloc(2*sizeof(double));
  double *DeltaG2=malloc(2*sizeof(double));
  double *DeltaG3=malloc(2*sizeof(double));

  v=v0;
  T=T0;

  EOS(v, T, 1, &P1, &E1, &S1, &F1, &G1);
  EOS(v, T, 2, &P2, &E1, &S1, &F1, &G2);
  EOS(v, T, 3, &P3, &E1, &S1, &F1, &G3);

  printf("Pression retenue pour chaque phase après selection d'un unique v : \n" );
  printf(" P0 = %.15lf\n",P0*1e-9 );
  printf(" P1 = %.15lf, P2 = %.15lf, P3 = %.15lf  (GPa)\n", P1*1e-9, P2*1e-9, P3*1e-9 );

  f=fabs(G1-G2)+fabs(G1-G3);

  double sgn1, sgn2;

  
  double df=1.;

  while((df>=epsilon)&&(n<=Nmax)){
    JacG(v, T, 1, DeltaG1);
  	JacG(v, T, 2, DeltaG2);
  	JacG(v, T, 3, DeltaG3);

  	EOS(v, T, 1, &P1, &E1, &S1, &F1, &G1);
    EOS(v, T, 2, &P1, &E1, &S1, &F1, &G2);
    EOS(v, T, 3, &P1, &E1, &S1, &F1, &G3);

    if(G1>G2){ sgn1=1; }
	else{	sgn1=-1; }

    if(G1>G3){ sgn2=1; }
    else{	sgn2=-1; }


    v-=pas*(sgn1*(DeltaG1[0]-DeltaG2[0]) + sgn2*(DeltaG1[0]-DeltaG3[0]));
    T-=pas*(sgn1*(DeltaG1[1]-DeltaG2[1]) + sgn2*(DeltaG1[1]-DeltaG3[1]));
    

    df=fabs(f-fabs(G1-G2)-fabs(G1-G3) );
    f=fabs(G1-G2)+fabs(G1-G3);

    printf("n=%d, df= %.15lf, f=%.15lf\n",n,df,f );

    n++;
  }

  printf("Nombre di'itérations = %d\n",n );
  printf("df = %.15lf\n",df );
  printf("f= %.15lf\n",f );

  *vres=v;
  *fres=f;
  *Tres=T;


}


// programme main pour trouver le point triple
int main(){
  printf("ok1 (main)\n");
  double pas=1e-5; 
  int Nmax=1e4; 
  double epsilon=1e-8;
  double P0=4e9;
  double T0=600; 
  
  double f, v, T;

  printf("debut triple (main)\n");
   
  triple(pas, Nmax, epsilon, P0, T0, &f, &v, &T);

  printf("sortie triple (main)\n");

  printf("f=%.15lf, v=%.15lf, T=%.15lf\n",f,v,T );


  double P1, P2, P3;
  double E1, E2, E3;
  double S1, S2, S3; 
  double F1, F2, F3;
  double G1, G2, G3;

  double v1,v2,v3;

  v1=fv(P03, T03, 1);
  v3=fv(P03, T03, 3);
 
  EOS(v1, T03, 1, &P1, &E1, &S1, &F1, &G1);
  EOS(v3, T03, 3, &P3, &E3, &S3, &F3, &G3);

  printf("P03= %.15lf, T03= %.15lf\n",P03*1e-9,T03 );
  printf("v1= %.15lf, v3= %.15lf\n", v1, v3);
  printf("P1= %.15lf, P3= %.15lf\n", P1*1e-9, P3*1e-9);
  printf("G1= %.15lf, G3= %.15lf\n", G1, G3);


  v1=fv(P02, T02, 1);
  v2=fv(P02, T02, 3);
 
  EOS(v1, T02, 1, &P1, &E1, &S1, &F1, &G1);
  EOS(v2, T02, 3, &P2, &E2, &S2, &F2, &G2);

  printf("P02= %.15lf, T02= %.15lf\n",P02*1e-9,T02 );
  printf("v1= %.15lf, v2= %.15lf\n", v1, v2);
  printf("P1= %.15lf, P2= %.15lf\n", P1*1e-9, P2*1e-9);
  printf("G1= %.15lf, G2= %.15lf\n", G1, G2);


  return 0;
}