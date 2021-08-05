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


// Phase 3
double P03=0;
double T03=505;
double Deltav13=0.004e-3;
double dPdT13=0.03125e9;
double K03=42e9;
double N03=5;
double sigma03=2.25;
double Cv3=200;




// Isotherme de Birch 
double Pkb(double v, int phase){
  double v0, K0, N0, P0;
  if(phase==1){
  	v0=v01; K0=K01; N0=N01; P0=P01;
  }
  else{printf("Erreur choix phase (Pkb)\n");}

  double res=P0; 
  double fac=(3./2.)*K0*( pow((v/v0),-7./3.) - pow((v/v0),-5./3.) );
  res+=fac*(1.0 - (3./4.)*(4.0-N0)*( pow((v/v0),-2./3.) - 1.0 ));
  return res;
}

// dérivée de l'Isotherme de Birch 
double Pkb_prime(double v, int phase){
  double v0, K0, N0, P0;
  if(phase==1){
  	v0=v01; K0=K01; N0=N01; P0=P01;
  }
  else{printf("Erreur choix phase (Pkb_prime)\n");}

  double fac2=(1 - (3./4)*(4-N0)*( pow((v/v0),-2./3)-1 ));
  double fac1=(3./2)*K0*( pow((v/v0),-7./3) - pow((v/v0),-5./3) );
  double res=fac2*(3./2)*K0*( -(7./3)*(1./v0)*pow((v/v0),-10./3) + (5./3)*(1./v0)*pow((v/v0),-8./3) );
  res+=fac1*(3./4)*(4-N0)*(-2./3)*(1./v0)*pow((v/v0),-5./3);
  return res;
}


// Volume spécifique
double fv(double P, double T, int phase){
  printf("fv :\n");
  double Cv, sigma0, T0, v0, K0, N0, P0;
  if(phase==1){
  	v0=v01; K0=K01; N0=N01; P0=P01, Cv=Cv1, sigma0=sigma01, T0=T01;
  }
  else{printf("Erreur choix phase (fv)\n");}

  int n=0;
  double epsilon=1e-8;
  int Nmax=5e2;
  double num, den;

  double v=v0; // initialisation


  /*
  // fonction à minimiser
  double f(double v){
    return Pkb(v, v0, K0, N0, P0) + Cv*sigma0*(T-T0)/v0 -P;
  }
  // fonction dérivée
  double fd(double v){
    return Pkb_prime(v, v0, K0, N0, P0);
  }
  */

  // racine de f
  // Algo de Newton
  double delta=1.0;
  while(( fabs(delta)>=epsilon )&&(n<Nmax)){
    num=Pkb(v, phase) + Cv*sigma0*(T-T0)/v0 - P;
    den=Pkb_prime(v, phase);
    delta = num/den;
    v-=delta;
    n++;
    printf(" n= %d, delta= %.10lf, v= %.10lf, f= %.15lf \n",n,delta,v, Pkb(v, phase) + Cv*sigma0*(T-T0)/v0 - P);
  }

  printf(" n= %d, delta= %.10lf, v= %.10lf, f= %.15lf \n",n,delta,v, Pkb(v, phase) + Cv*sigma0*(T-T0)/v0 - P);


  return v;
}

double Ek(double v, int phase){
  double v0, K0, N0, P0, E0;
  if(phase==1){
  	v0=v01; K0=K01; N0=N01; P0=P01, E0=E01;
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
  else{printf("Erreur choix phase (fS)\n");}

  return S0+Cv*log(T/T0)+Cv*sigma0*(v-v0)/v0;
}

// Pression
double fP(double v, double T, int phase){
  double Cv, sigma0, T0, v0, K0, N0, P0, E0, S0;
  if(phase==1){
  	v0=v01; K0=K01; N0=N01; P0=P01, Cv=Cv1, sigma0=sigma01, T0=T01, E0=E01, S0=S01;
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
  else{printf("Erreur choix phase (fE)\n");}

  double res=E0 + Ek(v, phase);
  res+=Cv*T0*sigma0*(v-v0)/v0;
  res+=Cv*(T-T0);
  return res;
}


// Fonction qui calcule les coefficient pour les phases 2 et 3
void Coeff(void){
  double v02=fv(P02, T02, 1) + Deltav12;
  double S02=fS(v02, T02, 1); 
  double E02=fE(v02, T02, 1) - P02*Deltav12 + T02*Deltav12*dPdT12;
  double F02=E02 - T02*S02;

  printf("Phase 2:\n");
  printf(" v02= %.15lf\n",v02 );
  printf(" S02= %.15lf\n",S02 );
  printf(" E02= %.15lf\n",E02 );
  printf(" F02= %.15lf\n",F02 );
 
  double v03=fv(P03, T03, 1) + Deltav13;
  double S03=fS(v03, T03, 1); 
  double E03=fE(v03, T03, 1) - P03*Deltav13 + T03*Deltav13*dPdT13;
  double F03=E03 - T03*S03;

  printf("Phase 3:\n");
  printf(" v03= %.15lf\n",v03 );
  printf(" S03= %.15lf\n",S03 );
  printf(" E03= %.15lf\n",E03 );
  printf(" F03= %.15lf\n",F03 );
  
  // Valeur de la fonction à annuler pour obtenir v03
  double f,v;

  printf("cst= %.15lf\n",Cv1*sigma01*(T03-T01)/v01 );


  FILE* Res;

  if((Res = fopen("pkb.txt", "w+")) != NULL){
    for(int i=-1e2; i<=1e4; i++){ 
      v=v01 + i*1e-6;
      f=Pkb(v, 1) + Cv1*sigma01*(T03-T01)/v01 - P03;
      fprintf(Res,"%.15lf %.15lf\n",v/v01,Pkb(v, 1) );
    }
  }
  else{printf("Erreur ouverture pkb.txt\n");}

}



// programme main pour trouver le point triple
int main(){
   
  Coeff();


  return 0;
}