#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "head.h"

/*  /////////////////////////////////////////////////////////////////////////////////

           Equation d'etat du GAZ PARFAIT

p = (gamma-1)*epsilon/tau
epsilon = tau/((gamma-1)*epsilon)
*/


//   fonction ENERGIE INTERNE SPECIFIQUE epsilon(tau, p)
double epsilonGP(double Gamma, double tau, double p){
  return tau*p/(Gamma-1) ;
}

//   fonction PRESSION p + VITESSE DU SON c(tau, p)
// renvoie ¨Pres et Cres
double egsGP(double Gamma, double tau, double epsilon, double *Pres, double* Cres){
  *Pres=(Gamma-1)*epsilon/tau;
  *Cres= sqrt(Gamma*(*Pres)*tau) ;
}


/* /////////////////////////////////////////////////////////////////////////////////
    
    Equation d'etat du BIZARIUM 

p = p_k0 + (sigma0/tau0).(epsilon - epsilon_k0)
epsilon = epsilon_k0 + (tau0/sigma0).(p - p_k0)  
*/

// données 
double rho0 = 1e4; double tau0=1./1e4;
double K0 = 1e11;
double Cv0 = 1e3;
double T0 = 300;
double eps0 = 0;
double sigma0 = 1.5;
double s = 1.5;
double q = -(42080895.0)/(14941154.0);
double r = (727668333.0)/(149411540.0);



//   fonction ENERGIE INTERNE SPECIFIQUE epsilon(tau, p)
double epsilonBIZ(double tau, double p){
  double x= tau0/tau - 1.0;

  double gEOS= sigma0*(1.0 - tau/tau0);

  double f0 = 1.0 + (s/3.0 - 2.0)*x + q*x*x + r*x*x*x;
  f0=f0/(1.0-s*x);

  double f1 = s/3.0 - 2.0 + 2.0*q*x + 3.0*r*x*x + s*f0;
  f1=f1/(1.0-s*x);


  double epsilon_k0= eps0;
  epsilon_k0-=Cv0*T0*(1.0+ gEOS);
  epsilon_k0+=K0*tau0*x*x*f0/2.0;


  double p_k0= -Cv0*T0*sigma0*rho0;
  p_k0+= K0*x*(1.0+x)*(1.0+x)*(2.0*f0 + x*f1)/2.0;


  return epsilon_k0 + (tau0/sigma0)*(p - p_k0);
}

//   fonction PRESSION p + VITESSE DU SON c(tau, p)
// renvoie ¨Pres et Cres
double egsBIZ(double tau, double epsilon, double* Pres, double* Cres){

  double x= tau0/tau - 1.0;

  double gEOS= sigma0*(1.0 - tau/tau0);

  double f0 = 1.0 + (s/3.0 - 2.0)*x + q*x*x + r*x*x*x;
  f0=f0/(1.0-s*x);

  double f1 = s/3.0 - 2.0 + 2.0*q*x + 3.0*r*x*x + s*f0;
  f1=f1/(1.0-s*x);

  double f2 = 2.0*q + 6.0*r*x + 2.0*s*f1;
  f2=f2/(1.0-s*x);

  double epsilon_k0= eps0;
  epsilon_k0-=Cv0*T0*(1.0+ gEOS);
  epsilon_k0+=K0*tau0*x*x*f0/2.0;


  double p_k0= -Cv0*T0*sigma0*rho0;
  p_k0+= K0*x*(1.0+x)*(1.0+x)*(2.0*f0 + x*f1)/2.0;

  double fac=-K0*rho0*(1.0+x)*(1.0+x)*(1.0+x)/2.0;
  double pd_k0= 2.0*(1.0+3.0*x)*f0;
  pd_k0+= 2.0*x*(2.0+3.0*x)*f1;
  pd_k0+= x*x*(1.0+x)*f2;
  pd_k0=pd_k0*fac;
 

  // resultats
  *Pres= p_k0 + (sigma0/tau0)*(epsilon - epsilon_k0);

  *Cres= tau*sqrt( (sigma0/tau0*(*Pres-p_k0)-pd_k0) ) ;
}





//////////////////////////////////////////////////////////////////////
//                     fonction qui "sortent" du fichier
/////////////////////////////////////////////////////////////////////


//   fonction ENERGIE INTERNE SPECIFIQUE epsilon(tau, p)
double epsilonEOS(int tst, double tau, double p){
  double Gamma;
  if(tst==0|| tst==4 || tst==3 || tst==5 || tst==6){
    Gamma=1.4;
    return epsilonGP(Gamma,tau,p);
  }
  else if(tst==2){
    return epsilonBIZ(tau,p);
  }
  else if(tst==1){
    Gamma=5./3;
    return epsilonGP(Gamma,tau,p);
  }
  else{
    printf("Erreur de choix de l'equation d'etat (epsilonEOS)\n");
  }
}

//   fonction PRESSION p + VITESSE DU SON c(tau, p)
// renvoie ¨Pres et Cres
void EGS(int tst, double tau, double epsilon, double* Pres, double* Cres){
  double Gamma; 
  if(tst==0 || tst==4 || tst==3 || tst==5 || tst==6){
    Gamma=1.4;
    egsGP(Gamma, tau,epsilon, Pres, Cres);
  }
  else if(tst==2){
    egsBIZ(tau,epsilon, Pres, Cres);
  }
  else if(tst==1){
    Gamma=5./3;
    egsGP(Gamma, tau,epsilon, Pres, Cres);
  }
  else{
    printf("Erreur de choix de l'equation d'etat (EGS)\n");
  }
}