#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "head.h"


// CALCUL DE P et U (p12, u12 || p^*, u^*) POUR LA SOLUTION ANALYTIQUE 
/*
 On lit un fichier qui contient les valeurs de wL et wR
 puis on écrit P et U dans un fichier
*/



// compilation :  compil_PU.sh 

double Gamma=1.4;


/*
Fonction qui lie P et U
lindice k corespond aux valeurs a gauche (k=1) ou à droite (k=2)
*/
double f(double P, double pk, double rhok){
  double res, ck, pik;
  pik=P/pk;
  ck=sqrt(Gamma*pk/rhok);
  if (P>= pk){
    res = (P-pk)/( rhok*ck*sqrt( (Gamma+1)*pik/(2*Gamma) + (Gamma-1)/(2*Gamma) ) );
  }
  else{
  	res = 2*ck*(pow(pik,(Gamma-1)/(2*Gamma)) -1 )/(Gamma-1);
  }
  return res;
}

/* Dérivee de f */
double fd(double P, double pk, double rhok){
  double num, den;
  double ck, pik;
  pik=P/pk;
  ck=sqrt(Gamma*pk/rhok);
  if (P>= pk){
    num =  (Gamma+1)*pik + 3*Gamma-1;
    den =  4*Gamma*rhok*ck*sqrt( pow( (Gamma+1)*pik/(2*Gamma) + (Gamma-1)/(2*Gamma), 3) ) ;
  }
  else{
  	num = ck*pow(pik,(Gamma-1)/(2*Gamma));
  	den = Gamma*P;
  }
  return num/den;
}


// Fonction flux massique a
double a(double P, double pk, double rhok){
  double num, den;
  double ck, pik;
  pik=P/pk;
  ck=sqrt(Gamma*pk/rhok);
  if (P>=pk){
  	num=sqrt( rhok*( (Gamma+1)*P/2 + (Gamma-1)*pk/2 ) );  // formule 13.4 / 13.6
  	den=1.0;
  }
  else{
  	num = ((Gamma-1)/(2*Gamma))*rhok*ck*(1 - pik);
  	den = 1 - pow( P/pk,(Gamma-1)/(2*Gamma) ) ;            // formule 13.11 / 13.12
  }
  return num/den;
}


int main(){
  double epsilon =1e-12;
  int n=0, Nmax=1e4;

// donnes initiale
  double rhoG, rhoD;
  double uG, uD;
  double pG, pD;

  double cG, cD;

  double P0, delta;
  double num, den;
  double aG, aD; // equivalent du flux massique 

// valeur à trouver 
  double P, U;
  
  rhoG=1.0; rhoD=0.125;
  uG=0; uD=0;
  pG=1.0; pD=0.1;

  cG=sqrt(Gamma*pG/rhoG); 
  cD=sqrt(Gamma*pD/rhoD);

// Initialisation P0 
  P0 = (pG*rhoD*cD + pD*rhoG*cG + (uG-uD)*rhoG*cG*rhoD*cD )/( rhoG*cG + rhoD*cD );
  P=P0;
  delta=1;
// Algo de Newton
  while(( fabs(delta)>=epsilon )&&(n<Nmax)){
    num = f(P,pG,rhoG) + f(P,pD,rhoD) - (uG - uD);
    den = fd(P,pG,rhoG) + fd(P,pD,rhoD);
    delta = num/den;
    P-=delta;
    n++;
  }
  /*
  printf("\n Newton :\n");
  printf("delta = %lf et n= %d\n",delta,n );
  */

  
// Calcul de U
  U = (a(P,pG, rhoG)*uG + a(P,pD, rhoD)*uD + pG - pD)/(a(P,pG, rhoG) + a(P,pD, rhoD));


  // ecriture dans un fichier des points (x,u(x)) de la solution
  FILE* Res;
  if((Res = fopen("PU.txt", "w+")) != NULL){
    //fprintf(Res,"%lf %lf",P, U);
    fprintf(Res,"%0.15f %0.15f",P, U);
    fclose(Res);
  } 
  else{printf("Erreur d'ouverture de sol.txt\n");}

}