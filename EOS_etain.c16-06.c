#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Matrice.h"


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

double **alloctabd(int dim1, int dim2) ;
void freetab(void *ptr);

/*  Valeur d'entree :
   Va,Ea, Vb,Eb, Vc,Ec

   corespodances :  1-a   2-b   3-c

   Valeur sortie
  VA,EA, VB,EB, VC,EC

*/
void Newton_triple(double Va, double Ea, double Vb, double Eb, double Vc, double Ec, double* VA, double* EA, double* VB, double* EB, double* VC, double* EC){
  
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

  
  double* Delta=malloc(6*sizeof(double));
  Delta[0]=Pb-Pa; Delta[1]=Pc-Pa;
  Delta[2]=Tb-Ta; Delta[3]=Tc-Ta;
  Delta[4]=Gb-Ga; Delta[5]=Gc-Ga; 
  
  int dim=6;
  double *A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));
  // on met le second menbre dans la dernière colonne
  //matrice A
  A[dim*0+0]=dPva;   A[dim*0+1]=dPea;  A[dim*0+2]=-dPvb;
  A[dim*0+3]=-dPeb;  A[dim*0+4]=0.;    A[dim*0+5]=0.;  

  A[dim*1+0]=dPva;  A[dim*1+1]=dPea;   A[dim*1+2]=0.;
  A[dim*1+3]=0.;    A[dim*1+4]=-dPvc;  A[dim*1+5]=-dPec; 

  A[dim*2+0]=dTva;   A[dim*2+1]=dTea; 
  A[dim*2+2]=-dTvb;  A[dim*2+3]=-dTeb; 
  A[dim*2+4]=0.;     A[dim*2+5]=0.; 

  A[dim*3+0]=dTva;  A[dim*3+1]=dTea; 
  A[dim*3+2]=0.;    A[dim*3+3]=0.; 
  A[dim*3+4]=-dTvc; A[dim*3+5]=-dTec; 

  A[dim*4+0]=-Sa*dTva +dPva*Va;  A[dim*4+1]=-Sa*dTea +dPea*Va; 
  A[dim*4+2]=Sb*dTvb - dPvb*Vb;  A[dim*4+3]=Sa*dTea -dPea*Va; 
  A[dim*4+4]=0.;                 A[dim*4+5]=0.; 

  A[dim*5+0]=-Sa*dTva +dPva*Va;  A[dim*5+1]=-Sa*dTea +dPea*Va; 
  A[dim*5+2]=0.;                 A[dim*5+3]=0.; 
  A[dim*5+4]=Sc*dTvc -dPvc*Vc;   A[dim*5+5]=Sc*dTec -dPec*Vc;
  

  double deta=inverse_matrice_pivot(A, dim, InvA);
  if ( fabs(deta)<1e-8){printf("A non inversible : detA= %.15lf\n", deta);}
  
  /*
  // But résoudre Delta=A*dX
  // Inversion de la matrice A :
  double Ja, Jb, Jc, SV, detA;
  double Geaba, Gvaba, Gebba, Gvbba, Gecba, Gvcba;
  double Geaca, Gvaca, Gebca, Gvbca, Gecca, Gvcca;

  Ja=dPva*dTea - dPea*dTva;
  Jb=dPvb*dTeb - dPeb*dTvb;
  Jc=dPvc*dTec - dPec*dTvc;
  SV=Sc*(-Va+Vb) + Sb*(Va-Vc) + Sa*(-Vb+Vc);
  detA=Ja*Jb*Jc*SV;
  
  printf("  Sc= %.15lf, Sb= %.15lf, Sa= %.15lf\n", Sc,Sb,Sa);
  printf("  Va= %.15lf, Vb= %.15lf, Vc= %.15lf\n", Va,Vb,Vc);
  
  printf("  Ja= %.15lf, Jb= %.15lf, Jc= %.15lf, SV= %.15lf\n", Ja,Jb,Jc,SV);
  printf("  detA= %.15lf\n", detA);


  Geaba=dPea*(Vb-Va) - dTea*(Sb-Sa);
  Gvaba=dPva*(Vb-Va) - dTva*(Sb-Sa);
  Gebba=dPeb*(Vb-Va) - dTeb*(Sb-Sa);
  Gvbba=dPvb*(Vb-Va) - dTvb*(Sb-Sa);
  Gecba=dPec*(Vb-Va) - dTec*(Sb-Sa);
  Gvcba=dPvc*(Vb-Va) - dTvc*(Sb-Sa);
  Geaca=dPea*(Vc-Va) - dTea*(Sc-Sa);
  Gvaca=dPva*(Vc-Va) - dTva*(Sc-Sa);
  Gebca=dPeb*(Vc-Va) - dTeb*(Sc-Sa);
  Gvbca=dPvb*(Vc-Va) - dTvb*(Sc-Sa);
  Gecca=dPec*(Vc-Va) - dTec*(Sc-Sa);
  Gvcca=dPvc*(Vc-Va) - dTvc*(Sc-Sa);

  // Creation de la matrice InvA
  double** InvA=alloctabd(6,6);
  InvA[0][0]=-Geaca*Jb*Jc*Vb/detA; InvA[0][1]=Geaba*Jb*Jc*Vc/detA;  InvA[0][2]=Geaca*Jb*Jc*Sb/detA;
  InvA[0][3]=-Geaba*Jb*Jc*Sc/detA; InvA[0][4]=Geaca*Jb*Jc/detA;     InvA[0][5]=-Geaba*Jb*Jc/detA;  

  InvA[1][0]=Gvaca*Jb*Jc*Vb/detA;  InvA[1][1]=-Gvaba*Jb*Jc*Vc/detA; InvA[1][2]=-Gvaca*Jb*Jc*Sb/detA;
  InvA[1][3]=Gvaba*Jb*Jc*Sc/detA;  InvA[1][4]=-Gvaca*Jb*Jc/detA;    InvA[1][5]=Gvaba*Jb*Jc/detA; 

  InvA[2][0]=-Ja*Jc*(SV*dTeb*Gebca*Vb)/detA; InvA[2][1]=Gebba*Ja*Jc*Vc/detA; 
  InvA[2][2]=Ja*Jc*(Gebca*Sb+dPeb*SV)/detA;  InvA[2][3]=-Gebba*Ja*Jc*Sc/detA; 
  InvA[2][4]=Gebca*Ja*Jc/detA;               InvA[2][5]=-Gebba*Ja*Jc/detA; 

  InvA[3][0]=Ja*Jc*(SV*dTvb+Gvbca*Vb)/detA;  InvA[3][1]=-Gvbba*Ja*Jc*Vc/detA; 
  InvA[3][2]=-Ja*Jc*(Gvbca*Sb+dPvb*SV)/detA; InvA[3][3]=Gvbba*Ja*Jc*Sc/detA; 
  InvA[3][4]=-Gvbca*Ja*Jc/detA;              InvA[3][5]=Gvbba*Ja*Jc/detA; 

  InvA[4][0]=-Gecca*Ja*Jb*Vb/detA; InvA[4][1]=-Ja*Jb*(SV*dTec*Gecba*Vc)/detA; 
  InvA[4][2]=Gecca*Ja*Jb*Sb/detA;  InvA[4][3]=Ja*Jb*(-Gecba*Sc+dPec*SV)/detA; 
  InvA[4][4]=Gecca*Ja*Jb/detA;     InvA[4][5]=-Gecba*Ja*Jb/detA; 

  InvA[5][0]=Gvcca*Ja*Jb*Vb/detA;  InvA[5][1]=Ja*Jb*(SV*dTvc-Gvcba*Vc)/detA; 
  InvA[5][2]=-Gvcca*Ja*Jb*Sb/detA; InvA[5][3]=-Ja*Jb*(-Gvcba*Sc + dPvc*SV)/detA; 
  InvA[5][4]=-Gvcca*Ja*Jb/detA;    InvA[5][5]=Gvcba*Ja*Jb/detA;
  */


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



  // Mutiplcation matricielle
  double* dX=malloc(6*sizeof(double));

  for(int i=0; i<6; i++){
    dX[i]=0;
    for(int j=0; j<6; j++){
      dX[i]+=InvA[dim*i+j]*Delta[j];
    }
    printf("  dX[%d]= %.15lf\n",i,dX[i] );
  }
  printf("\n");
   
  double correc=0.15;

  *VA=Va+correc*dX[0];
  *EA=Ea+correc*dX[1];
  *VB=Vb+correc*dX[2];
  *EB=Eb+correc*dX[3];
  *VC=Vc+correc*dX[4];
  *EC=Ec+correc*dX[5];



  //freetab(InvA);
  free(Delta);
  free(dX);
}






void triple(int Nmax, double epsilon, double* Va0, double* Ea0, double* Vb0, double* Eb0, double* Vc0, double* Ec0, double* Gares, double* Gbres, double* Gcres){
  int n=0;

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

  ///////////////////////////////////
  printf("triple init :\n");
  printf(" Sa= %.15lf, Sb= %.15lf,  Sc= %.15lf\n",Sa, Sb, Sc );
  printf(" Ga= %.15lf, Gb= %.15lf,  Gc= %.15lf\n",Ga, Gb, Gc );
  printf(" Pa= %.15lf, Pb= %.15lf,  Pc= %.15lf  (GPa)\n",Pa*1e-9, Pb*1e-9, Pc*1e-9 );
  printf(" Ta= %.15lf, Tb= %.15lf,  Tc= %.15lf\n",Ta, Tb, Tc );
  printf(" critere= %.15lf \n", fabs(Ga-Gb)+fabs(Ga-Gc));
  //////////////////////////////////

  double critere=fabs(Ga-Gb)+fabs(Ga-Gc);

  while(critere>=epsilon && n<Nmax){
    printf("n= %d\n",n );
    Newton_triple(Va, Ea, Vb, Eb, Vc, Ec, &VA, &EA, &VB, &EB, &VC, &EC);

    printf(" dVa= %.15lf, VA= %.15lf (cm^3/kg)\n",(VA-Va)*1e6, VA*1e6 );
    printf(" dEa= %.15lf, EA= %.15lf (kJ/kg)\n",(EA-Ea)*1e-3, EA*1e-3 );
    printf(" dVb= %.15lf, VB= %.15lf (cm^3/kg)\n",(VB-Vb)*1e6, VB*1e6 );
    printf(" dEb= %.15lf, EB= %.15lf (kJ/kg)\n",(EB-Eb)*1e-3, EB*1e-3 );
    printf(" dVc= %.15lf, VC= %.15lf (cm^3/kg)\n",(VC-Vc)*1e6, VC*1e6 );
    printf(" dEc= %.15lf, EC= %.15lf (kJ/kg)\n",(EC-Ec)*1e-3, EC *1e-3);

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
 
 /////////////////////////////////////////////////////////////////////
    printf(" Iteration (triple) :\n");
    printf("  Sa= %.15lf, Sb= %.15lf,  Sc= %.15lf\n",Sa, Sb, Sc );
    printf("  Ga= %.15lf, Gb= %.15lf,  Gc= %.15lf\n",Ga, Gb, Gc );
    printf("  Pa= %.15lf, Pb= %.15lf,  Pc= %.15lf  (GPa)\n",Pa*1e-9, Pb*1e-9, Pc*1e-9 );
    printf("  Ta= %.15lf, Tb= %.15lf,  Tc= %.15lf\n",Ta, Tb, Tc );
/////////////////////////////////////////////////////////////////////


    critere=fabs(Ga-Gb)+fabs(Ga-Gc);
    printf(" critere= %.14lf\n",critere );
    printf(" Ga= %.14lf, Gc= %.14lf, Gc= %.14lf\n",Ga,Gb,Gc );


    n++;
  }


  *Va0=Va; *Ea0=Ea; 
  *Vb0=Vb; *Eb0=Eb; 
  *Vc0=Vc; *Ec0=Ec; 

  *Gares=Ga; *Gbres=Gb; *Gcres=Gc;

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


int main(){
  int Nv=200;
  int Ne=200;

  double Vmin=90e-6;
  double Vmax=150e-6;
  double Emin=50e3;
  double Emax=200e3;
  
  /*
  int res=trace_phase(Nv, Ne, Vmin, Vmax, Emin, Emax);
  if (res){
    return res;
  }
  */
  
  
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
*/

  /*
  int phase=3;
  V=130e-6;  
  double G;
  for(int i=0; i<Ne; i++){
    E=Emin+i*he;
    G=E+fP(V,E,phase)*V-fS(V,E,phase)*fT(V,E,phase);
    printf("phase= %d, E= %.15lf, T= %.15lf, P= %.15lf, G= %.15lf, S= %.15lf\n",phase, E ,fT(V,E,phase), fP(V,E,phase)*1e-9,G,fS(V,E,phase) );
  }
  
*/



  int Nmax=1e2;
  double epsilon=1e-8;

  double Gares, Gbres, Gcres;

  double Pa, Ta, Sa, Ga, Ea, Va;
  double Pb, Tb, Sb, Gb, Eb, Vb;
  double Pc, Tc, Sc, Gc, Ec, Vc;

  double Va0=129e-6;
  double Ea0=100e3;

  double Vb0=0.98*Va0;
  double Eb0=1.01*Ea0;
  
  double Vc0=1.02*Va0;
  double Ec0=1.05*Ea0;
  

/*
  double Vtest=110e-6; //97.41757866e-6;
  double Etest=300e3; //157.45449603e3; 
  double Ptest=fP(Vtest,Etest, 3); double Ttest=fT(Vtest,Etest, 3);
  double Stest=fS(Vtest,Etest, 3); double Gtest=Etest+Ptest*Vtest-Ttest*Stest;

  printf("main test phase 3:\n");
  printf(" Vtest= %.10lf (cm^3/kg), Etest= %.10lf (kJ/kg)\n",Vtest*1e6, Etest*1e-3 );
  printf(" Ptest= %.10lf, Ttest= %.10lf, Stest= %.10lf, Gtest= %.10lf\n",Ptest*1e-9, Ttest, Stest, Gtest );
*/

 
//////////////////////////////
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

  printf("main init :\n");
  printf(" Sa= %.15lf, Sb= %.15lf,  Sc= %.15lf\n",Sa, Sb, Sc );
  printf(" Ga= %.15lf, Gb= %.15lf,  Gc= %.15lf\n",Ga, Gb, Gc );
  printf(" Pa= %.15lf, Pb= %.15lf,  Pc= %.15lf  (GPa)\n",Pa*1e-9, Pb*1e-9, Pc*1e-9 );
  printf(" Ta= %.15lf, Tb= %.15lf,  Tc= %.15lf\n",Ta, Tb, Tc );
////////////////////////////////////////
  

  
  triple(Nmax, epsilon, &Va0, &Ea0, &Vb0, &Eb0, &Vc0, &Ec0, &Gares, &Gbres, &Gcres);
  

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
  printf(" Ga= %.15lf, Gb= %.15lf,  Gc= %.15lf\n",Gares, Gbres, Gcres );
  printf(" Pa= %.15lf, Pb= %.15lf,  Pc= %.15lf  (GPa)\n",Pa*1e-9, Pb*1e-9, Pc*1e-9 );
  printf(" Ta= %.15lf, Tb= %.15lf,  Tc= %.15lf\n",Ta, Tb, Tc );
  

}

