#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "head.h"
#include "head_etain.h"



// SCHEMA vNR  JCP 2009

// données 
//   fonction ENERGIE INTERNE SPECIFIQUE epsilon(tau, p)
double fPn1(int tst, double TAUn, double TAUn1, double Pn, double Qn, double EPSn){
  double res;
  double dTAU=TAUn1-TAUn; 

  if (tst==0 || tst==3 || tst==4|| tst==5|| tst==6){  
    double Gamma=1.4;

    res=( EPSn - dTAU*(Pn/2 + Qn) )/( TAUn1/(Gamma-1) + dTAU/2 ) ;
  }
  else if (tst==2) {  // Bizarrium
    double rho0 = 1e4; double tau0=1./1e4;
    double K0 = 1e11;
    double Cv0 = 1e3;
    double T0 = 300;
    double eps0 = 0;
    double sigma0 = 1.5;
    double s = 1.5;
    double q = -(42080895.0)/(14941154.0);
    double r = (727668333.0)/(149411540.0);
    double x= tau0/TAUn1 - 1.0;

    double gEOS= sigma0*(1.0 - TAUn1/tau0);

    double f0 = 1.0 + (s/3.0 - 2.0)*x + q*x*x + r*x*x*x;
    f0=f0/(1.0-s*x);

    double f1 = s/3.0 - 2.0 + 2.0*q*x + 3.0*r*x*x + s*f0;
    f1=f1/(1.0-s*x);

    double epsilon_k0= eps0;
    epsilon_k0-=Cv0*T0*(1.0+ gEOS);
    epsilon_k0+=K0*tau0*x*x*f0/2.0;

    double p_k0= -Cv0*T0*sigma0*rho0;
    p_k0+= K0*x*(1.0+x)*(1.0+x)*(2.0*f0 + x*f1)/2.0;

    // eps_n = epsilon_k0 + (tau0/sigma0)*(p - p_k0);

    res=( -epsilon_k0 + tau0*p_k0/sigma0 + EPSn - (Pn/2 + Qn)*dTAU )/( tau0/sigma0 + dTAU/2 );
  }
  else if (tst==1){  // LeBlanc
    double Gamma=5./3;

    res=( EPSn - dTAU*(Pn/2 + Qn) )/( TAUn1/(Gamma-1) + dTAU/2 ) ;
  }
  else { printf("Erreur choix tst dans vNR/fPn1\n");}
  return res;
}




/////////////////////  MAIN  ////////////////////
int funcvNR_etain(int tst, double Tf, double a, double b, int nx, int Nmax, double CFL, int aff, int ind_cond, int iu, int Ntau, double dt_first, double Ttau, int cin_phase, int sch_cin_phase, int tst_multi, int q, double Cq, double Cl,
                  int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, int** tabZONE, int* tab_Nb_ZONE, 
                  double* tabVx, double* tabEx, double* tabVy, double* tabEy,
                  double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC){


  int etain=tst/100;


	// variables
	double t=0;
	double dt,dt1,dt_new;
	int n=0,err=0;
  double maxRCsM, maxDTAUsRT; 
	char* W0;
  double sum1, sum2, sum3, sum4, sum5;
  int N;
  int nbghosts=10; // nombre de mailles fantomes pour les variables du schéma (\rho, \tau, u, \epsilon, ekin)
  int nbghosts2=4; 
  // l'idee  : appliquer les cond aux lim sur les var. sur nbghosts mailles qui calculer la valeur des fcts sur les nbghosts2 mailles
  int ideb, ifin;
  double dmm, dmp;
  double TAUtemp;
  double Pn1, EPSn1;
  double dX;
  int zone;

  double dx=(b-a)/nx; // pas d'espace

  // Calcul de la taille des tableaux
  N=nx+2*nbghosts;
  ideb=nbghosts;
  ifin=nbghosts+nx-1;
  int N2=nbghosts+nx+nbghosts2; // cond au lim de Q et P sur i=nbghosts- nbghosts2, N2

  // Allocation des tableaux
  double* X=malloc((N+1)*sizeof(double)); // maillage
  double* Xc=malloc((N)*sizeof(double)); // maillage
  // PW
	double* RHO0c=malloc((N)*sizeof(double));   //  rho0  (centrées aux mailles)
	double* RHO0d=malloc((N+1)*sizeof(double)); //  rho0  (dé-centrées | frontières)
   
	// solution 
	double* TAU=malloc((N)*sizeof(double)); //  tau  (centrées aux mailles)
	double* U=malloc((N+1)*sizeof(double)); //  u (frontières)
  double* EPS=malloc((N)*sizeof(double)); //  epsilon (centrées aux mailles)
  double* E=malloc((N)*sizeof(double)); //  energie totale (centrées aux mailles)

  // solution intermedaire
  double* TAU1=malloc((N)*sizeof(double)); //  tau au temsp n+1/2 (centrées aux mailles)

  double* P=malloc((N)*sizeof(double));
  double* Q=malloc((N)*sizeof(double));   // artificial viscosity
  double* Qt=malloc((N)*sizeof(double));   // artificial viscosity
  double* C=malloc((N)*sizeof(double));

  double* P1=malloc((N)*sizeof(double));
  double* C1=malloc((N)*sizeof(double));

  double* RCsM=malloc((N)*sizeof(double));  // rho.c/Dm   centres aux mailles
  double* DTAUsRT=malloc((N)*sizeof(double));  // dtau/rho.dt   centres aux mailles

  double* DMc=malloc((N)*sizeof(double));  // Delta m = (X[i+1]-X[i])*(RHOC[i]/RTAU[i])   centres aux mailles
  double* DMd=malloc((N+1)*sizeof(double));  // Delta m = (Dm + Dm)/2   dé-centres aux mailles
  
  double* T=malloc((N)*sizeof(double));  // température 
  double* FM=malloc((N)*sizeof(double)); // fraction massique
  
  double** matFM=alloctabd(N,3);          // fraction massique
  double* lambda_out_unuse=malloc(3*sizeof(double));
  
	// Maillage
	for(int i=0; i<=N; i++){
	  X[i]=(a-nbghosts*dx)+dx*i;
    //printf("i= %d, X=%lf\n",i,X[i]);
	}
  
  double* W=malloc(4*sizeof(double));
  
  // Initialisation 
  double x;
  // variables défini aux centres des mailles
  if(etain){
  	for(int i=0; i<=N-1; i++){
  	  x=(X[i]+X[i+1])/2.0;
      err=u0_etain(tst,a,b, x, W);
      if(err){return err;}
      RHO0c[i]=W[0];
      TAU[i]=1./W[0];
      T[i]=W[3];
      //EPS[i]=epsilonEOS(tst, 1./W[0], W[2]);
      err=epsilon_etain(1./W[0], W[2], T[i], &EPS[i], nb_points, SA, SB, SC, SAB, SAC, SBC );
      if(err){return err;}
  	  //printf("i= %d, tau=%f, eps=%f, rho0c=%f\n",i,TAU[i],EPS[i],RHO0c[i] );
  	}
    // variables défini aux frontières des mailles
    for(int i=0; i<=N; i++){
      x=X[i];
      err=u0_etain(tst,a,b, x, W);
      if(err){return err;}
      RHO0d[i]=W[0];
      U[i]=W[1];
      //printf("u=%.8lf\n",W[1] );
      //printf("i= %d,  rho0.u=%f,  rho0d=%f\n",i,RUpw[i],RHO0d[i] );
    }
  }
  else{
    for(int i=0; i<=N-1; i++){
      x=(X[i]+X[i+1])/2.0;
      err=u0(tst,a,b, x, W);
      if(err){return err;}
      RHO0c[i]=W[0];
      TAU[i]=1./W[0];
      EPS[i]=epsilonEOS(tst, 1./W[0], W[2]);
      //printf("i= %d, tau=%f, eps=%f, rho0c=%f\n",i,TAU[i],EPS[i],RHO0c[i] );
    }
    // variables défini aux frontières des mailles
    for(int i=0; i<=N; i++){
      x=X[i];
      err=u0(tst,a,b, x, W);
      if(err){return err;}
      RHO0d[i]=W[0];
      U[i]=W[1];
      //printf("u=%.8lf\n",W[1] );
      //printf("i= %d,  rho0.u=%f,  rho0d=%f\n",i,RUpw[i],RHO0d[i] );
    }
  }


  //Init des tableaux 
  if(etain){ 
    for(int i=0; i<=N-1; i++){
      //EGS(tst,TAU[i], EPS[i], &P[i], &C[i]);
      /*
      err=fPTC(TAU[i], EPS[i], &P[i], &C[i], &T[i], &FM[i],
               nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, 
               SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
      */
      err=fPTC_cin_phase(0, 0, sch_cin_phase, TAU[i], EPS[i], &P[i], &C[i], &T[i], matFM[i], matFM[i], Ntau, dx, dt, Ttau, &zone,
                   nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, 
                   SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
      if(err){return err;}
      DMc[i]=(X[i+1]-X[i])*RHO0c[i];
      RCsM[i]=C[i]/(DMc[i]*TAU[i]);
      //printf("i=%d P= %lf, C=%lf, U=%lf, EPS=%lf, XsC=%lf \n",i,P[i],C[i],U[i],REPS[i]/RHO0c[i], XsC[i]);           
    }
  }
  else{
    for(int i=0; i<=N-1; i++){
      EGS(tst,TAU[i], EPS[i], &P[i], &C[i]);
      DMc[i]=(X[i+1]-X[i])*RHO0c[i];
      RCsM[i]=C[i]/(DMc[i]*TAU[i]);
      //printf("i=%d P= %lf, C=%lf, U=%lf, EPS=%lf, XsC=%lf \n",i,P[i],C[i],U[i],REPS[i]/RHO0c[i], XsC[i]);           
    }
  }

  for(int i=1; i<=N; i++){
    DMd[i]=(DMc[i-1]+DMc[i])/2;          
  }


	// maximum de DMsRC  (calcul de dt)
  maxRCsM=RCsM[ideb];
  for(int j=ideb+1; j<=ifin; j++){
    if(maxRCsM<RCsM[j]){
      maxRCsM=RCsM[j];
    }
  }

  dt = dt_first*CFL /maxRCsM; // pas de temps
	printf("nx= %d et dt= %g\n", nx, dt);
  
  // energie totale du système
  int nbiter=5*ceil(Tf/dt);
  //printf("nbiter :%d\n",nbiter);
  double* Etot=malloc(nbiter*sizeof(double));
  double* IMPUL=malloc(nbiter*sizeof(double)); // impulsion
  double* Mtot=malloc(nbiter*sizeof(double));  // masse totale
  double* VOL=malloc(nbiter*sizeof(double));   // volume total
  

  sum1=0; sum3=0.; sum4=0.; sum5=0.;
  for(int i=ideb; i<=ifin; i++){
    sum1+=dx*RHO0c[i]*EPS[i]+dx*RHO0d[i]*U[i]*U[i]/2;;
    sum3+=dx*RHO0c[i]*U[i];
    sum4+=dx*RHO0c[i];
    sum5+=dx*RHO0c[i]*TAU[i];
  }
  sum1+=dx*RHO0d[ifin+1]*U[ifin+1]*U[ifin+1]/2;
  sum3+=dx*RHO0c[ifin+1]*U[ifin+1];
  Etot[0]=sum1;
  IMPUL[0]=sum3;
  Mtot[0]=sum4;
  VOL[0]=sum5;

  /*
  for (int i = ideb; i <=ifin; ++i){
    printf("i=%d, RTAU=%lf, RU=%lf, REPS=%lf\n",i,RTAU[i],RU[i],REPS[i] );
  }
  */

// debut de la boucle =============================================================================
	while((t<Tf)&&(n<Nmax)){
    // cond aux limites
    err=condlim(ind_cond, ideb,ifin, nbghosts, TAU, 0);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin, nbghosts, EPS, 0);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin+1, nbghosts, U, iu);
    if(err){return err;}


    // 1 - Computation of the velocity at time n+1/2 and node positions
    for(int i=nbghosts-nbghosts2; i<=N2; i++){
      U[i]-= (dt/DMd[i])*( (P[i]+Q[i]) - (P[i-1]+Q[i-1]) );
    }
    err=condlim(ind_cond, ideb, ifin+1, nbghosts, U, iu);
    if(err){return err;}

    for(int i=0; i<=N; i++){
      X[i]+= dt*U[i];
    }

    // 2 - Computation of the density (and specific volume) at time n+1
    for(int i=ideb; i<=ifin; i++){
      TAU1[i]=(X[i+1]-X[i])/DMc[i];
    }
    err=condlim(ind_cond, ideb, ifin, nbghosts, TAU1, 0);
    if(err){return err;}

    // 3 - Approximation of the artificial viscosity at time n+1/2
    for(int i=nbghosts-nbghosts2; i<=N2-1; i++){
      Qt[i]=qvis(q,Cq,Cl, &U[i], (TAU1[i]+TAU[i])/2, C[i], &DMc[i]);
    }

    // 4 - Computation of the internal energy at time n+1
    if(etain){
      for(int i=ideb; i<=ifin; i++){
        //printf("-i= %d   ifin= %d\n",i,ifin );
        //Pn1=fPn1(tst, TAU[i], TAU1[i], P[i], Qt[i], EPS[i]);
        /*
        err= fPn1_etain(tst, TAU[i], TAU1[i], P[i], Qt[i], EPS[i], &EPSn1, &Pn1,
                 nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, 
                 SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
        */
        err= fPn1_etain_cin_phase(tst, cin_phase, sch_cin_phase, TAU[i], TAU1[i], P[i], Qt[i], EPS[i], &EPSn1, &Pn1, matFM[i], Ntau, dX, dt, Ttau,
                                  nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, 
                                  SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
        if(err){return err;}
        //EPS[i]-= ((Pn1+P[i])/2 + Qt[i])*(TAU1[i]-TAU[i]);
        EPS[i]= EPSn1;

      }
    }
    else{
      for(int i=ideb; i<=ifin; i++){
        Pn1=fPn1(tst, TAU[i], TAU1[i], P[i], Qt[i], EPS[i]);
        EPS[i]-= ((Pn1+P[i])/2 + Qt[i])*(TAU1[i]-TAU[i]);
      }
    }


    err=condlim(ind_cond, ideb, ifin, nbghosts, EPS, 0);
    if(err){return err;}
    
    if(etain){
      for(int i=0; i<=N-1; i++){
        dX=X[i+1]-X[i];
        //EGS(tst, TAU1[i], EPS[i], &P1[i], &C1[i]);
        /*
        err=fPTC(TAU1[i], EPS[i], &P1[i], &C1[i], &T[i], &FM[i],
                 nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, 
                 SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
        */
        err=fPTC_cin_phase(0, cin_phase, sch_cin_phase, TAU1[i], EPS[i], &P1[i], &C1[i], &T[i], matFM[i], matFM[i], Ntau, dX, dt, Ttau, &zone,
                   nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, 
                   SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
        if(err){return err;}
      }
    }
    else{
      for(int i=0; i<=N-1; i++){
        EGS(tst, TAU1[i], EPS[i], &P1[i], &C1[i]);
      }
    }

    // 5 - Computation of the artificial viscosity at time n+1/2
    for(int i=nbghosts-nbghosts2; i<=N2-1; i++){
      Q[i]=qvis(q,Cq,Cl, &U[i], TAU1[i], C1[i], &DMc[i]);
    }

    // energie totale


    ///////////////
    // Energie du système
    sum1=0; sum3=0.; sum4=0.; sum5=0.;
    for(int i=ideb; i<=ifin; i++){
      sum1+=dx*RHO0c[i]*EPS[i]+dx*RHO0d[i]*U[i]*U[i]/2;
      sum3+=(X[i+1]-X[i])*U[i]/TAU[i];
      sum4+=(X[i+1]-X[i])/TAU[i];
      sum5+=(X[i+1]-X[i])*RHO0c[i]*TAU[i];
    }
    sum1+=dx*RHO0d[ifin+1]*U[ifin+1]*U[ifin+1]/2;
    sum3+=(X[ifin+1]-X[ifin])*U[ifin+1]/TAU[ifin+1];
    Etot[n+1]=sum1;
    IMPUL[n+1]=sum3;
    Mtot[n+1]=sum4;
    VOL[n+1]=sum5;

    /*
    for (int i = ideb; i <= ifin; i++){
      printf("i=%d, TAU_n-1= %lf, TAU= %lf, U=%lf, EPS=%lf\n",i,TAU[i],TAU1[i],U[i],EPS[i] );
    }
    printf("\n");
    */


    // Mise a jour de RCsM
    for(int i=ideb; i<=ifin; i++){
      RCsM[i]=C1[i]/(DMc[i]*TAU1[i]);
      DTAUsRT[i]=(1./TAU1[i] - 1./TAU[i])*TAU1[i]/(dt);
    }

    // Calcul du pas de temps
    maxRCsM=RCsM[ideb];
    maxDTAUsRT=DTAUsRT[ideb];
    for(int j=ideb+1; j<=ifin; j++){
      if (maxRCsM<RCsM[j]){
        maxRCsM=RCsM[j];
      }
      if (maxDTAUsRT<DTAUsRT[j]){
        maxDTAUsRT=DTAUsRT[j];
      }
    }
    dt1 = CFL /fmax(maxRCsM, maxDTAUsRT) ; // pas de temps
    dt_new=fmin(dt1, (dt1+dt)/2);
    //if ((T-t)<dt){ dt=T-t; }
    dt=pastemps(Tf,t,dt_new,dt);

    // On met a jour la solution
    for (int i=0; i<=N-1; i++){
      TAU[i]=TAU1[i];
      P[i]=P1[i];
      C[i]=C1[i];
    }


	  t=t+dt;
	  n=n+1;
    if(aff==1){
	    printf("n= %d, t= %g, dt= %g\n",n,t,dt);
    }
	}
//==========================================================================

  printf("Nombre d'itérations : %d\n",n);
  printf("Temps final : %g\n",t);

  // energie totale
  for(int i=ideb; i<=ifin; i++){
    dmm=(DMc[i-1]+DMc[i])*(U[i]*U[i])/DMc[i];
    dmp=(DMc[i+1]+DMc[i])*(U[i+1]*U[i+1])/DMc[i];
    E[i]=EPS[i] + (1./8)*(dmm+dmp);
    Xc[i]=(X[i+1]+X[i])/2;
  }


    // printf_sol
  err=print_sol(0, ideb, ifin, X, Xc, TAU, U, E, P, EPS);
  if(err){return err;}

  err=print_nrj(n, Etot, IMPUL, Mtot, VOL);
  if(err){return err;}


  if(etain){
    FILE* Res;
    /*
    L1 : xA
    L2 : xB
    L3 : xC
    */
    if((Res = fopen("frac_mass.txt", "w+")) != NULL){
      // xA
      for(int i=ideb; i<=ifin; i++){  fprintf(Res,"%.15lf ",matFM[i][0]);  }
      fprintf(Res, "1e20\n \n");
      // xB
      for(int i=ideb; i<=ifin; i++){  fprintf(Res,"%.15lf ",matFM[i][1]);  }
      fprintf(Res, "1e20\n \n");
      // xC
      for(int i=ideb; i<=ifin; i++){  fprintf(Res,"%.15lf ",matFM[i][2]);  }
      fprintf(Res, "1e20\n \n");
      fclose(Res);
      return 0;
    } 
    else{printf("erreur impression\n"); return 1;}
  }



  free(W); free(lambda_out_unuse);

  free(X); free(Xc);
  free(RHO0d); free(RHO0c);
  free(DMc); free(DMd);
  
  free(TAU); free(U); free(EPS);

  free(P); free(Q); free(C); 
  free(RCsM); free(DTAUsRT);
  free(E); 

  free(P1); free(C1); free(Qt);
  free(TAU1); 
  
  //  
  free(Etot); 
  free(IMPUL);
  free(Mtot);
  free(VOL);

  return 0;
  
}
