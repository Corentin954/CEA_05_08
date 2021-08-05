#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "head.h"
#include "head_etain.h"



// SCHEMA CAUCHY-KOVALESKAYA FORMULE EN ENERGIE INTERNE

/////////////////////  MAIN  ////////////////////
int funcKOVA_etain(int tst, double Tf, double a, double b, int nx, int Nmax, double CFL, int aff, int ind_cond, int iu, int Ntau, double dt_first, double Ttau, int cin_phase, int sch_cin_phase, int tst_multi, int z,
                 int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, int** tabZONE, int* tab_Nb_ZONE, 
                   double* tabVx, double* tabEx, double* tabVy, double* tabEy,
                   double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC){
 
  int etain=tst/100;

	// variables
  double xi=1./4;
	double t=0;
	double dt,dt_new;
	int n=0,err=0;
  double minXsC; 
	char* W0;
  double sum1, sum3, sum4, sum5;
  int N;
  int nbghosts=10; // nombre de mailles fantomes pour les variables du schéma (\rho, \tau, u, \epsilon, ekin)
  int nbghosts2=6; // nombre de mailles fantomes pour les fonctions qui dependent des variables ci-dessus
  // l'idee  : appliquer les ocnd aux lim sur les var. sur nbghosts mailles qui calculer la valeur des fcts sur les nbghosts2 mailles
  int ideb, ifin;
  double rhoc;
  double dX;


  double dx=(b-a)/nx; // pas d'espace


  // Calcul de la taille des tableaux
  N=nx+2*nbghosts;
  ideb=nbghosts;
  ifin=nbghosts+nx-1;
  int N2=nbghosts+nx+nbghosts2; // cond au lim de Q et P sur i=nbghosts- nbghosts2, N2
  int ideb2=nbghosts-nbghosts2;

  // Allocation des tableaux
  double* X=malloc((N+1)*sizeof(double)); // maillage
  double* Xc=malloc((N)*sizeof(double)); // maillage
  double* E=malloc((N)*sizeof(double));
  // PW
	double* RHO0c=malloc((N)*sizeof(double));   //  rho0  (centrées aux mailles)
	double* RHO0d=malloc((N+1)*sizeof(double)); //  rho0  (dé-centrées | frontières)
   
	// solution 
	double* RTAU=malloc((N)*sizeof(double)); //  rho0.tau  (centrées aux mailles)
	double* RU=malloc((N+1)*sizeof(double)); //  rho0.u (frontières)
  double* REPS=malloc((N)*sizeof(double)); //  rho0.epsilon (centrées aux mailles) 


  // AVERAGE (plot)
  double* TAU=malloc((N)*sizeof(double)); //  tau  (centrées aux mailles)
  double* EPS=malloc((N)*sizeof(double)); //  epsilon (centrées aux mailles)  
  double* U=malloc((N+1)*sizeof(double)); //  u  (frontières aux mailles)
  double* R2C2dU=malloc((N)*sizeof(double)); //  epsilon (centrées aux mailles)  

  // Point-Wise
  //double* TAU=malloc((N)*sizeof(double)); //  tau  (centrées aux mailles)
  //double* EPS=malloc((N)*sizeof(double)); //  epsilon (centrées aux mailles)
  double* P=malloc((N)*sizeof(double));
  double* C=malloc((N)*sizeof(double));
  double* XsC=malloc((N)*sizeof(double));  // dX/c   centres aux mailles

  
  // tableaux etoiles
  double* Ustar=malloc((N+1)*sizeof(double));
  double* Pstar=malloc((N)*sizeof(double));

  // kinetic fix  
  double* REKIN=malloc((N+1)*sizeof(double));
  double* DK=malloc((N+1)*sizeof(double));
  double* R2U2=malloc((N+1)*sizeof(double)); //(point-wise)

  double* RU2=malloc((N+1)*sizeof(double)); //(point-wise)

  double* T=malloc((N)*sizeof(double));  // température 
  double* FM=malloc((N)*sizeof(double)); // fraction massique


	// Maillage
	for(int i=0; i<=N; i++){
	  X[i]=(a-nbghosts*dx)+dx*i;
    //printf("i= %d, X=%lf\n",i,X[i]);
	}
  
  double* W=malloc(4*sizeof(double));
  // Initialisation POINT-WISE
  double x, EPStemp;
  // variables défini aux centres des mailles
  if(etain){
  	for(int i=0; i<=N-1; i++){
      RTAU[i]=1.0;
  	  x=(X[i]+X[i+1])/2.0;
      err=u0_etain(tst,a,b, x, W);
      if(err){return err;}
      RHO0c[i]=W[0];
      //REPS[i]=W[0]*epsilonEOS(tst, 1./W[0], W[2]);
      err=epsilon_etain(1./W[0], W[2], T[i], &EPStemp, nb_points, SA, SB, SC, SAB, SAC, SBC );
      if(err){return err;}
      REPS[i]=W[0]*EPStemp;
      //printf("rho=%lf, p=%lf\n",W[0], W[2] );
  	  //printf("i= %d, rho0.tau=%f, rho0.eps=%f, rho0c=%f\n",i,RTAUpw[i],REPSpw[i],RHO0c[i] );
  	}
    // variables défini aux frontières des mailles
    for(int i=0; i<=N; i++){
      x=X[i];
      err=u0_etain(tst,a,b, x, W);
      if(err){return err;}
      RHO0d[i]=W[0];
      RU[i]=W[0]*W[1];
      REKIN[i]=W[0]*W[1]*W[1]/2.0;
      //printf("u=%.8lf\n",W[1] );
      //printf("i= %d,  rho0.u=%f,  rho0d=%f\n",i,RU[i],RHO0d[i] );
    }
  }
  else{
    for(int i=0; i<=N-1; i++){
      RTAU[i]=1.0;
      x=(X[i]+X[i+1])/2.0;
      err=u0(tst,a,b, x, W);
      if(err){return err;}
      RHO0c[i]=W[0];
      REPS[i]=W[0]*epsilonEOS(tst, 1./W[0], W[2]);
      //printf("rho=%lf, p=%lf\n",W[0], W[2] );
      //printf("i= %d, rho0.tau=%f, rho0.eps=%f, rho0c=%f\n",i,RTAUpw[i],REPSpw[i],RHO0c[i] );
    }
    // variables défini aux frontières des mailles
    for(int i=0; i<=N; i++){
      x=X[i];
      err=u0(tst,a,b, x, W);
      if(err){return err;}
      RHO0d[i]=W[0];
      RU[i]=W[0]*W[1];
      REKIN[i]=W[0]*W[1]*W[1]/2.0;
      //printf("u=%.8lf\n",W[1] );
      //printf("i= %d,  rho0.u=%f,  rho0d=%f\n",i,RU[i],RHO0d[i] );
    }
  }



  //Init des tableaux 
  if(etain){ 
    for(int i=0; i<=N-1; i++){
      //EGS(tst,RTAU[i]/RHO0c[i], REPS[i]/RHO0c[i], &P[i], &C[i]);
      /*
      err=fPTC(RTAU[i]/RHO0c[i], REPS[i]/RHO0c[i], &P[i], &C[i], &T[i], &FM[i],
               nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, 
               SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
      */
      if(err){return err;}
      XsC[i]=(X[i+1]-X[i])/C[i];
      //printf("i=%d P= %lf, C=%lf, U=%lf, EPS=%lf, XsC=%lf \n",i,P[i],C[i],U[i],REPS[i]/RHO0c[i], XsC[i]);           
    }
  }
  else{
    for(int i=0; i<=N-1; i++){
      EGS(tst,RTAU[i]/RHO0c[i], REPS[i]/RHO0c[i], &P[i], &C[i]);
      XsC[i]=(X[i+1]-X[i])/C[i];
      //printf("i=%d P= %lf, C=%lf, U=%lf, EPS=%lf, XsC=%lf \n",i,P[i],C[i],U[i],REPS[i]/RHO0c[i], XsC[i]);           
    }
  }


  /*
  printf("Initialisation :\n");
  for (int i = ideb; i <= ifin; i++){
    printf("i=%d, Rtau= %lf, RU=%lf, REPS=%lf REKIN= %lf\n",i,RTAU[i],RU[i],REPS[i],REKIN[i] );
  } 
  printf("\n");
  */


	// minimum de XsC  (calcul de dt)
  minXsC=XsC[ideb];
  for(int j=ideb+1; j<=ifin; j++){
    if(minXsC>XsC[j]){
      minXsC=XsC[j];
    }
  }

  printf("minXsC= %g\n",minXsC );
  dt = dt_first*CFL * minXsC; // pas de temps
	printf("nx= %d et dt= %g\n", nx, dt);
  
  // energie totale du système
  int nbiter=5*ceil(Tf/dt);
  printf("nbiter :%d\n",nbiter);
  double* Etot=malloc(nbiter*sizeof(double));
  double* IMPUL=malloc(nbiter*sizeof(double)); // impulsion
  double* Mtot=malloc(nbiter*sizeof(double)); // masse totale
  double* VOL=malloc(nbiter*sizeof(double)); // volume total
  
  // Energie du système

  sum1=0; sum3=0.; sum4=0.; sum5=0.;
  for(int i=ideb; i<=ifin; i++){
    sum1+=dx*(REPS[i]+REKIN[i]);
    sum3+=dx*RU[i];
    sum4+=dx*RHO0c[i];
    sum5+=dx*RTAU[i];
  }
  sum1+=dx*REKIN[ifin+1];
  sum3+=dx*RU[ifin+1];
  //sum4+=dx*phi(so,&RHO0d[ifin+1],Ckbar);
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

    // Calcul du pas de temps
    minXsC=XsC[ideb];
    for(int j=ideb+1; j<=ifin; j++){
      if (minXsC>XsC[j]){
        minXsC=XsC[j];
      }
    }
    dt_new = CFL * minXsC; // pas de temps
    //if ((T-t)<dt){ dt=T-t; }
    dt=pastemps(Tf,t,dt_new,dt);


    for(int i=0; i<=N; i++){
      U[i]=RU[i]/RHO0d[i];
    }

    for(int i=ideb2; i<=N2; i++){
      R2C2dU[i]=RHO0c[i]*RHO0c[i]*(C[i]/RTAU[i])*(C[i]/RTAU[i])*(U[i+1]-U[i]);
    }


    // Calcul de Ustar et Pstar
    for(int i=ideb2+1; i<=N2-1; i++){
      Pstar[i]=P[i] - (dt/(2*dx))*RHO0c[i]*(C[i]/RTAU[i])*(C[i]/RTAU[i])*(U[i+1]-U[i]);
      //Pstar[i]+=xi*((dt/dx)*(dt/dx)/6)*(C[i]/RTAU[i])*(C[i]/RTAU[i])*(2*P[i]-P[i-1]-P[i+1]);
    }

    for(int i=ideb2+1; i<=N2; i++){
      //rhoc=( C[i]+C[i-1] )/( RTAU[i]/RHO0c[i] + RTAU[i-1]/RHO0c[i-1] );
      Ustar[i]=U[i] - (dt/(2*dx))*(1./RHO0d[i])*(P[i]-P[i-1]);
      //Ustar[i]+=xi*((dt/dx)*(dt/dx)/6)*(1./RHO0d[i])*(R2C2dU[i]-R2C2dU[i-1]);
    }

    // Calcul de la solution
    for(int i=ideb; i<=ifin; i++){
      RTAU[i]+=(dt/dx)*(Ustar[i+1]-Ustar[i]);
      RU[i]-=(dt/dx)*(Pstar[i]-Pstar[i-1]);  
      REPS[i]-=(dt/dx)*Pstar[i]*(Ustar[i+1]-Ustar[i]);
      REKIN[i]-=(dt/dx)*Ustar[i]*(Pstar[i]-Pstar[i-1]);
    }
    RU[ifin+1]-=(dt/dx)*(Pstar[ifin+1]-Pstar[ifin]);
    REKIN[ifin+1]-=(dt/dx)*Ustar[ifin+1]*(Pstar[ifin+1]-Pstar[ifin]);


    // cond aux limites  symétrie
    err=condlim(ind_cond, ideb,ifin, nbghosts, RTAU, 0);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin, nbghosts, REPS, 0);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin+1, nbghosts, RU, iu);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin+1, nbghosts, REKIN, 0); // ATTENTION frontières aux mailles



    // Energy kinetic fix
    // CRAS 2016 == JCP 2019 Dakin & al.  (a l'ordre 2 d'espace)
    if(z==1 || z==2){ 

      for (int i=ideb2; i<=N2; i++){
        DK[i]=REKIN[i]-RU[i]*RU[i]/(2.0*RHO0d[i]);
        //printf("i: %d, DK= %.8lf\n",i,DK[i] );
      }  
      for (int i=ideb; i<=ifin; i++){
        REKIN[i]-=DK[i];
        REPS[i]+=(DK[i]+DK[i+1])/2; 
      }
      REKIN[ifin+1]-=DK[ifin+1];

      // On recalcule les variables dépendant de REPS
      // cond aux limites  symétrie  
      err=condlim(ind_cond, ideb,ifin, nbghosts, REPS, 0);
      if(err){return err;}
      err=condlim(ind_cond, ideb,ifin, nbghosts, REKIN, 0);
      if(err){return err;}

    } //// fin kin fix


    // Mise à jour des tableaux
    if(etain){
      for(int i=ideb2; i<=N2-1; i++){ ////
        //EGS(tst,RTAU[i]/RHO0c[i], REPS[i]/RHO0c[i], &P[i], &C[i]);
        /*
        err=fPTC(RTAU[i]/RHO0c[i], REPS[i]/RHO0c[i], &P[i], &C[i], &T[i], &FM[i],
               nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, 
               SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
        */
        if(err){return err;}
      }
    }
    else{
      for(int i=ideb2; i<=N2-1; i++){ ////
        EGS(tst,RTAU[i]/RHO0c[i], REPS[i]/RHO0c[i], &P[i], &C[i]);
      }
    }


    ///////////////
    // Energie du système
    sum1=0; sum3=0.; sum4=0.; sum5=0.;
    for(int i=ideb; i<=ifin; i++){
      sum1+=dx*(REPS[i]+REKIN[i]);
      sum3+=dx*RU[i];
      sum4+=(X[i+1]-X[i])*(RHO0c[i]/RTAU[i]);
      //sum4+=dx*phi(so,&RHO0c[i], Ckbar);
      sum5+=dx*RTAU[i];
    }
    sum1+=dx*REKIN[ifin+1];
    sum3+=dx*RU[ifin+1];
    Etot[n+1]=sum1;
    //IMPUL[n+1]=sum3;
    IMPUL[n+1]=sum3;
    /* 
    if (tst==0){
      IMPUL[n+1]-=0.9*t; 
    }
    */
    Mtot[n+1]=sum4;
    VOL[n+1]=sum5;
    
    /*
    for (int i = ideb; i <= ifin; i++){
      printf("i=%d, Rtau= %lf, RU=%lf, REPS=%lf REKIN= %lf\n",i,RTAU[i],RU[i],REPS[i],REKIN[i] );
    } 
    printf("\n");
    */

    // Mise a jour de XsC
    for(int i=ideb; i<=ifin; i++){
      XsC[i]=(X[i+1]-X[i])/C[i];
    }
     
        
	  t=t+dt;
	  n=n+1;
    if(aff==1){  printf("n= %d, t= %g, dt= %g\n",n,t,dt);  }
  }
//==========================================================================

  printf("Nombre d'itérations : %d\n",n);
  printf("Temps final : %g\n",t);



  // Tableaux pour plot
  for (int i=ideb; i<=ifin; i++){
    TAU[i]=RTAU[i]/RHO0c[i];
    EPS[i]=REPS[i]/RHO0c[i];
    U[i]=RU[i]/RHO0d[i];
    Xc[i]=(X[i+1]+X[i])/2;

    //printf("i=%d, TAUav=%lf, TAU=%lf, RHO0c=%lf\n",i,TAUav[i],RTAU[i],RHO0c[i] );
  }
  U[ifin+1]=RU[ifin+1]/RHO0d[ifin+1];

  for (int i=ideb; i<=ifin; i++){
    E[i]=EPS[i]+((U[i]+U[i+1])/2)*((U[i]+U[i+1])/2)/2;
  }

  // printf_sol
  err=print_sol(0, ideb, ifin, X, Xc, TAU, U, E, P, EPS);
  if(err){return err;}

  err=print_nrj(n, Etot, IMPUL, Mtot, VOL);
  if(err){return err;}
  


  //free(WG); free(WD);
  free(W); free(E);

  free(X); free(Xc);
  free(RHO0d); free(RHO0c);
  
  free(RTAU); free(RU); free(REPS); 

  free(P); free(C); free(XsC);
  
  free(Pstar); free(Ustar); 
  
  // Plot
  free(TAU); free(EPS); free(U);  

  // Energy kinetic fix
  free(R2U2); free(DK);
  free(REKIN); 
  free(Etot);
  free(RU2);
  free(IMPUL);
  free(Mtot);
  free(VOL);

  return 0;
  
}
