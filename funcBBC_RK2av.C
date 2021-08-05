#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "head.h"



// SCHEMA BBC PREDICTOR-CORRECTOR FORMULE EN ENERGIE INTERNE   JCP 2016

/////////////////////  MAIN  ////////////////////
int funcBBC_RK2av(int tst, double T, double a, double b, int nx, int Nmax, double CFL, int aff, int ind_cond, int iu, int q , double Cq, double Cl){

  double* W=malloc(3*sizeof(double));

	// variables
	double t=0;
	double dt,dt_new;
	int n=0,err=0;
  double minDMsRC; 
	char* W0;
  double sum1, sum2, sum3, sum4, sum5;
  int N;
  int nbghosts=10; // nombre de mailles fantomes pour les variables du schéma (\rho, \tau, u, \epsilon, ekin)
  int nbghosts2=4; 
  // l'idee  : appliquer les cond aux lim sur les var. sur nbghosts mailles qui calculer la valeur des fcts sur les nbghosts2 mailles
  int ideb, ifin;
  double temp;

  double dx=(b-a)/nx; // pas d'espace


  // Calcul de la taille des tableaux
  N=nx+2*nbghosts;
  ideb=nbghosts;
  ifin=nbghosts+nx-1;
  int N2=nbghosts+nx+nbghosts2; // cond au lim de Q et P sur i=nbghosts- nbghosts2, N2

  // Allocation des tableaux
  double* X=malloc((N+1)*sizeof(double)); // maillage
  double* Xc=malloc((N)*sizeof(double)); // maillage
  double* E=malloc((N)*sizeof(double)); 
  // PW
	double* RHO0c=malloc((N)*sizeof(double));   //  rho0  (centrées aux mailles)
	double* RHO0d=malloc((N+1)*sizeof(double)); //  rho0  (dé-centrées | frontières)
   
	// solution 
  double* TAU=malloc((N)*sizeof(double)); //  tau  (centrées aux mailles)
	double* RHO=malloc((N)*sizeof(double)); //  tau  (centrées aux mailles)
	double* U=malloc((N+1)*sizeof(double)); //  u (frontières)
  double* EPS=malloc((N)*sizeof(double)); //  epsilon (centrées aux mailles)
  double* EKIN=malloc((N+1)*sizeof(double)); //  energie totale (centrées aux mailles)

  // first stage
  double* RHO2=malloc((N)*sizeof(double)); //  tau au temsp n+1/2 (centrées aux mailles)
  double* U2=malloc((N+1)*sizeof(double)); //  u au temsp n+1/2 (frontières)
  double* Umoy=malloc((N+1)*sizeof(double)); //  u au temsp n+1/4 (frontières)
  double* EPS2=malloc((N)*sizeof(double)); //  epsilon au temsp n+1/2 (centrées aux mailles)
  double* P2=malloc((N)*sizeof(double));
  double* C2=malloc((N)*sizeof(double));
  double* X2=malloc((N+1)*sizeof(double)); // maillage



  double* P=malloc((N)*sizeof(double));
  double* Q=malloc((N)*sizeof(double));   // artificial viscosity
  double* C=malloc((N)*sizeof(double));
  double* DMsRC=malloc((N)*sizeof(double));  // dX/c   centres aux mailles
  double* DMc=malloc((N)*sizeof(double));  // Delta m = (X[i+1]-X[i])*(RHOC[i]/RTAU[i])   centres aux mailles
  double* DMd=malloc((N+1)*sizeof(double));  // Delta m = (Dm + Dm)/2   dé-centres aux mailles
  
  
  

	// Maillage
	for(int i=0; i<=N; i++){
	  X[i]=(a-nbghosts*dx)+dx*i;
    //printf("i= %d, X=%lf\n",i,X[i]);
	}
  

  // Initialisation 
  double x;
  // variables défini aux centres des mailles
	for(int i=0; i<=N-1; i++){
	  x=(X[i]+X[i+1])/2.0;
    err=u0(tst,a,b, x, W);
    if(err){return err;}
    RHO0c[i]=W[0];
    RHO[i]=W[0];
    EPS[i]=epsilonEOS(tst, 1./W[0], W[2]);
    //printf("rho=%lf, p=%lf\n",W[0], W[2] );
	  //printf("i= %d, RHO=%lf, EPS=%lf\n",i,RHO[i],EPS[i]);
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



  //Init des tableaux  
  for(int i=0; i<=N-1; i++){
    EGS(tst,1./RHO[i], EPS[i], &P[i], &C[i]);
    DMc[i]=dx*RHO0c[i];
    DMsRC[i]=DMc[i]/(RHO[i]*C[i]);
    //printf("i=%d P= %lf, C=%lf, U=%lf, EPS=%lf, XsC=%lf \n",i,P[i],C[i],U[i],REPS[i]/RHO0c[i], XsC[i]);           
  }

  for(int i=1; i<=N; i++){
    DMd[i]=(DMc[i-1]+DMc[i])/2;          
  }

  // energie cinetique
  for(int i=ideb; i<=ifin+1; i++){ 
    EKIN[i]=DMd[i]*U[i]*U[i]/2;
  }


	// minimum de DMsRC  (calcul de dt)
  minDMsRC=DMsRC[ideb];
  for(int j=ideb+1; j<=ifin; j++){
    if(minDMsRC>DMsRC[j]){
      minDMsRC=DMsRC[j];
    }
  }

  printf("minDMsRC= %.8lf\n",minDMsRC );
  dt = CFL * minDMsRC; // pas de temps
	printf("nx= %d et dt= %.13lf\n", nx, dt);
  
  // energie totale du système
  int nbiter=5*ceil(T/dt);
  double* Etot=malloc(nbiter*sizeof(double));
  double* IMPUL=malloc(nbiter*sizeof(double)); // impulsion
  double* Mtot=malloc(nbiter*sizeof(double));  // masse totale
  double* VOL=malloc(nbiter*sizeof(double));   // volume total
  

  sum1=0; sum3=0.; sum4=0.; sum5=0.;
  for(int i=ideb; i<=ifin; i++){
    sum1+=DMc[i]*EPS[i] + DMd[i]*EKIN[i];
    sum3+=dx*RHO0c[i]*U[i];
    sum4+=dx*RHO0c[i];
    sum5+=dx*RHO0c[i]/RHO[i];
  }
  sum1+=DMd[ifin+1]*EKIN[ifin+1];
  sum3+=dx*RHO0c[ifin+1]*U[ifin+1];
  Etot[0]=sum1;
  IMPUL[0]=sum3;
  Mtot[0]=sum4;
  VOL[0]=sum5;
  printf("nbiter :%d\n",nbiter);

  /*
  for (int i = ideb; i <=ifin; ++i){
    printf("i=%d, RTAU=%lf, RU=%lf, REPS=%lf\n",i,RTAU[i],RU[i],REPS[i] );
  }
  */

// debut de la boucle ============================================================================= 
	while((t<T)&&(n<Nmax)){
    // cond aux limites     
    err=condlim(ind_cond, ideb,ifin, nbghosts, RHO, 0);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin, nbghosts, EPS, 0);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin+1, nbghosts, U, iu);
    if(err){return err;}

    // Calcul du pas de temps
    minDMsRC=DMsRC[ideb];
    for(int j=ideb+1; j<=ifin; j++){
      if (minDMsRC>DMsRC[j]){
        minDMsRC=DMsRC[j];
      }
    }
    dt_new = CFL * minDMsRC; // pas de temps
    //if ((T-t)<dt){ dt=T-t; }
    dt=pastemps(T,t,dt_new,dt);

    for(int i=nbghosts-nbghosts2; i<=N2-1; i++){
      Q[i]=qvis(q,Cq,Cl,&U[i], 1./RHO[i], C[i], &DMc[i]);           
    }


    // First STAGE
    for(int i=ideb; i<=ifin+1; i++){
      U2[i]=U[i] - (dt/(2*DMd[i]))*( (P[i]+Q[i]) - (P[i-1]+Q[i-1]) );           
    }
    err=condlim(ind_cond, ideb,ifin+1, nbghosts, U2, iu);
    if(err){return err;}

    for(int i=0; i<=N; i++){
      X2[i]=X[i] + (dt/2)*U2[i];           
    }

    for(int i=ideb; i<=ifin; i++){
      RHO2[i]=DMc[i]/(X2[i+1]-X2[i]); 
      EPS2[i]=EPS[i] -(dt/(2*DMc[i]))*(P[i]+Q[i])*(U2[i+1]-U2[i]);           
    }
    err=condlim(ind_cond, ideb,ifin, nbghosts, RHO2, 0);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin, nbghosts, EPS2, 0);
    if(err){return err;} 
    for(int i=0; i<=N-1; i++){
      EGS(tst,1./RHO2[i], EPS2[i], &P2[i], &C2[i]);
    }
    

    // 2nd STAGE

    for(int i=ideb; i<=ifin+1; i++){
      temp=(dt/DMd[i])*( (P2[i]+Q[i]) - (P2[i-1]+Q[i-1]) );
      Umoy[i]=U[i] - temp/2;
      U[i]-=temp;           
    }
    err=condlim(ind_cond, ideb,ifin+1, nbghosts, U, iu);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin+1, nbghosts, Umoy, iu);
    if(err){return err;}

    for(int i=0; i<=N; i++){
      X[i]+=dt*Umoy[i];           
    }

    for(int i=ideb; i<=ifin; i++){
      RHO[i]=DMc[i]/(X[i+1]-X[i]); 
      EPS[i]-=(dt/DMc[i])*(P2[i]+Q[i])*(Umoy[i+1]-Umoy[i]);           
    }
    err=condlim(ind_cond, ideb,ifin, nbghosts, RHO, 0);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin, nbghosts, EPS, 0);
    if(err){return err;}
    
    for(int i=0; i<=N-1; i++){
      EGS(tst,1./RHO[i], EPS[i], &P[i], &C[i]);
    }

    // energie cinetique
    for(int i=ideb; i<=ifin+1; i++){ 
      EKIN[i]=U[i]*U[i]/2;
    }



    ///////////////
    // Energie du système
    sum1=0; sum3=0.; sum4=0.; sum5=0.;
    for(int i=ideb; i<=ifin; i++){
      sum1+=DMc[i]*EPS[i]+DMd[i]*EKIN[i];
      sum3+=(X[i+1]-X[i])*U[i]*RHO[i];
      sum4+=(X[i+1]-X[i])*RHO[i];
      sum5+=(X[i+1]-X[i])*RHO0c[i]/RHO[i];
    }
    sum1+=DMd[ifin+1]*EKIN[ifin+1];
    sum3+=(X[ifin+1]-X[ifin])*U[ifin+1]*RHO[ifin+1];
    Etot[n+1]=sum1;
    IMPUL[n+1]=sum3;
    Mtot[n+1]=sum4;
    VOL[n+1]=sum5;
    

    /*
    for (int i = ideb; i <= ifin; i++){
      printf("i=%d, RHO= %lf, U=%lf, EPS=%lf\n",i,RHO[i],U[i],EPS[i] );
    } 
    printf("\n");
    */

    // Mise a jour de XsC
    for(int i=ideb; i<=ifin; i++){
      DMsRC[i]=DMc[i]/(RHO[i]*C[i]);
    }
     
        
	  t=t+dt;
	  n=n+1;
    if(aff==1){
	    printf("n= %d, t= %.15lf, dt= %.15lf\n",n,t,dt);
	  }
  }
//==========================================================================

  printf("Nombre d'itérations : %d\n",n);
  printf("Temps final : %.15lf\n",t);


  for(int i=ideb;i<=ifin;i++){
    Xc[i]=(X[i+1]+X[i])/2;
    E[i]=EPS[i] + (EKIN[i+1]+EKIN[i])/2;
    TAU[i]=1./RHO[i];
  }

  // printf_sol
  err=print_sol(0, ideb, ifin, X, Xc, TAU, U, E, P, EPS);
  if(err){return err;}

  err=print_nrj(n, Etot, IMPUL, Mtot, VOL);
  if(err){return err;}


  free(W); free(E);

  free(X); free(X2); free(Xc);
  free(RHO0d); free(RHO0c);
  free(DMc); free(DMd);

  free(TAU);
  
  free(RHO); free(U); free(EPS);

  free(P); free(Q); free(C); free(DMsRC);
  free(EKIN);
 
  free(P2); free(C2);

  free(RHO2); free(U2); free(EPS2); 
  free(Umoy);
  
  //  
  free(Etot); 
  free(IMPUL);
  free(Mtot);
  free(VOL);

  return 0;
  
}
