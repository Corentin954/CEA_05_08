#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "head.h"



// SCHEMA DE TYPE RUNGE-KUTTA FORMULE EN ENERGIE INTERNE

/////////////////////  MAIN  ////////////////////
int funcRKint(int sch, int tst, double T, double a, double b, int nx, int Nmax, double CFL, int aff, int ind_cond, int iu, int q , double Cq, double Cl, int z, int dpi, int so){

  printf("CFL= %lf\n",CFL );

  double* W=malloc(3*sizeof(double));


	// variables
	double t=0;
	double dt, dt_new;
	int n=0, err=0;
  double minXsC; 
	char* W0;
  double sum1, sum2, sum3, sum4, sum5;
  int N;
  int nbghosts=10; // nombre de mailles fantomes pour les variables du schéma (\rho, \tau, u, \epsilon, ekin)
  int nbghosts2=6; // nombre de mailles fantomes pour les fonctions qui dependent des variables ci-dessus
  // l'idee  : appliquer les ocnd aux lim sur les var. sur nbghosts mailles qui calculer la valeur des fcts sur les nbghosts2 mailles
  int ideb, ifin;
  double** Ck; double** Ckbar; double** dk; double** Qbar; double** rk;
  int sodPi; // spatial order for dPi


  
  // spatial order for dPi
  sodPi=so;


  // Coefficient
  Ck=alloctabd(5,5); Ckbar=alloctabd(5,5);
  dk=alloctabd(5,5); Qbar=alloctabd(5,5);
  rk=alloctabd(5,5);
  initCoef(Ck, Ckbar, dk, Qbar, rk);
  
  // Coefficient de Btucher
  int to;  // time order
  double** A; double* THETA; double* ALPHA;
  // RK1 et RK2
  if(sch==0||sch==1){   to=1;  }
  // SSPRK3(3,3) Spiteri-Ruuth
  else if (sch==2){    to=2; }
  // RK4 Classic
  else if (sch==3){   to=3;  }
  // RK5 Cash-Karp 
  else if (sch==4){   to=5;  }
  // Dormand-Prince 
  else if (sch==5){   to=6;  }
  // HEUN  ordre 3
  else if (sch==6){   to=2;  }
  // RK3 RALSTON
  else if (sch==7){   to=2;  }
  // RK2 ou RK3 Bogaki-Shampine
  else if (sch==8 || sch==9){    to=3;  }
  // KUTTA odre 3      
  else if (sch==10){   to=2;  }
  // KUTTA ordre 2 (Explicit Midpoint Method) 
  else if (sch==11){   to=1;  }
  // KUTTA ordre 2 Ralston
  else if (sch==12){   to=1;  }
  // EULER forward
  else if (sch==13){   to=1;  }
  // SSP RK4 OPTIMUM
  else if (sch==14){   to=4;  }
  // SSPRK3(4,3) Spiteri-Ruuth
  else if (sch==15){   to=3;  }
  // SSPRK3(5,3) Spiteri-Ruuth
  else if (sch==16){   to=4;  }
  // SPPRK1(3,1) Spiteri-Ruuth
  else if (sch==17){   to=2;  }
  // SPPRK2(3,2) Spiteri-Ruuth
  else if (sch==18){   to=2;  }
  // SPPRK2(4,2) Spiteri-Ruuth
  else if (sch==19){   to=3;  }
  // SPPRK1(2,1) Spiteri-Ruuth
  else if (sch==20){   to=1;  }
  // SPPRK3(4,3) Spiteri-Ruuth
  else if (sch==21){   to=3;  }
  else{ printf("Erreur dans le choix de sch (RKint)\n");}
  

  A=alloctabd(to, to);
  THETA=malloc((to+1)*sizeof(double));
  ALPHA=malloc(to*sizeof(double));
	initButcher(sch, A, THETA,ALPHA);
  /*
  for(int i=0; i<to; i++){
    for(int j=0; j<=i; j++){
      printf("A[%d][%d]= %lf ",i,j,A[i][j] );
    }
    printf("\n");
  }
  for (int i = 0; i < to+1; ++i){
    printf("theta[%d]= %lf ",i,THETA[i] );
  }
  printf("\n");
  */

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
   
	// solution AVERAGE
	double* RTAU=malloc((N)*sizeof(double)); //  rho0.tau  (centrées aux mailles)
	double* RU=malloc((N+1)*sizeof(double)); //  rho0.u (frontières)
  double* REPS=malloc((N)*sizeof(double)); //  rho0.epsilon (centrées aux mailles)

  // solution POINT-WISE
  double* RTAUpw=malloc((N)*sizeof(double)); //  rho0.tau  (centrées aux mailles)
  double* RUpw=malloc((N+1)*sizeof(double)); //  rho0.u (frontières)
  double* REPSpw=malloc((N)*sizeof(double)); //  rho0.epsilon (centrées aux mailles)

  // AVERAGE (plot)
  double* TAUav=malloc((N)*sizeof(double)); //  tau  (centrées aux mailles)
  double* EPSav=malloc((N)*sizeof(double)); //  epsilon (centrées aux mailles)  
  double* Uav=malloc((N+1)*sizeof(double)); //  u  (frontières aux mailles)
  double* Pav=malloc((N)*sizeof(double));   //  p (centrées aux mailles)  
  double* E=malloc((N)*sizeof(double));   //  p (centrées aux mailles)  

  // Point-Wise
  //double* TAU=malloc((N)*sizeof(double)); //  tau  (centrées aux mailles)
  //double* EPS=malloc((N)*sizeof(double)); //  epsilon (centrées aux mailles)
  double* P=malloc((N)*sizeof(double));
  double* Q=malloc((N)*sizeof(double));   // artificial viscosity
  double* U=malloc((N+1)*sizeof(double)); // vitesse donc frontières aux mailles
  double* C=malloc((N)*sizeof(double));
  double* XsC=malloc((N)*sizeof(double));  // dX/c   centres aux mailles
  double* DM=malloc((N)*sizeof(double));  // Delta m = (X[i+1]-X[i])*(RHOC[i]/RTAU[i])   centres aux mailles

  double* dPI=malloc((N)*sizeof(double)); // delta(P+Q)
  double** dPIstar=alloctabd(2,N); // delta(P+Q)
  
  // tableaux temporaires
  // AVERAGE
  double** RTAUstar=alloctabd(to,N); //  rho0.tau  (centrées aux mailles)
  double** RUstar=alloctabd(to,N+1); //  rho0.u (frontières)
  double** REPSstar=alloctabd(to,N); //  rho0.tau (centrées aux mailles)
  // Point-Wise
  double** Pstar=alloctabd(to,N);
  double** Qstar=alloctabd(to,N);  // artificial viscosity
  double** Cstar=alloctabd(to,N);
  double** Ustar=alloctabd(to,N+1); // vitesse donc frontières aux mailles
  double** TAUstar=alloctabd(to,N); //  rho0.tau  (centrées aux mailles)
  double** EPSstar=alloctabd(to,N); //  rho0.tau  (centrées aux mailles)
  double** DMstar=alloctabd(to,N);  // Delta m = (X[i+1]-X[i])*(RHOC[i]/RTAU[i])   centres aux mailles
  double** Xstar=alloctabd(to,N+1); // maillage

  // kinetic fix  
  double* REKINpw=malloc((N+1)*sizeof(double));
  double** REKINstar=alloctabd(to,N+1);
  double* REKIN=malloc((N+1)*sizeof(double));
  double* DK=malloc((N+1)*sizeof(double));
  double* R2U2pw=malloc((N+1)*sizeof(double)); //(point-wise)

  double* RU2pw=malloc((N+1)*sizeof(double)); //(point-wise)


	// Maillage
	for(int i=0; i<=N; i++){
	  X[i]=(a-nbghosts*dx)+dx*i;
    //printf("i= %d, X=%lf\n",i,X[i]);
	}
  

  // Initialisation POINT-WISE
  double x;
  // variables défini aux centres des mailles
	for(int i=0; i<=N-1; i++){
    RTAUpw[i]=1.0;
    RTAU[i]=1.0;
	  x=(X[i]+X[i+1])/2.0;
    err=u0(tst,a,b, x, W);
    if(err){return err;}
    RHO0c[i]=W[0];
    REPSpw[i]=W[0]*epsilonEOS(tst, 1./W[0], W[2]);
    //printf("rho=%lf, p=%lf\n",W[0], W[2] );
	  //printf("i= %d, rho0.tau=%f, rho0.eps=%f, rho0c=%f\n",i,RTAUpw[i],REPSpw[i],RHO0c[i] );
	}
  // variables défini aux frontières des mailles
  for(int i=0; i<=N; i++){
    x=X[i];
    err=u0(tst,a,b, x, W);
    if(err){return err;}
    RHO0d[i]=W[0];
    RUpw[i]=W[0]*W[1];
    REKINpw[i]=W[0]*W[1]*W[1]/2.0;
    //printf("u=%.8lf\n",W[1] );
    //printf("i= %d,  rho0.u=%f,  rho0d=%f\n",i,RUpw[i],RHO0d[i] );
  }



  //Init des tableaux   PW
  for(int i=0; i<=N-1; i++){
    U[i]=RUpw[i]/RHO0d[i];
    EGS(tst,RTAUpw[i]/RHO0c[i], REPSpw[i]/RHO0c[i], &P[i], &C[i]);
    XsC[i]=(X[i+1]-X[i])/C[i];
    DM[i]=(X[i+1]-X[i])*(RHO0c[i]/RTAUpw[i]);
    //printf("i=%d P= %lf, C=%lf, U=%lf, EPS=%lf, XsC=%lf \n",i,P[i],C[i],U[i],REPS[i]/RHO0c[i], XsC[i]);           
  }
  U[N]=RUpw[N]/RHO0d[N];

  for(int i=nbghosts-nbghosts2; i<=N2-1; i++){
    Q[i]=qvis(q,Cq,Cl,&U[i], RTAUpw[i]/RHO0c[i], C[i], &DM[i]);
    //printf("i=%d Q= %lf \n",i,Q[i]);           
  }


  // Initialisation de REPS et RU (average)
  for(int i=ideb; i<=ifin; i++){
    REPS[i]=phi(so, &REPSpw[i], Ckbar);
    RU[i]=phi(so, &RUpw[i], Ckbar);
    REKIN[i]=phi(so, &REKINpw[i], Ckbar);
    //printf("i=%d  REPS=%lf, RU=%lf\n",i,REPS[i],RU[i]);           
  }
  RU[ifin+1]=phi(so, &RUpw[ifin+1], Ckbar);
  REKIN[ifin+1]=phi(so, &REKINpw[ifin+1], Ckbar);



	// minimum de XsC  (calcul de dt)
  minXsC=XsC[ideb];
  for(int j=ideb+1; j<=ifin; j++){
    if(minXsC>XsC[j]){
      minXsC=XsC[j];
    }
  }

  printf("minXsC= %.8lf\n",minXsC );
  dt = CFL * minXsC; // pas de temps
  printf("nx= %d et dt= %.13lf\n", nx, dt);
  
  // energie totale du système
  int nbiter=5*ceil(T/dt);
  printf("nbiter :%d\n",nbiter);
  double* Etot=malloc(nbiter*sizeof(double));
  double* IMPUL=malloc(nbiter*sizeof(double)); // impulsion
  double* Mtot=malloc(nbiter*sizeof(double)); // masse totale
  double* VOL=malloc(nbiter*sizeof(double)); // volume total
  
  // Energie du système
  for(int i=0; i<=N; i++){
    R2U2pw[i]=RUpw[i]*RUpw[i];
  }

  sum1=0; sum2=0; sum3=0.; sum4=0.; sum5=0.;
  for(int i=ideb; i<=ifin; i++){
    sum1+=dx*REPS[i];
    sum2+=(dx/2.0)*phidiv(so, &R2U2pw[i], &RHO0d[i], Ckbar);
    sum3+=dx*RU[i];
    sum4+=dx*phi(so,&RHO0c[i],Ckbar);
    sum5+=dx*RTAU[i];
  }
  sum2+=(dx/2.0)*phidiv(so, &R2U2pw[ifin+1], &RHO0d[ifin+1], Ckbar);
  sum3+=dx*RU[ifin+1];
  //sum4+=dx*phi(so,&RHO0d[ifin+1],Ckbar);
  Etot[0]=sum1+sum2;
  IMPUL[0]=sum3;
  Mtot[0]=sum4;
  VOL[0]=sum5;
  

  /*
  for (int i = ideb; i <=ifin; ++i){
    printf("i=%d, RTAU=%lf, RU=%lf, REPS=%lf\n",i,RTAU[i],RU[i],REPS[i] );
  }
  */

// debut de la boucle ============================================================================= 
	while((t<T)&&(n<Nmax)){
    // cond aux limites  symétrie    
    err=condlim(ind_cond, ideb,ifin, nbghosts, RTAU, 0);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin, nbghosts, REPS, 0);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin+1, nbghosts, RU, iu);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin+1, nbghosts, REKIN, 0); // attention frontières aux mailles
    if(err){return err;}

    // Calcul du pas de temps
    minXsC=XsC[ideb];
    for(int j=ideb+1; j<=ifin; j++){
      if (minXsC>XsC[j]){
        minXsC=XsC[j];
      }
    }
    dt_new = CFL * minXsC; // pas de temps
    //if ((T-t)<dt){ dt=T-t; }
    dt=pastemps(T,t,dt_new,dt);

    // dPI=delta(P+Q)
    if(dpi==1){
      for(int i=nbghosts-nbghosts2; i<=N2-1; i++){ ////
        dPI[i]=delta(sodPi,0,&P[i], dk)+delta(sodPi,0,&Q[i], dk);
      }
    }
    
	  // Boucle RUNGE-KUTTA
    for(int k=0; k<to; k++){  // ordre en temps
  	  // 
      for(int i=ideb; i<=ifin; i++){
        RTAUstar[k][i]=RTAU[i]+A[k][0]*(dt/dx)*(U[i+1]-U[i]);
        if(dpi==0){  RUstar[k][i]=RU[i]-A[k][0]*(dt/dx)*((P[i]+Q[i])-(P[i-1]+Q[i-1]));  }
        else if(dpi==1){  RUstar[k][i]=RU[i]-A[k][0]*(dt/dx)*phi(sodPi, &dPI[i], Ckbar);  }
        REPSstar[k][i]=REPS[i]-A[k][0]*(dt/dx)*YdZ(so,0,&P[i],&Q[i],&U[i],dk,Ckbar);
        REKINstar[k][i]=REKIN[i]-A[k][0]*(dt/dx)*YdZ(so,1,&P[i],&Q[i],&U[i],dk,Ckbar);
      }
      if(dpi==0){  RUstar[k][ifin+1]=RU[ifin+1]-A[k][0]*(dt/dx)*((P[ifin+1]+Q[ifin+1])-(P[ifin]+Q[ifin]));  }
      else if(dpi==1){  RUstar[k][ifin+1]=RU[ifin+1]-A[k][0]*(dt/dx)*phi(sodPi, &dPI[ifin+1], Ckbar);  }
      REKINstar[k][ifin+1]=REKIN[ifin+1]-A[k][0]*(dt/dx)*YdZ(so,1,&P[ifin+1],&Q[ifin+1],&U[ifin+1],dk,Ckbar); 

      for(int i=nbghosts - nbghosts2; i<=N2; i++){
        Xstar[k][i]= X[i] + A[k][0]*dt*U[i];
      }
      

      for(int j=1; j<=k; j++){ // calcul des sous pas de temps
        //
        for(int i=ideb; i<=ifin; i++){
    	    RTAUstar[k][i]+=A[k][j]*(dt/dx)*(Ustar[j-1][i+1]-Ustar[j-1][i]);
          if(dpi==0){   RUstar[k][i]-=A[k][j]*(dt/dx)*((Pstar[j-1][i]+Qstar[j-1][i])-(Pstar[j-1][i-1]+Qstar[j-1][i-1]));   }
          else if(dpi==1){    RUstar[k][i]-=A[k][j]*(dt/dx)*phi(sodPi, &dPIstar[j-1][i], Ckbar);   }
          REPSstar[k][i]-=A[k][j]*(dt/dx)*YdZ(so,0,&Pstar[j-1][i],&Qstar[j-1][i],&Ustar[j-1][i],dk,Ckbar);
          REKINstar[k][i]-=A[k][j]*(dt/dx)*YdZ(so,1,&Pstar[j-1][i],&Qstar[j-1][i],&Ustar[j-1][i],dk,Ckbar); 
    	  }
        if(dpi==0){   RUstar[k][ifin+1]-=A[k][j]*(dt/dx)*((Pstar[j-1][ifin+1]+Qstar[j-1][ifin+1])-(Pstar[j-1][ifin]+Qstar[j-1][ifin]));  }
        else if(dpi==1){   RUstar[k][ifin+1]-=A[k][j]*(dt/dx)*phi(sodPi, &dPIstar[j-1][ifin+1], Ckbar);  }
        REKINstar[k][ifin+1]-=A[k][j]*(dt/dx)*YdZ(so,1,&Pstar[j-1][ifin+1],&Qstar[j-1][ifin+1],&Ustar[j-1][ifin+1],dk,Ckbar);

        for(int i=nbghosts - nbghosts2; i<=N2; i++){
          Xstar[k][i]+=A[k][j]*dt*Ustar[j-1][i];
        }
      }

      // cond aux limites
      err=condlim(ind_cond, ideb,ifin, nbghosts, RTAUstar[k], 0);
      if(err){return err;}
      err=condlim(ind_cond, ideb,ifin, nbghosts, REPSstar[k], 0);
      if(err){return err;}
      err=condlim(ind_cond, ideb,ifin+1, nbghosts, RUstar[k], iu);
      if(err){return err;}
      err=condlim(ind_cond, ideb,ifin+1, nbghosts, REKINstar[k], 0);
      if(err){return err;}

      // Energy kinetic fix
      // CRAS 2016
      if(z==11 || z==12){ 
        for(int i=nbghosts-nbghosts2; i<=N2; i++){ 
          RUpw[i]=phi(so, &RUstar[k][i],Ck);
        }
        if(z==11){
          for(int i=nbghosts-nbghosts2; i<=N2; i++){ ////  
            REKINpw[i]=phi(so, &REKINstar[k][i],Ck);
          }
          for (int i=nbghosts-nbghosts2; i<=N2; i++){
            DK[i]=REKINpw[i]-RUpw[i]*RUpw[i]/(2.0*RHO0d[i]);
            //printf("i: %d, DK= %.8lf\n",i,DK[i] );
          }
          for (int i=ideb; i<=ifin; i++){
            REKINstar[k][i]-=phi(so, &DK[i],Ckbar);
            REPSstar[k][i]+=phiQK(so, &DK[i],Qbar); 
          }
          REKINstar[k][ifin+1]-=phi(so,&DK[ifin+1],Ckbar);
        }
        // JCP 2019 Dakin & al.
        else if(z==12){
          for(int i=nbghosts-nbghosts2; i<=N2; i++){ ////
            RU2pw[i]=RUpw[i]*RUpw[i]/(2.0*RHO0d[i]);
            //printf("i: %d, RU2pw= %.8lf\n",i,RU2pw[i] );
          }
          for (int i=ideb; i<=ifin+1; i++){
            DK[i]=REKINstar[k][i]-phi(so, &RU2pw[i], Ckbar);
            //printf("i: %d, DK= %.8lf\n",i,DK[i] );
          } 
          for (int i=ideb; i<=ifin; i++){
            REKINstar[k][i]-= DK[i];
            REPSstar[k][i]+= (DK[i]+DK[i+1])/2.0; 
          }
          REKINstar[k][ifin+1]-= DK[ifin+1]; 
        }

        // On recalcule les variables dépendant de REPS
        // cond aux limites  symétrie  
        err=condlim(ind_cond, ideb,ifin, nbghosts, REPSstar[k], 0);
        if(err){return err;}
        err=condlim(ind_cond, ideb,ifin, nbghosts, REKINstar[k], 0);
        if(err){return err;}
      }
      

      // On veut que U soit PW or PUstar est AV
      for(int i=nbghosts-nbghosts2; i<=N2-1; i++){ ////
        Ustar[k][i]=phi(so, &RUstar[k][i], Ck)/RHO0d[i]; 
        TAUstar[k][i]=phi(so, &RTAUstar[k][i], Ck)/RHO0c[i]; 
        EPSstar[k][i]=phi(so, &REPSstar[k][i],  Ck)/RHO0c[i]; 
      }
      Ustar[k][N2]=phi(so, &RUstar[k][N2], Ck)/RHO0d[N2]; 
      
      // DMstar pour le calcul Qstar
      for(int i=nbghosts-nbghosts2; i<=N2-1; i++){  ///
        DMstar[k][i]=(Xstar[k][i+1]-Xstar[k][i])/TAUstar[k][i];
      }

  	  for(int i=nbghosts-nbghosts2; i<=N2-1; i++){ ////
  	    EGS(tst,TAUstar[k][i], EPSstar[k][i], &Pstar[k][i], &Cstar[k][i]); // On a besoin d'argument PW
        Qstar[k][i]=qvis(q,Cq,Cl,&Ustar[k][i], TAUstar[k][i], Cstar[k][i], &DMstar[k][i]);
  	  }


      // dPIstar
      if(dpi==1){
        for(int i=nbghosts-nbghosts2; i<=N2; i++){ 
          dPIstar[k][i]=delta(sodPi,0,&Pstar[k][i], dk)+delta(sodPi,0,&Qstar[k][i], dk);
        }
      }

    } // fin de boucle sur k ///////////


     
    // Solution au temps n+1
    for(int i=ideb; i<=ifin; i++){
      RTAU[i]+=THETA[0]*(dt/dx)*(U[i+1]-U[i]);
      if(dpi==0){  RU[i]-=THETA[0]*(dt/dx)*((P[i]+Q[i])-(P[i-1]+Q[i-1]));  }
      else if(dpi==1){  RU[i]-=THETA[0]*(dt/dx)*phi(sodPi, &dPI[i], Ckbar);  }
      REPS[i]-=THETA[0]*(dt/dx)*YdZ(so,0,&P[i],&Q[i],&U[i],dk,Ckbar);
      REKIN[i]-=THETA[0]*(dt/dx)*YdZ(so,1,&P[i],&Q[i],&U[i],dk,Ckbar);    
    }
    if(dpi==0){  RU[ifin+1]-=THETA[0]*(dt/dx)*((P[ifin+1]+Q[ifin+1])-(P[ifin]+Q[ifin]));  }
    else if(dpi==1){  RU[ifin+1]-=THETA[0]*(dt/dx)*phi(sodPi, &dPI[ifin+1], Ckbar);  }
    REKIN[ifin+1]-=THETA[0]*(dt/dx)*YdZ(so,1,&P[ifin+1],&Q[ifin+1],&U[ifin+1],dk,Ckbar);

    for(int i=nbghosts-nbghosts2; i<=N2; i++){
      X[i]+=dt*THETA[0]*U[i];
    }


    for (int j=1; j<=to; j++){ // calcul des sous pas de temps
      //
      for(int i=ideb; i<=ifin; i++){
        RTAU[i]+=THETA[j]*(dt/dx)*(Ustar[j-1][i+1]-Ustar[j-1][i]);
        if(dpi==0){   RU[i]-=THETA[j]*(dt/dx)*((Pstar[j-1][i]+Qstar[j-1][i])-(Pstar[j-1][i-1]+Qstar[j-1][i-1]));   }
        else if(dpi==1){    RU[i]-=THETA[j]*(dt/dx)*phi(sodPi, &dPIstar[j-1][i], Ckbar);   }
        REPS[i]-=THETA[j]*(dt/dx)*YdZ(so,0,&Pstar[j-1][i],&Qstar[j-1][i],&Ustar[j-1][i],dk,Ckbar); 
        REKIN[i]-=THETA[j]*(dt/dx)*YdZ(so,1,&Pstar[j-1][i],&Qstar[j-1][i],&Ustar[j-1][i],dk,Ckbar);  
      }
      if(dpi==0){   RU[ifin+1]-=THETA[j]*(dt/dx)*((Pstar[j-1][ifin+1]+Qstar[j-1][ifin+1])-(Pstar[j-1][ifin]+Qstar[j-1][ifin]));  }
      else if(dpi==1){   RU[ifin+1]-=THETA[j]*(dt/dx)*phi(sodPi, &dPIstar[j-1][ifin+1], Ckbar);  }
      REKIN[ifin+1]-=THETA[j]*(dt/dx)*YdZ(so,1,&Pstar[j-1][ifin+1],&Qstar[j-1][ifin+1],&Ustar[j-1][ifin+1],dk,Ckbar);  
      
      for(int i=nbghosts-nbghosts2; i<=N2; i++){ // permet la gestion des cond aux lim du maillage
        X[i]+=dt*THETA[j]*Ustar[j-1][i];
      }

    }

  

    // cond aux limites  symétrie

    err=condlim(ind_cond, ideb,ifin, nbghosts, RTAU, 0);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin, nbghosts, REPS, 0);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin+1, nbghosts, RU, iu);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin+1, nbghosts, REKIN, 0); // ATENTION frontières aux mailles
    if(err){return err;}

    // Mise a jour des rho.Z  Point-Wise
    for(int i=nbghosts-nbghosts2; i<=N2-1; i++){ //// 
      RTAUpw[i]=phi(so, &RTAU[i], Ck);
      REPSpw[i]=phi(so, &REPS[i], Ck);
      RUpw[i]=phi(so, &RU[i], Ck);
    }
    RUpw[N2]=phi(so, &RU[N2], Ck);


    // Energy kinetic fix
    // CRAS 2016
    if(z==1 || z==2 || z==11 || z==12){ 
      if(z==1 || z==11){

        for(int i=nbghosts-nbghosts2; i<=N2; i++){ ////  
          REKINpw[i]=phi(so, &REKIN[i],Ck);
        }
        for (int i=nbghosts-nbghosts2; i<=N2; i++){
          DK[i]=REKINpw[i]-RUpw[i]*RUpw[i]/(2.0*RHO0d[i]);
          //printf("i: %d, DK= %.8lf\n",i,DK[i] );
        }  
        for (int i=ideb; i<=ifin; i++){
          REKIN[i]-=phi(so, &DK[i],Ckbar);
          REPS[i]+=phiQK(so, &DK[i],Qbar); 
        }
        REKIN[ifin+1]-=phi(so,&DK[ifin+1],Ckbar);
      }
      // JCP 2019 Dakin & al.
      else if(z==2 || z==12){
        for(int i=nbghosts-nbghosts2; i<=N2; i++){ ////
          RU2pw[i]=RUpw[i]*RUpw[i]/(2.0*RHO0d[i]);
          //printf("i: %d, RU2pw= %.8lf\n",i,RU2pw[i] );
        }
        for (int i=ideb; i<=ifin+1; i++){
          DK[i]=REKIN[i]-phi(so, &RU2pw[i], Ckbar);
          //printf("i: %d, DK= %.8lf\n",i,DK[i] );
        } 
        for (int i=ideb; i<=ifin; i++){
          REKIN[i]-= DK[i];
          REPS[i]+= (DK[i]+DK[i+1])/2.0; 
        }
        REKIN[ifin+1]-= DK[ifin+1]; 
      }

      // On recalcule les variables dépendant de REPS
      // cond aux limites  symétrie  
      err=condlim(ind_cond, ideb,ifin, nbghosts, REPS, 0);
      if(err){return err;}
      
      // Mise a jour des rho.Z  Point-Wise
      for(int i=nbghosts-nbghosts2; i<=N2-1; i++){  
        REPSpw[i]=phi(so, &REPS[i], Ck);
      }

    } //// fin kin fix


    for(int i=nbghosts-nbghosts2; i<=N2-1; i++){ ////
      DM[i]=(X[i+1]-X[i])*(RHO0c[i]/RTAUpw[i]);
    }

    // Mise à jour des tableaux
    for(int i=nbghosts-nbghosts2; i<=N2-1; i++){ ////  
      EGS(tst,RTAUpw[i]/RHO0c[i], REPSpw[i]/RHO0c[i], &P[i], &C[i]);
      U[i]=RUpw[i]/RHO0d[i];       
    }
    U[N2]=RUpw[N2]/RHO0d[N2];

    for(int i=nbghosts-nbghosts2; i<=N2-1; i++){ //// 
      Q[i]=qvis(q,Cq,Cl,&U[i], RTAUpw[i]/RHO0c[i], C[i],&DM[i]);        
    }

    ///////////////
    // Energie du système
    for(int i=nbghosts-nbghosts2; i<=N2; i++){
      R2U2pw[i]=RUpw[i]*RUpw[i];  // (rho.u)^2
    }

    sum1=0; sum2=0; sum3=0.; sum4=0.; sum5=0.;
    for(int i=ideb; i<=ifin; i++){
      sum1+=dx*REPS[i];
      sum2+=(dx/2.0)*phidiv(so, &R2U2pw[i], &RHO0d[i], Ckbar);
      sum3+=dx*RU[i];
      sum4+=(X[i+1]-X[i])/phidiv(so,&RTAUpw[i],&RHO0c[i], Ckbar);
      //sum4+=dx*phi(so,&RHO0c[i], Ckbar);
      sum5+=dx*RTAU[i];
    }
    sum2+=(dx/2.0)*phidiv(so, &R2U2pw[ifin+1], &RHO0d[ifin+1], Ckbar);
    sum3+=dx*RU[ifin+1];
    Etot[n+1]=sum1+sum2;
    //IMPUL[n+1]=sum3;
    IMPUL[n+1]=sum3; 
    if (tst==0){
      IMPUL[n+1]-=0.9*t; 
    }
    Mtot[n+1]=sum4;
    VOL[n+1]=sum5;
    

    /*
    for (int i = ideb; i <= ifin; i++){
      printf("i=%d, Rtau= %lf, RU=%lf, REPS=%lf\n",i,RTAU[i],RU[i],REPS[i] );
    } 
    printf("\n");
    */

    // Mise a jour de XsC
    for(int i=ideb; i<=ifin; i++){
      XsC[i]=(X[i+1]-X[i])/C[i];
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



  // Tableaux pour plot
  for (int i=ideb; i<=ifin; i++){
    TAUav[i]=phidiv(so, &RTAUpw[i], &RHO0c[i], Ckbar);
    EPSav[i]=phidiv(so, &REPSpw[i], &RHO0c[i], Ckbar);
    Uav[i]=phidiv(so, &RUpw[i], &RHO0d[i], Ckbar);
    Pav[i]=phi(so, &P[i], Ckbar);

    //printf("i=%d, TAUav=%lf, TAUpw=%lf, RHO0c=%lf\n",i,TAUav[i],RTAUpw[i],RHO0c[i] );
  }
  Uav[ifin+1]=phidiv(so, &RUpw[ifin+1], &RHO0d[ifin+1], Ckbar);

  // energie totale
  for(int i=ideb; i<=ifin; i++){ 
    Xc[i]=fxc(so, &X[i], rk);
    E[i]=EPSav[i]+((Uav[i]+Uav[i+1])/2)*((Uav[i]+Uav[i+1])/2)/2;
  }

 
  // printf_sol
  err=print_sol(0, ideb, ifin, X, Xc, TAUav, Uav, E, Pav, EPSav);
  if(err){return err;}
 
  err=print_nrj(n, Etot, IMPUL, Mtot, VOL);
  if(err){return err;}



  free(dPI); freetab(dPIstar);

  free(W);

  free(X);  free(Xc); free(E);
  free(RHO0d); free(RHO0c);
  
  free(RTAU); free(RU); free(REPS);
  free(RTAUpw); free(RUpw); free(REPSpw); 

  freetab(RTAUstar); freetab(RUstar); freetab(REPSstar);
  freetab(TAUstar); freetab(Ustar); freetab(EPSstar); 
  freetab(Pstar); freetab(Qstar); freetab(Cstar);
  freetab(DMstar); freetab(Xstar);

  free(P); free(Q); free(C); free(XsC);
  
  free(Ck); free(Ckbar); free(dk); free(Qbar); free(rk);
  
  // Plot
  free(TAUav); free(EPSav); 
  free(Uav); free(Pav);  

  // Energy kinetic fix
  free(R2U2pw); free(DK);
  free(REKIN); free(REKINpw);
  freetab(REKINstar); 
  free(Etot); 
  free(RU2pw);
  free(IMPUL);
  free(Mtot);
  free(VOL);

  return 0;
  
}
