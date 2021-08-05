#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "head_etain.h"
#include "head.h"


/*
  Schéma lagrangien de Godunov 
  ( c.f JCP 2009)
*/

// Met a jour les condition aux limites pour la variables Z
// Différent de condlim car U est centré aux mailles
/* Si Z est centré aux maille  jfin= ifin
      Z est frontières aux mailles jfin=ifin+1

  int iu = 0 si cond sym
      iu = 1 si cond anti-sym
*/
int condlim_Gtot(int ind_cond, int jdeb,int jfin, int nbghosts, double* Z, int ind){
  if(ind_cond==0){
    if (ind==0){
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=Z[jdeb+j];
        // à droite
        Z[jfin+1+j]=Z[jfin-j];
      }
    }
    else if (ind==1){
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=-Z[jdeb+j];
        // à droite
        Z[jfin+1+j]=-Z[jfin-j];
      }
    }
    return 0;
  }
  // Onde Acoustique 
  else if(ind_cond==1){
    if(ind==0){
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=Z[jfin-j];
        // à droite
        Z[jfin+1+j]=Z[jdeb+j];
      }
    }
    else if(ind==1){ // pour l'energie (variable au carré)
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=Z[jfin+j];
        // à droite
        Z[jfin+1+j]=Z[jdeb-j];
      }
    }
    return 0;
  }
  if(ind_cond==2){
    if (ind==0){
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=Z[jdeb+j];
        // à droite
        Z[jfin+1+j]=Z[jfin-j];
      }
    }
    else if (ind==1){
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=-Z[jdeb+j];
        // à droite
        Z[jfin+1+j]=Z[jfin-j];
      }
    }
    return 0;
  }
  else{
    printf("Erreur choix ind_cond (condlim)\n");
    return 5;
  }
}




/////////////////////  MAIN  ////////////////////
int funcGtot_etain(int sch, int tst, double Tf, double a, double b, int nx, int Nmax, double CFL, int aff, int ind_cond, int iu, int Ntau, double dt_first, double Ttau, int cin_phase, int sch_cin_phase,
	                 int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, int** tabZONE, int* tab_Nb_ZONE, 
                   double* tabVx, double* tabEx, double* tabVy, double* tabEy,
                   double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC){
	

  int etain=tst/100;
	
	// variables
	double t, dt, dt_new;
	int n=0;
	int err=0;
	double rhoc, rhoc1;
	double minXsC;
	char* W0;
	double mu;
	double phiplusU, phimoinsU;
	double phiplusP, phimoinsP;
	double rhoc2, rhoc21;
	double ms, mi;
	double sum1, sum2, sum3, sum4, sum5;
	double dX;
  int zone;

	double dx=(b-a)/nx; // pas d'espace

	// Calcul de la taille des tableaux
	int nbghosts=2; // nombre de mailles fantomes pour les variables du schéma (\rho, \tau, u, \epsilon, ekin)
  int N=nx+2*nbghosts;
  int ideb=nbghosts;
  int ifin=nbghosts+nx-1;


  // Allocation des tableaux
  double* X=malloc((N+1)*sizeof(double));   // maillage
  double* Xc=malloc((N)*sizeof(double));    // maillage
	double* M=malloc((N)*sizeof(double));     // masse donc centrées aux mailles
	double* RHO0c=malloc((N)*sizeof(double)); // masse donc centrées aux mailles
	//double** U=alloctabd(nx+2,3);           // solution  (tau, u, e)
	
	double* TAU=malloc((N)*sizeof(double)); //  rho0.tau  (centrées aux mailles)
	double* U=malloc((N)*sizeof(double));   //  rho0.u (frontières)
  double* E=malloc((N)*sizeof(double));   //  rho0.epsilon (centrées aux mailles)

	double* XsC=malloc((N)*sizeof(double));  // dX/c   centres aux mailles
  double* P=malloc((N)*sizeof(double));
  double* EPS=malloc((N)*sizeof(double));
  double* C=malloc((N)*sizeof(double));

  double* Ustar=malloc((N+1)*sizeof(double));  // aux frontières
  double* Pstar=malloc((N+1)*sizeof(double));
  
  double* T=malloc((N)*sizeof(double));  // température 
  double* FM=malloc((N)*sizeof(double)); // fraction massique

  double** matFM=alloctabd(N,3);          // fraction massique
  double* lambda_out_unuse=malloc(3*sizeof(double));

  
	// Maillage
	for(int i=0; i<=N; i++){
	  X[i]=(a-nbghosts*dx)+dx*i;
    //printf("i= %d, X=%lf\n",i,X[i]);
	}
  
  
	// Initialisation 
	double x;
  // variables défini aux centres des mailles

  if(etain){
  	double* W=malloc(4*sizeof(double));
		for(int i=0; i<=N-1; i++){
			//printf("  --i= %d\n", i);
		  x=(X[i]+X[i+1])/2.0;

		  err=u0_etain(tst,a,b, x, W);
		  if(err){return err;}
		  
		  RHO0c[i]=W[0];
		  TAU[i]=1./W[0]; 
		  U[i]=W[1];
		  T[i]=W[3];
		  
		  err=epsilon_etain(1./W[0], W[2], T[i], &EPS[i], nb_points, SA, SB, SC, SAB, SAC, SBC );
		  if(err){return err;}

		  E[i]=W[1]*W[1]/2 + EPS[i];
		  M[i]=(X[i+1]-X[i])/TAU[i];
      /*
		  err=fPTC(TAU[i], EPS[i], &P[i], &C[i], &T[i], &FM[i],
		           nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, 
		           SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
      */
      //printf("   INIT\n");
		  err=fPTC_cin_phase(0, 0, sch_cin_phase, TAU[i], EPS[i], &P[i], &C[i], &T[i], matFM[i], matFM[i], Ntau, dx, dt, Ttau, &zone,
                   nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, 
		               SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
    
		  if(err){return err;}

		  XsC[i]=(X[i+1]-X[i])/C[i]; 
		}
		free(W);

  }
  else{
  	double* W=malloc(3*sizeof(double));
	  for(int i=0; i<=N-1; i++){
	    x=(X[i]+X[i+1])/2.0;
		  err=u0(tst,a,b, x, W);
		  if(err){return err;}

		  RHO0c[i]=W[0];
		  TAU[i]=1./W[0]; 
		  U[i]=W[1];
		  EPS[i]= epsilonEOS(tst, 1./W[0], W[2]);  // energie interne
		  E[i]=W[1]*W[1]/2 + EPS[i];

		  M[i]=(X[i+1]-X[i])/TAU[i];
		  
	    EGS(tst,TAU[i], EPS[i], &P[i], &C[i]);
	    XsC[i]=(X[i+1]-X[i])/C[i]; 
		}
		free(W);
	}
  

  // minimum de XsC  (calcul de dt)
	minXsC=XsC[ideb];
	for (int j=ideb+1; j<=ifin; j++){
	  if (minXsC>XsC[j]){
  	  minXsC=XsC[j];
  	}
  }
	dt = dt_first* CFL * minXsC; // pas de temps
	//printf("minXsC= %.15f\n", minXsC);
	printf("nx= %d et dt= %g\n", nx, dt);

  // energie totale du système
  int nbiter=10*ceil(Tf/dt);
  //printf("nbiter= %d\n", nbiter);
  double* Etot=malloc(nbiter*sizeof(double));
  double* IMPUL=malloc(nbiter*sizeof(double)); // impulsion
  double* Mtot=malloc(nbiter*sizeof(double));  // masse totale
  double* VOL=malloc(nbiter*sizeof(double));   // volume total

  sum1=0;  sum3=0.; sum4=0.; sum5=0.;
  for(int i=ideb; i<=ifin; i++){
    sum1+=dx*RHO0c[i]*E[i];
    sum3+=dx*U[i];
    sum4+=(X[i+1]-X[i])*TAU[i];
    //sum4+=dx*phi(so,&RHO0c[i], Ckbar);
    sum5+=dx*RHO0c[i]*TAU[i];
  }
  Etot[0]=sum1;
  //IMPUL[n+1]=sum3;
  IMPUL[0]=sum3; 
  if (tst==0){
    IMPUL[0]-=0.9*t; 
  }
  Mtot[0]=sum4;
  VOL[0]=sum5;
  

// debut de la boucle ============================================================================= 
  while((t<Tf)&&(n<Nmax)){

	  if(sch == 0){ // Despres
	  //printf("Despres\n");
	    for(int i=ideb; i<=ifin+1; i++){
	   	  rhoc=(C[i-1]/TAU[i-1]+ C[i]/TAU[i])/2; 
		    Ustar[i]=(U[i-1]+U[i])/2 + (P[i-1]-P[i])/(2*rhoc) ;
		    Pstar[i]=(P[i]+P[i-1])/2 + (rhoc/2)*(U[i-1]-U[i]) ;
		    //printf("  i=%d, Ustar=%g, Pstar=%g, rhoc=%g, TAU= %g , C= %g\n",i, Ustar[i],Pstar[i],rhoc,TAU[i],C[i]);
	    }
	  }	
	  else if (sch == 1){ // Jaouen 2001
	  	//printf("Jaouen\n");
			for(int i=ideb; i<=ifin+1; i++){
			  rhoc2=C[i-1]*C[i-1]/TAU[i-1]; 
			  rhoc21=C[i]*C[i]/TAU[i];
			  ms=fmax(rhoc21,rhoc2);
			  mi=fmax(TAU[i],TAU[i-1]);
			  rhoc=sqrt(ms/mi);
			  Ustar[i]=(U[i-1]+U[i])/2 + (P[i-1]-P[i])/(2*rhoc) ;
			  Pstar[i]=(P[i]+P[i-1])/2 + (rhoc/2)*(U[i-1]-U[i]) ;
			  //printf("i=%d, Ustar=%f, Pstar=%f, rhoc=%f, rhoc1=%f P[i-1]=%f, P[i]=%f\n",i, Ustar[i],Pstar[i],rhoc,rhoc1,P[i-1], P[i] );
			}
	  }
	  else if (sch == 2){ // Solveur Acoustique aux faces ordre 1
	  	//printf("Sol acous\n");
			for(int i=ideb; i<=ifin+1; i++){
			  rhoc=C[i-1]/TAU[i-1]; 
			  rhoc1=C[i]/TAU[i];
			  Ustar[i]=(rhoc*U[i-1]+rhoc1*U[i])/(rhoc+rhoc1) + (P[i-1]-P[i])/(rhoc+rhoc1) ;
			  Pstar[i]=(P[i]*rhoc+P[i-1]*rhoc1)/(rhoc+rhoc1) + (rhoc*rhoc1/(rhoc+rhoc1))*(U[i-1]-U[i]) ;
			  //printf("i=%d, Ustar=%f, Pstar=%f, rhoc=%f, rhoc1=%f P[i-1]=%f, P[i]=%f\n",i, Ustar[i],Pstar[i],rhoc,rhoc1,P[i-1], P[i] );
			}
	  }
	  else if ( (sch == 10) || (sch == 11) || (sch == 12) ){ // Solveur Acoustique aux faces d'odre 2
	  	//printf("Sol acous odre 2\n");
			for(int i=ideb; i<=ifin+1; i++){
			  rhoc=C[i-1]/TAU[i-1]; 
			  rhoc1=C[i]/TAU[i];
			  // 1er ordre
			  Ustar[i]=(rhoc*U[i-1]+rhoc1*U[i])/(rhoc+rhoc1) + (P[i-1]-P[i])/(rhoc+rhoc1) ;
			  Pstar[i]=(P[i]*rhoc+P[i-1]*rhoc1)/(rhoc+rhoc1) + (rhoc*rhoc1/(rhoc+rhoc1))*(U[i-1]-U[i]) ;
			  // 2nd ordre
	      mu= 0.5*(rhoc1+rhoc)*dt/(0.5*(M[i-1]+M[i]));

			  if (sch == 10 ){  // ordre 2 sans limiteur
			    phimoinsU=1.0;
			    phimoinsP=1.0;
			    phiplusU=1.0;
			    phiplusP=1.0;
			  }
			  else if( sch == 11 ){ // ordre 2 avec limiteur minmod
			    phimoinsU=MinMod((Ustar[i+1]-U[i])/(Ustar[i]-U[i-1]));
			    phimoinsP=MinMod((Pstar[i+1]-P[i])/(Pstar[i]-P[i-1]));

			    phiplusU=MinMod((U[i-1]-Ustar[i-1])/(U[i]-Ustar[i]));
			    phiplusP=MinMod((P[i-1]-Pstar[i-1])/(P[i]-Pstar[i]));
			  }
			  else if( sch == 12 ){ // ordre 2 avec limiteur superbee
			    // printf("superbee\n");
			    phimoinsU=Superbee((Ustar[i+1]-U[i])/(Ustar[i]-U[i-1]));
			    phimoinsP=Superbee((Pstar[i+1]-P[i])/(Pstar[i]-P[i-1]));

			    phiplusU=Superbee((U[i-1]-Ustar[i-1])/(U[i]-Ustar[i]));
			    phiplusP=Superbee((P[i-1]-Pstar[i-1])/(P[i]-Pstar[i]));
			  }
			  Ustar[i]+=(0.5)*(1.0-mu)*( phiplusU*(U[i]-Ustar[i]) - phimoinsU*(Ustar[i]-U[i-1]) );
			  Pstar[i]+=(0.5)*(1.0-mu)*( phiplusP*(P[i]-Pstar[i]) - phimoinsP*(Pstar[i]-P[i-1]) );
			  //printf("i=%d, Ustar=%f, Pstar=%f, rhoc=%f, rhoc1=%f P[i-1]=%f, P[i]=%f\n",i, Ustar[i],Pstar[i],rhoc,rhoc1,P[i-1], P[i] );
		    //printf("delta Ustar=%lf\n",(1.0/2)*(1-mu)*( (U[i][1]-Ustar[i]) - (Ustar[i]-U[i-1][1]) ) );
		    //printf("delta Pstar=%lf\n",(1.0/2)*(1-mu)*( (P[i]-Pstar[i]) - (Pstar[i]-P[i-1]) ) );
	    }
	  }
	  else{printf("Erreur dans le choix de Ustar et Pstar (Gtot) sch=%d\n",sch); return 1;}

    // Calcul du pas de temps
    minXsC=XsC[ideb];
	  for(int j=ideb+1; j<=ifin; j++){
	    if(minXsC>XsC[j]){
	  	  minXsC=XsC[j];
	  	}
	  }
	  dt_new = CFL * minXsC; // pas de temps	  
	  dt=pastemps(Tf,t,dt_new,dt);
      
      
    //dt = CFL * minXsC;
    //if ((T-t)<dt){ dt=T-t; }

	  // Mise a jour de la sol U
	  for(int i=ideb; i<=ifin; i++){
	    TAU[i]+=dt*(Ustar[i+1]-Ustar[i])/M[i];
	    U[i]-=dt*(Pstar[i+1]-Pstar[i])/M[i];
	    E[i]-=dt*(Ustar[i+1]*Pstar[i+1]-Ustar[i]*Pstar[i])/M[i];
	    //printf(" --i= %d   TAU= %g , U= %g , E= %g\n",i,TAU[i],U[i],E[i] );
	  }
	
    
	  // Mise à jour du maillage
	  for(int i=ideb; i<=ifin+1; i++){
	  	X[i]+=Ustar[i]*dt;
	  }

	  // cond aux limites de type mur (reflexion)
	  // cond aux limites     
    err=condlim_Gtot(ind_cond, ideb,ifin, nbghosts, TAU, 0);
    if(err){return err;}
    err=condlim_Gtot(ind_cond, ideb,ifin, nbghosts, E, 0);
    if(err){return err;}
    err=condlim_Gtot(ind_cond, ideb,ifin, nbghosts, U, iu);
	  if(err){return err;}

    // Mise à jour des tableaux
    if(etain){
	    for(int i=0; i<=N-1; i++){ 
	    	dX=X[i+1]-X[i];
	      EPS[i]= E[i]-U[i]*U[i]/2;  // energie interne 
	      /*
	      err=fPTC(TAU[i], EPS[i], &P[i], &C[i], &T[i], &FM[i],
	               nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, 
	               SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
        */
	      err=fPTC_cin_phase(0, cin_phase, sch_cin_phase, TAU[i], EPS[i], &P[i], &C[i], &T[i], matFM[i], matFM[i], Ntau, dX, dt, Ttau, &zone,
                   nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, 
		               SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
	      if(err){return err;}

	      XsC[i]=(X[i+1]-X[i])/C[i];          // dx/c  
	      M[i]=(X[i+1]-X[i])/TAU[i];
	      //printf(" --i= %d   TAU= %g , EPS= %g , P= %g , T= %g , C= %g\n",i,TAU[i],EPS[i],P[i],T[i],C[i] );
	    }
	  }
	  else{
	  	for(int i=0; i<=N-1; i++){ 
	      EPS[i]= E[i]-U[i]*U[i]/2;  // energie interne   
	      EGS(tst,TAU[i], EPS[i], &P[i], &C[i]);
	      XsC[i]=(X[i+1]-X[i])/C[i];          // dx/c  
	      M[i]=(X[i+1]-X[i])/TAU[i];
	    }
	  }


    sum1=0;  sum3=0.; sum4=0.; sum5=0.;
    for(int i=ideb; i<=ifin; i++){
      sum1+=dx*RHO0c[i]*E[i];
      sum3+=dx*U[i];
	    sum4+=(X[i+1]-X[i])*TAU[i];
	    //sum4+=dx*phi(so,&RHO0c[i], Ckbar);
	    sum5+=dx*RHO0c[i]*TAU[i];
	  }
	  Etot[n+1]=sum1;
	  //IMPUL[n+1]=sum3;
	  IMPUL[n+1]=sum3; 
	  if (tst==0){
	    IMPUL[n+1]-=0.9*t; 
	  }
	  Mtot[n+1]=sum4;
	  VOL[n+1]=sum5;

    /*
    if (n==1000){
    for(int i=0; i<=N-1; i++){
       printf("i=%d, TAU= %lf, U= %lf, E= %lf, P= %lf\n",i,TAU[i],U[i],E[i],P[i]);
    }
    for(int i=0; i<=N; i++){
      printf("i=%d, Ustar= %lf, Pstar=%lf\n",i,Ustar[i],Pstar[i]);
    } 
    }
    */
      
      
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
    Xc[i]=(X[i+1]+X[i])/2;
  }


  // printf_sol
  err=print_sol(1, ideb, ifin, X, Xc, TAU, U, E, P, EPS);
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

  
  free(U); free(TAU); free(E); 
  free(X); free(Xc);
  free(Ustar);    free(Pstar);
  free(P);  free(C); 
  free(EPS);
  free(XsC);

  freetab(matFM); free(T); free(lambda_out_unuse);
  
  free(Etot); 
  free(IMPUL);
  free(Mtot);
  free(VOL);

  return 0;
}
