#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "head2.h"

// Modification de la montée en ordre spatiale


// SCHEMA DE TYPE RUNGE-KUTTA FORMULE EN ENERGIE INTERNE

/////////////////////  MAIN  ////////////////////
int main(int argc, char* argv[]){
  double T;
	int nx=200;
	double CFL=0.3; 
	double a, b;
	int Nmax=9e8;
	int tst=0; 
	int sch=0;
  int z=0;
  int q=0;
  int so=2; // spatial order
  int dpi=0; //choix du dPI pour l'intégration de rho.u  
  double Cq=1.5; //--> artificial vicosity q
  double Cl=0.15;


  double* W=malloc(3*sizeof(double));


	// variables
	double t=0;
	double dt;
	int n=0;
  double minXsC; 
	char* W0;
  double sum1, sum2, sum3, sum4, sum5;
  int N;
  int nbghosts=20; // nombre de mailles fantomes pour les variables du schéma (\rho, \tau, u, \epsilon, ekin)
  int nbghosts2=10; // nombre de mailles fantomes pour les fonctions qui dependent des variables ci-dessus
  // l'idee  : appliquer les ocnd aux lim sur les var. sur nbghosts mailles qui calculer la valeur des fcts sur les nbghosts2 mailles
  int ideb, ifin;
  double** Ck; double** Ckbar; double** dk; double** Qbar; double** rk;
  int iu, ie;
  int sodPi; // spatial order for dPi

	for (int i = 0; i < argc; i++){
	  if (argv[i][0] == '-'){
	  	switch(argv[i][1]){
	  	  case 'h' : printf("Usage : \n"); 
	  	  			     printf(" -t test problem (0 : Sod , 1 : Bizarrium)\n");

                   printf(" -s hydro scheme Runge-Kutta energie interne\n");
                   printf("        0: RK1 (Matsuno forward-backward)\n");
                   printf("        1: RK2 (Heun's method)\n");
                   printf("        2: RK3 SSP\n");
                   printf("        3: RK4 Classique\n");
                   printf("        4: Cash-Karp ordre 5\n");
                   printf("        5: Dormand-Prince ordre 5\n");
                   printf("        6: HEUN  ordre 3\n");
                   printf("        7: RALSTON  ordre 3\n");
                   printf("        8: Bogaki-Shampine ordre 3\n");
                   printf("        9: Bogaki-Shampine ordre 2 \n");
                   printf("       10: KUTTA odre 3\n");
                   printf("       11: KUTTA ordre 2 (Explicit Midpoint Method)\n");
                   printf("       12: KUTTA ordre 2 Ralston\n");
                   printf("       13: EULER forward\n");

                   printf(" -q pseudo viscosity \n");
                   printf("      Ordre 1 : \n");
                   printf("        00 : von Neumann-Ritchmyer\n");
                   printf("        01 : Rosenbluth\n");
                   printf("        02 : Landschoff\n");
                   printf("        03 : Magical q combination\n");
                   printf("        04 : off\n");
                   printf("      Ordre 2 : \n");
                   printf("        10 : von Neumann-Ritchmyer\n");
                   printf("        11 : Rosenbluth\n");
                   printf("        12 : Landschoff\n");
                   printf("        13 : Magical q combination\n");
                   printf("        14 : off\n");

	  	  		    	 printf(" -o spatial order \n");
                   printf("        2 : 2nd order\n");
                   printf("        3 : 3rd order\n");
                   printf("        4 : 4th and 5th order\n");
                   printf("        5 : 6th and 7th order\n");
                   printf("        6 : 8th and 9th order\n");
                   
                   printf(" -n number cells \n");
	  	  	    		 
                   printf(" -f final time \n");
	  	       			 
                   printf(" -c final cycle \n");

                   printf(" -z kinétic energy fix\n");
                   printf("        0 : off\n");
                   printf("        1 : CRAS 2016\n");
                   printf("        2 : JCP Dakin et al. 2019\n");
                   printf("      A chaque sous-pas :\n");
                   printf("        11 : CRAS 2016\n");
                   printf("        21 : JCP Dakin et al. 2019\n");
                   
                   printf(" -p options\n");
                   printf("        0 : dPI\n");
                   printf("        1 : AV(delta PI)\n");
                   
                   printf(" -d options\n");
                   printf("        valeur du coefficent pseudo-quadratique\n");
                   
                   printf(" -l options\n");
                   printf("        valeur du coefficent pseudo-lineaire\n");
	  	  	    		 return 0;
	  	  case 't' : tst=atoi(argv[++i]);
	  	  			     break;
	  	  case 's' : sch=atoi(argv[++i]);
	  	  			     break;
        case 'o' : so=atoi(argv[++i]);
                   break;
	  	  case 'n' : nx=atoi(argv[++i]);
	  	  			     break;
	  	  case 'f' : T=strtod(argv[++i],NULL);
	  	  			     break;
	  	  case 'c' : Nmax=atoi(argv[++i]);
	  	  			     break;
        case 'z' : z=atoi(argv[++i]);
                   break;
        case 'p' : dpi=atoi(argv[++i]);
                   break;
        case 'q' : q=atoi(argv[++i]);
                   break;
        case 'd' : Cq=strtod(argv[++i],NULL);
                   break;
        case 'l' : Cl=strtod(argv[++i],NULL);
                   break;
	  	  default :  printf("Not recognized option : %c\n",argv[i][1] );
	  	  			     return 1;
	  	}
	  }
	}

  if (tst==0){
    a=1.;  b=2.0;  ////////////////////////////////////////////           !!!!!!
    T=0.2;
  }
  else if (tst==1){
    a=0.;  b=1.0;  
    T=80e-6; 
  }
  else if(tst==2){
    a=0.;  b=1.0;  
    T=1.0;
  }
  else if(tst==3){
    a=-1.0;  b=1.0;  
    //T=0.2;
  }

	printf("sch= %d\n",sch );
  printf("tst= %d\n",tst );
  printf("q= %d\n",q );

  // params pour les cond aux limites
  // changement de signe pour u si iu==1
  if (tst==0|| tst==3){
    iu=1;  // indice pour la vitesse
    ie=0;  // indice pour l'energie
  }
  else if (tst==1){
    iu=0; 
    ie=0; 
  }
  else if(tst==2){
    iu=0;
    ie=0; 
  }
  
  // spatial order for dPi
  sodPi=so;


  // Coefficient
  Ck=alloctabd(5,5); Ckbar=alloctabd(5,5);
  dk=alloctabd(5,5); Qbar=alloctabd(5,5);
  rk=alloctabd(5,5);
  initCoef(Ck, Ckbar, dk, Qbar, rk);
  
  // Coefficient de Btucher
  int to; // time order
  double** A; double* THETA; double* ALPHA;
  // RK1 et RK2
  if(sch==0||sch==1){   to=1;  }
  // SSPRK3(3,3) Spiteri-Ruuth
  else if (sch==2){    to=2;  }
  // RK4 Classic
  else if (sch==3){    to=3;  }
  // RK5 Cash-Karp 
  else if (sch==4){    to=5;  }
  // Dormand-Prince 
  else if (sch==5){    to=6;  }
  // HEUN  ordre 3
  else if (sch==6){    to=2;  }
  // RK3 RALSTON
  else if (sch==7){    to=2;  }
  // RK2 ou RK3 Bogaki-Shampine
  else if (sch==8 || sch==9){    to=3;  }
  // KUTTA odre 3
  else if (sch==10){    to=2;  }
  // KUTTA ordre 2 (Explicit Midpoint Method) 
  else if (sch==11){    to=1;  }
  // KUTTA ordre 2 Ralston
  else if (sch==12){    to=1;  }
  // EULER forward
  else if (sch==13){    to=1;  }
  // SSP RK4 OPTIMUM
  else if (sch==14){    to=4;  }
  // SSPRK3(4,3) Spiteri-Ruuth
  else if (sch==15){    to=3;  }
  // SSPRK3(5,3) Spiteri-Ruuth
  else if (sch==16){    to=4;  }
  // SPPRK1(3,1) Spiteri-Ruuth
  else if (sch==17){    to=2;  }
  // SPPRK2(3,2) Spiteri-Ruuth
  else if (sch==18){    to=2;  }
  // SPPRK2(4,2) Spiteri-Ruuth
  else if (sch==19){    to=3;  }
  // SPPRK1(2,1) Spiteri-Ruuth
  else if (sch==20){    to=1;  }
  // SPPRK3(4,3) Spiteri-Ruuth
  else if (sch==21){    to=3;  }
  else{ printf("Erreur dans le choix de sch (main)\n");}
  

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
  double* TAUpw=malloc((N)*sizeof(double)); //  rho0.tau  (centrées aux mailles)
  double* Upw=malloc((N+1)*sizeof(double)); //  rho0.u (frontières)
  double* EPSpw=malloc((N)*sizeof(double)); //  rho0.epsilon (centrées aux mailles)

  // AVERAGE (plot)
  double* TAUav=malloc((N)*sizeof(double)); //  tau  (centrées aux mailles)
  double* EPSav=malloc((N)*sizeof(double)); //  epsilon (centrées aux mailles)  
  double* Uav=malloc((N+1)*sizeof(double)); //  u  (frontières aux mailles)
  double* Pav=malloc((N)*sizeof(double));   //  p (centrées aux mailles)  

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
    u0(tst,a,b, x, W);
    RHO0c[i]=W[0];
    REPSpw[i]=W[0]*epsilonEOS(tst, 1./W[0], W[2]);
    //printf("rho=%lf, p=%lf\n",W[0], W[2] );
	  //printf("i= %d, rho0.tau=%f, rho0.eps=%f, rho0c=%f\n",i,RTAUpw[i],REPSpw[i],RHO0c[i] );
	}
  // variables défini aux frontières des mailles
  for(int i=0; i<=N; i++){
    x=X[i];
    u0(tst,a,b, x, W);
    RHO0d[i]=W[0];
    RUpw[i]=W[0]*W[1];
    REKINpw[i]=W[0]*W[1]*W[1]/2.0;
    //printf("i= %d,  RUpw=%f,  REKINpw=%f\n",i,RUpw[i],REKINpw[i] );
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
    REPS[i]=phibar(so,0, &REPSpw[i], &X[i], rk);
    RU[i]=phibar(so,1, &RUpw[i], &X[i], rk);
    REKIN[i]=phibar(so,1, &REKINpw[i], &X[i], rk);
    //printf("i=%d  REPS=%lf, RU=%lf, REKIN=%lf\n",i,REPS[i],RU[i],REKIN[i]);           
  }
  RU[ifin+1]=phibar(so,1, &RUpw[ifin+1], &X[ifin+1], rk);
  REKIN[ifin+1]=phibar(so,1, &REKINpw[ifin+1], &X[ifin+1], rk);
  //printf("i=ifin+1   RU=%lf, REKIN=%lf\n",RU[ifin+1],REKIN[ifin+1]);    
  
  
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
  printf("nbiter= %d\n",nbiter );
  double* Eint=malloc(nbiter*sizeof(double));
  double* Ekin=malloc(nbiter*sizeof(double));
  double* IMPUL=malloc(nbiter*sizeof(double)); // impulsion
  double* Mtot=malloc(nbiter*sizeof(double)); // masse totale
  double* VOL=malloc(nbiter*sizeof(double)); // volume total
  
  // Energie du système
  for(int i=0; i<=N; i++){
    R2U2pw[i]=RUpw[i]*RUpw[i];
    RU2pw[i]=RUpw[i]*RUpw[i]/RHO0d[i];
  }

  sum1=0; sum2=0; sum3=0.; sum4=0.; sum5=0.;
  for(int i=ideb; i<=ifin; i++){
    sum1+=dx*REPS[i];
    sum2+=(dx/2.0)*phibar(so,1, &RU2pw[i], &X[i], rk);
    sum3+=dx*RU[i];
    sum4+=dx*phibar(so,0,&RHO0c[i],&X[i], rk);
    sum5+=dx*RTAU[i];
  }
  sum2+=(dx/2.0)*phibar(so,1, &RU2pw[ifin+1], &X[ifin+1], rk);
  sum3+=dx*RU[ifin+1];
  //sum4+=dx*phi(so,&RHO0d[ifin+1],Ckbar);
  Eint[0]=sum1;
  Ekin[0]=sum2;
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
    condlim(tst, ideb,ifin, nbghosts, RTAU, 0);
    condlim(tst, ideb,ifin, nbghosts, REPS, 0);
    condlim(tst, ideb,ifin+1, nbghosts, RU, iu);
    condlim(tst, ideb,ifin+1, nbghosts, REKIN, 0); // attention frontières aux mailles
    
    // Calcul du pas de temps
    minXsC=XsC[ideb];
    for(int j=ideb+1; j<=ifin; j++){
      if (minXsC>XsC[j]){
        minXsC=XsC[j];
      }
    }
    dt = CFL * minXsC; // pas de temps
    if ((T-t)<dt){ dt=T-t; }

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
        else if(dpi==1){  RUstar[k][i]=RU[i]-A[k][0]*(dt/dx)*phibar(sodPi,1, &dPI[i], &X[i], rk);  }
        REPSstar[k][i]=REPS[i]-A[k][0]*(dt/dx)*YdZ(so,0,&P[i],&Q[i],&U[i],&X[i],rk);
        REKINstar[k][i]=REKIN[i]-A[k][0]*(dt/dx)*YdZ(so,1,&P[i],&Q[i],&U[i],&X[i],rk);
      }
      if(dpi==0){  RUstar[k][ifin+1]=RU[ifin+1]-A[k][0]*(dt/dx)*((P[ifin+1]+Q[ifin+1])-(P[ifin]+Q[ifin]));  }
      else if(dpi==1){  RUstar[k][ifin+1]=RU[ifin+1]-A[k][0]*(dt/dx)*phibar(sodPi,1, &dPI[ifin+1], &X[ifin+1], rk);  }
      REKINstar[k][ifin+1]=REKIN[ifin+1]-A[k][0]*(dt/dx)*YdZ(so,1,&P[ifin+1],&Q[ifin+1],&U[ifin+1],&X[ifin+1],rk); 

      for(int i=nbghosts - nbghosts2; i<=N2; i++){
        Xstar[k][i]= X[i] + A[k][0]*dt*U[i];
      }
      

      for(int j=1; j<=k; j++){ // calcul des sous pas de temps
        //
        for(int i=ideb; i<=ifin; i++){
    	    RTAUstar[k][i]+=A[k][j]*(dt/dx)*(Ustar[j-1][i+1]-Ustar[j-1][i]);
          if(dpi==0){   RUstar[k][i]-=A[k][j]*(dt/dx)*((Pstar[j-1][i]+Qstar[j-1][i])-(Pstar[j-1][i-1]+Qstar[j-1][i-1]));   }
          else if(dpi==1){    RUstar[k][i]-=A[k][j]*(dt/dx)*phibar(sodPi,1, &dPIstar[j-1][i], &Xstar[k][i], rk);   }
          REPSstar[k][i]-=A[k][j]*(dt/dx)*YdZ(so,0,&Pstar[j-1][i],&Qstar[j-1][i],&Ustar[j-1][i],&Xstar[k][i], rk);
          REKINstar[k][i]-=A[k][j]*(dt/dx)*YdZ(so,1,&Pstar[j-1][i],&Qstar[j-1][i],&Ustar[j-1][i],&Xstar[k][i], rk); 
    	  }
        if(dpi==0){   RUstar[k][ifin+1]-=A[k][j]*(dt/dx)*((Pstar[j-1][ifin+1]+Qstar[j-1][ifin+1])-(Pstar[j-1][ifin]+Qstar[j-1][ifin]));  }
        else if(dpi==1){   RUstar[k][ifin+1]-=A[k][j]*(dt/dx)*phibar(sodPi,1, &dPIstar[j-1][ifin+1], &Xstar[k][ifin+1], rk);  }
        REKINstar[k][ifin+1]-=A[k][j]*(dt/dx)*YdZ(so,1,&Pstar[j-1][ifin+1],&Qstar[j-1][ifin+1],&Ustar[j-1][ifin+1],&Xstar[k][ifin+1], rk);

        for(int i=nbghosts - nbghosts2; i<=N2; i++){
          Xstar[k][i]+=A[k][j]*dt*Ustar[j-1][i];
        }
      }
      // cond aux limites 

      condlim(tst, ideb,ifin, nbghosts, RTAUstar[k], 0);
      condlim(tst, ideb,ifin, nbghosts, REPSstar[k], 0);
      condlim(tst, ideb,ifin+1, nbghosts, RUstar[k], iu);
      condlim(tst, ideb,ifin+1, nbghosts, REKINstar[k], 0);

      // Energy kinetic fix
      // CRAS 2016
      if(z==11 || z==12){ 
        for(int i=nbghosts-nbghosts2; i<=N2; i++){ 
          RUpw[i]=phibar(so,1, &RUstar[k][i],&Xstar[k][i],rk);
        }

        if(z==11){
          for(int i=nbghosts-nbghosts2; i<=N2; i++){ ////  
            REKINpw[i]=phipw(so,1, &REKINstar[k][i],&Xstar[k][i],rk);
          }
          for (int i=nbghosts-nbghosts2; i<=N2; i++){
            DK[i]=REKINpw[i]-RUpw[i]*RUpw[i]/(2.0*RHO0d[i]);
            //printf("i: %d, DK= %.8lf\n",i,DK[i] );
          }
          for (int i=ideb; i<=ifin; i++){
            REKINstar[k][i]-=phibar(so,1, &DK[i],&Xstar[k][i],rk);
            REPSstar[k][i]+=phiQK(so, &DK[i],&Xstar[k][i],rk); 
          }
          REKINstar[k][ifin+1]-=phibar(so,1,&DK[ifin+1],&Xstar[k][ifin+1],rk);
        }
        // JCP 2019 Dakin & al.
        else if(z==12){
          for(int i=nbghosts-nbghosts2; i<=N2; i++){ ////
            RU2pw[i]=RUpw[i]*RUpw[i]/(2.0*RHO0d[i]);
            //printf("i: %d, RU2pw= %.8lf\n",i,RU2pw[i] );
          }
          for (int i=ideb; i<=ifin+1; i++){
            DK[i]=REKINstar[k][i]-phibar(so,1, &RU2pw[i], &Xstar[k][i], rk);
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
        condlim(tst, ideb,ifin, nbghosts, REPSstar[k], 0);
        condlim(tst, ideb,ifin, nbghosts, REKINstar[k], 0);
        
      }
      

      // On veut que U soit PW or PUstar est AV
      for(int i=nbghosts-nbghosts2; i<=N2-1; i++){ ////
        Ustar[k][i]=phipw(so,1, &RUstar[k][i],&Xstar[k][i], rk)/RHO0d[i]; 
        TAUstar[k][i]=phipw(so,0, &RTAUstar[k][i],&Xstar[k][i], rk)/RHO0c[i]; 
        EPSstar[k][i]=phipw(so,0, &REPSstar[k][i], &Xstar[k][i], rk)/RHO0c[i]; 
      }
      Ustar[k][N2]=phipw(so,1, &RUstar[k][N2],&Xstar[k][N2], rk)/RHO0d[N2]; 
      
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
      else if(dpi==1){  RU[i]-=THETA[0]*(dt/dx)*phibar(sodPi,1, &dPI[i], &X[i], rk);  }
      REPS[i]-=THETA[0]*(dt/dx)*YdZ(so,0,&P[i],&Q[i],&U[i],&X[i],rk);
      REKIN[i]-=THETA[0]*(dt/dx)*YdZ(so,1,&P[i],&Q[i],&U[i],&X[i],rk);    
    }
    if(dpi==0){  RU[ifin+1]-=THETA[0]*(dt/dx)*((P[ifin+1]+Q[ifin+1])-(P[ifin]+Q[ifin]));  }
    else if(dpi==1){  RU[ifin+1]-=THETA[0]*(dt/dx)*phibar(sodPi,1, &dPI[ifin+1],&X[ifin+1],rk);  }
    REKIN[ifin+1]-=THETA[0]*(dt/dx)*YdZ(so,1,&P[ifin+1],&Q[ifin+1],&U[ifin+1],&X[ifin+1],rk);

    for(int i=nbghosts-nbghosts2; i<=N2; i++){
      X[i]+=dt*THETA[0]*U[i];
    }


    for (int j=1; j<=to; j++){ // calcul des sous pas de temps
      //
      for(int i=ideb; i<=ifin; i++){
        RTAU[i]+=THETA[j]*(dt/dx)*(Ustar[j-1][i+1]-Ustar[j-1][i]);
        if(dpi==0){   RU[i]-=THETA[j]*(dt/dx)*((Pstar[j-1][i]+Qstar[j-1][i])-(Pstar[j-1][i-1]+Qstar[j-1][i-1]));   }
        else if(dpi==1){    RU[i]-=THETA[j]*(dt/dx)*phibar(sodPi,1, &dPIstar[j-1][i],&Xstar[j-1][i],rk);   }
        REPS[i]-=THETA[j]*(dt/dx)*YdZ(so,0,&Pstar[j-1][i],&Qstar[j-1][i],&Ustar[j-1][i],&Xstar[j-1][i],rk); 
        REKIN[i]-=THETA[j]*(dt/dx)*YdZ(so,1,&Pstar[j-1][i],&Qstar[j-1][i],&Ustar[j-1][i],&Xstar[j-1][i],rk);  
      }
      if(dpi==0){   RU[ifin+1]-=THETA[j]*(dt/dx)*((Pstar[j-1][ifin+1]+Qstar[j-1][ifin+1])-(Pstar[j-1][ifin]+Qstar[j-1][ifin]));  }
      else if(dpi==1){   RU[ifin+1]-=THETA[j]*(dt/dx)*phibar(sodPi,1, &dPIstar[j-1][ifin+1],&Xstar[j-1][ifin+1],rk);  }
      REKIN[ifin+1]-=THETA[j]*(dt/dx)*YdZ(so,1,&Pstar[j-1][ifin+1],&Qstar[j-1][ifin+1],&Ustar[j-1][ifin+1],&Xstar[j-1][ifin+1],rk);  
      
      for(int i=nbghosts-nbghosts2; i<=N2; i++){ // permet la gestion des cond aux lim du maillage
        X[i]+=dt*THETA[j]*Ustar[j-1][i];
      }

    }

  
    // cond aux limites  symétrie
    condlim(tst, ideb,ifin, nbghosts, RTAU, 0);
    condlim(tst, ideb,ifin, nbghosts, REPS, 0);
    condlim(tst, ideb,ifin+1, nbghosts, RU, iu);
    condlim(tst, ideb,ifin+1, nbghosts, REKIN, iu); // ATENTION frontières aux mailles


    // Mise a jour des rho.Z  Point-Wise
    for(int i=nbghosts-nbghosts2; i<=N2-1; i++){ //// 
      RTAUpw[i]=phipw(so,0, &RTAU[i], &X[i],rk);
      REPSpw[i]=phipw(so,0, &REPS[i], &X[i],rk);
      RUpw[i]=phipw(so,1, &RU[i], &X[i],rk);
    }
    RUpw[N2]=phipw(so,1, &RU[N2], &X[N2],rk);


    // Energy kinetic fix
    // CRAS 2016
    if(z==1 || z==2 || z==11 || z==12){ 
      if(z==1 || z==11){

        for(int i=nbghosts-nbghosts2; i<=N2; i++){ ////  
          REKINpw[i]=phipw(so,1, &REKIN[i], &X[i], rk);
        }
        for (int i=nbghosts-nbghosts2; i<=N2; i++){
          DK[i]=REKINpw[i]-RUpw[i]*RUpw[i]/(2.0*RHO0d[i]);
          //printf("i: %d, DK= %.8lf\n",i,DK[i] );
        }  
        for (int i=ideb; i<=ifin; i++){
          REKIN[i]-=phibar(so,1, &DK[i],&X[i],rk);
          REPS[i]+=phiQK(so, &DK[i],&X[i],rk); 
        }
        REKIN[ifin+1]-=phibar(so,1,&DK[ifin+1],&X[ifin+1],rk);
      }
      // JCP 2019 Dakin & al.
      else if(z==2 || z==12){
        for(int i=nbghosts-nbghosts2; i<=N2; i++){ ////
          RU2pw[i]=RUpw[i]*RUpw[i]/(2.0*RHO0d[i]);
          //printf("i: %d, RU2pw= %.8lf\n",i,RU2pw[i] );
        }
        for (int i=ideb; i<=ifin+1; i++){
          DK[i]=REKIN[i]-phibar(so,1, &RU2pw[i], &X[i],rk);
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
      condlim(tst, ideb,ifin, nbghosts, REPS, 0);
      
      // Mise a jour des rho.Z  Point-Wise
      for(int i=nbghosts-nbghosts2; i<=N2-1; i++){  
        REPSpw[i]=phipw(so,0, &REPS[i], &X[i],rk);
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
    for(int i=nbghosts-nbghosts2; i<=N2-1; i++){
      R2U2pw[i]=RUpw[i]*RUpw[i];  // (rho.u)^2
      RU2pw[i]=RUpw[i]*RUpw[i]/RHO0d[i];
      TAUpw[i]=RTAUpw[i]/RHO0c[i];
    }
    R2U2pw[N2]=RUpw[N2]*RUpw[N2];  // (rho.u)^2
    RU2pw[N2]=RUpw[N2]*RUpw[N2]/RHO0d[N2];

    sum1=0; sum2=0; sum3=0.; sum4=0.; sum5=0.;
    for(int i=ideb; i<=ifin; i++){
      sum1+=dx*REPS[i];
      sum2+=(dx/2.0)*phibar(so,1, &RU2pw[i], &X[i], rk);
      sum3+=dx*RU[i];
      sum4+=(X[i+1]-X[i])/phibar(so,0, &TAUpw[i],&X[i], rk);
      //sum4+=dx*phi(so,&RHO0c[i], Ckbar);
      sum5+=dx*RTAU[i];
    }
    sum2+=(dx/2.0)*phibar(so,1, &RU2pw[ifin+1], &X[ifin+1], rk);
    sum3+=dx*RU[ifin+1];
    Eint[n+1]=sum1;
    Ekin[n+1]=sum2;
    //IMPUL[n+1]=sum3;
    IMPUL[n+1]=sum3; 
    if (tst==0){
      IMPUL[n+1]-=0.9*t; 
    }
    Mtot[n+1]=sum4;
    VOL[n+1]=sum5;
    

    /*
    for (int i = 1; i <= nx; i++){
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
	  printf("n= %d, t= %.15lf, dt= %.15lf\n",n,t,dt);
	}
//==========================================================================

  printf("Nombre d'itérations : %d\n",n);
  printf("Temps final : %.15lf\n",t);

  for(int i=nbghosts-nbghosts2; i<=N2-1; i++){ //// 
    TAUpw[i]=RTAUpw[i]/RHO0c[i];
    EPSpw[i]=REPSpw[i]/RHO0c[i];
    Upw[i]=RUpw[i]/RHO0d[i];
  }
  Upw[N2]=RUpw[N2]/RHO0d[N2];


  // Tableaux pour plot
  for (int i=ideb; i<=ifin; i++){
    TAUav[i]=phibar(so,0, &TAUpw[i], &X[i], rk);
    EPSav[i]=phibar(so,0, &EPSpw[i], &X[i], rk);
    Uav[i]=phibar(so,1, &Upw[i], &X[i], rk);
    Pav[i]=phibar(so,0, &P[i], &X[i], rk);
    //printf("i=%d, TAUav=%lf, TAUpw=%lf, RHO0c=%lf\n",i,TAUav[i],RTAUpw[i],RHO0c[i] );
  }
  Uav[ifin+1]=phibar(so,1, &RUpw[ifin+1], &X[ifin+1], rk);
  

  // ecriture dans un fichier des points (x,u(x)) de la solution
  FILE* Res;
  /*
  L1 : X  maillage primal
  L2 : X  maillage duale
  L3 : tau  volume spécifique
  L4 : u  vitesse
  L5 : e   energie interne
  L6 : p   pression
  L7 : epsilon   energie interne spécifique
  */

  if((Res = fopen("sol.txt", "w+")) != NULL){
    // X frontière
    for(int i=ideb; i<=ifin+1; i++){ 
      fprintf(Res,"%.15lf ",X[i]);
    }
    fprintf(Res, "\n");
    // X centrée
    for (int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15lf ", fxc(so, &X[i], rk) );
    }
    fprintf(Res, "1e20\n");
    fprintf(Res, "\n");
    // tau
    for(int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15lf ",TAUav[i]);
    }
    fprintf(Res, "1e20\n");
    fprintf(Res, "\n");
    // u
    for(int i=ideb; i<=ifin+1; i++){
      fprintf(Res,"%.15lf ",Uav[i]);
    }
    fprintf(Res, "\n");
    // e
    for(int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15lf ",EPSav[i]+((Uav[i]+Uav[i+1])/2)*((Uav[i]+Uav[i+1])/2)/2);
    }
    fprintf(Res, "1e20\n");
    fprintf(Res, "\n");
    // p
    for(int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15lf ",Pav[i]);
    }
    fprintf(Res, "1e20\n");
    fprintf(Res, "\n");
    // epsilon
    for(int i=ideb; i<=ifin; i++){ 
      fprintf(Res,"%.15lf ",EPSav[i]);
      // fprintf(Res,"%.15lf ",REPS[i]/RHO0cAV[i]);
    }
    fprintf(Res, "1e20\n");
    fclose(Res);
  } 
  else{printf("Erreur d'ouverture de sol.txt\n");}
  
  // Conservation au cours du temps
  if((Res = fopen("NRJ.txt", "w+")) != NULL){
    // energie interne
    for(int i=0; i<=n-1; i++){
      fprintf(Res,"%.15lf ",Eint[i]);
    }
    fprintf(Res, "\n");
    // energie cinétique
    for(int i=0; i<=n-1; i++){ 
      fprintf(Res,"%.15lf ",Ekin[i]);
    }
    fprintf(Res, "\n");
    // impulsion
    for(int i=0; i<=n-1; i++){ 
      fprintf(Res,"%.15lf ",IMPUL[i]);
    }
    fprintf(Res, "\n");
    // masse totale
    for(int i=0; i<=n-1; i++){ 
      fprintf(Res,"%.15lf ",Mtot[i]);
    }
    fprintf(Res, "\n");
    // volume total
    for(int i=0; i<=n-1; i++){ 
      fprintf(Res,"%.15lf ",VOL[i]);
    }
    fprintf(Res, "\n");
    fclose(Res);
  } 
  else{printf("Erreur d'ouverture de NRJ.txt\n");}
  

  
  free(dPI); freetab(dPIstar);

  //free(WG); free(WD);
  free(W);

  free(X); 
  free(RHO0d); free(RHO0c);
  
  free(RTAU); free(RU); free(REPS);
  free(RTAUpw); free(RUpw); free(REPSpw); 
  free(TAUpw); free(Upw); free(EPSpw); 

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
  free(Eint); free(Ekin);
  free(RU2pw);
  free(IMPUL);
  free(Mtot);
  free(VOL);
}
