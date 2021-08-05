#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "head2.h"


//  -->  Faire un cas pour les variables sur maille primal
//          et un cas pour les variables sur maille dual


// Z :  variable AVERAGE dont on veut la valeur ponctuelle (&Z[i])
// X :  maillage primal   (&X[i])
// ind : 0-centré  1-décentré
double phipw(int so, int ind, double* Z, double* X, double** rk){
  // entier utile :
  int o=so-2; //ordre de la méthode
  int n=2*(so-2)+2; // taille du sytème
  double akk, aik;
  double res=0.0;
  
  // La matrice et le second membre sont stockés sur le même tableau
  double** A=alloctabd(n,n+1);
 
  //////////////////////////////////////////////////////////
  // variable sur maille duale
  if (ind==0){
	// Contruction de la matrice de Cramer A
	// sur le carré de taille n dans A
	for(int i=0; i<=n-1; i++){
	  for(int j=0; j<=n-1; j++){
	    A[i][j]=1.0;
	    // A=x^?
	    for(int k=1; k<=n-1-j; k++){ // quand j=n-1, on entre pas dans la boucle A[i][n-1]=1
	      A[i][j]*=X[(-o)+i]; //?????
	    } 
      }
	}


    // Construction du second menbre 
    // sur la (n+1)ième colonne de A
    for(int j=0; j<=n-1; j++){
      A[j][n]=0.;
      for(int k=1; k<=j; k++){ // quand j=0, on entre pas dans la boucle Y[0]=0
	    A[j][n]+=Z[(-o)+k-1]*(X[(-o)+k]-X[(-o)+k-1]);
	  } 
	}


	// Résolution du sytème linèaire
	// ref fonction d'internet
	/* algorithme de Gauss-Jordan */
	for(int k=0;k<n;k++){
	  akk=A[k][k];
	  for(int j=k;j<n+1;j++){
	    A[k][j]=A[k][j]/akk;
	  }
	  for(int i=0;i<n;i++){
	    aik=A[i][k]; 
	    if(i!=k){
	      for(int j=k;j<n+1;j++){
	        A[i][j]=A[i][j]-aik*A[k][j];
	      }
	    }
	  }
	}

	// Calcul du point de calcul x
	double x=fxc(so, X, rk);

	// Calcul de la valeur ponctuelle
	double monome_x=1.0;
	for(int i=0; i<=n-2; i++){
	  res+=(i+1)*monome_x*A[n-2-i][n];
	  monome_x*=x;
	}
  }

  //////////////////////////////////////////////////////////
  // variable sur maille primale
  else if(ind==1){
  	double *Xc=malloc(n*sizeof(double));
  	for(int i=0; i<=n-1; i++){
  	  Xc[i]=fxc(so,&X[-(o+1)+i], rk); /// o <- o+1
  	}
    // Contruction de la matrice de Cramer A
	// sur le carré de taille n dans A
	for(int i=0; i<=n-1; i++){
	  for(int j=0; j<=n-1; j++){
	    A[i][j]=1.0;
	    // A=x^?
	    for(int k=1; k<=n-1-j; k++){ // quand j=n-1, on entre pas dans la boucle A[i][n-1]=1
	      A[i][j]*=Xc[i]; 
	    } 
      }
	}


    // Construction du second menbre 
    // sur la (n+1)ième colonne de A
    for(int j=0; j<=n-1; j++){
      A[j][n]=0.;
      for(int k=1; k<=j; k++){ // quand j=0, on entre pas dans la boucle Y[0]=0
	    A[j][n]+=Z[(-o)+k-1]*(Xc[k]-Xc[k-1]);
	  } 
	}


	// Résolution du sytème linèaire
	// ref fonction d'internet
	/* algorithme de Gauss-Jordan */
	for(int k=0;k<n;k++){
	  akk=A[k][k];
	  for(int j=k;j<n+1;j++){
	    A[k][j]=A[k][j]/akk;
	  }
	  for(int i=0;i<n;i++){
	    aik=A[i][k]; 
	    if(i!=k){
	      for(int j=k;j<n+1;j++){
	        A[i][j]=A[i][j]-aik*A[k][j];
	      }
	    }
	  }
	}
    

	// Calcul de la valeur ponctuelle
	double monome_x=1.0;
	for (int i=0; i<=n-2; i++){
	  res+=(i+1)*monome_x*A[n-2-i][n];
	  monome_x*=X[0];
	}


    free(Xc);
  }
  
  return res;
  
  freetab(A);
  
}


// Calcul la maille dual (centrée) à partir de la maille primal (frontière) 
double fxc(int so,double* X, double** rk){
  // voir CRAS 2016 p214
  double res=0.0;
  int k=so-2;
  for(int l=0; l<=k; l++){
    res+=rk[k][l]*(X[l+1] + X[-l]);
  }
  return res;
}


// Z :  variable POINT-WISE dont on veut la valeur moyenne (&Z[i])
// X :  maillage primal   (&X[i])
double phibar(int so, int ind, double* Z, double* X, double** rk){
  int o=so-2; //ordre de la méthode
  int n=2*(so-2)+1; // taille du sytème 
  // La taille du système n'est pas la même car on a besoin d'un poly d'ordre so-1 uniquement 
  // car il n'y a pas de dérivation contraiement a phipw
  double akk, aik;
  double res=0.0;

  // La matrice et le second membre sont stockée sur le même tableau
  double** A=alloctabd(n,n+1);

  //////////////////////////////////////////////////////////
  // variable sur maille duale
  if (ind==0){
  	double *Xc=malloc(n*sizeof(double));
  	for(int i=0; i<=n-1; i++){
  	  Xc[i]=fxc(so,&X[(-o)+i], rk);
  	}
	// Contruction de la matrice de Cramer A
	// sur le carré de taille n dans A
	for(int i=0; i<=n-1; i++){
	  for(int j=0; j<=n-1; j++){
	    A[i][j]=1.0;
	    // A=x^?
	    for(int k=1; k<=n-1-j; k++){ // quand j=n-1, on entre pas dans la boucle A[i][n-1]=1
	      A[i][j]*=Xc[i]; 
	    } 
      } 
	}


    // Construction du second menbre 
    // sur la (n+1)ième colonne de A
    for(int j=0; j<=n-1; j++){
	  A[j][n]=Z[(-o)+j];
	}


	// Résolution du sytème linèaire
	// ref fonction d'internet
	/* algorithme de Gauss-Jordan */
	for(int k=0;k<n;k++){
	  akk=A[k][k];
	  for(int j=k;j<n+1;j++){
	    A[k][j]=A[k][j]/akk;
	  }
	  for(int i=0;i<n;i++){
	    aik=A[i][k]; 
	    if(i!=k){
	      for(int j=k;j<n+1;j++){
	        A[i][j]=A[i][j]-aik*A[k][j];
	      }
	    }
	  }
	}


	// Calcul de la valeur moyenne
	double xp, xm;
	xp=X[1];
	xm=X[0];
	double monome_xp=xp, monome_xm=xm;
	// intégrale 
	for (int i=0; i<=n-1; i++){
	  res+=(monome_xp-monome_xm)*A[n-1-i][n]/(i+1);
	  monome_xp*=xp;
	  monome_xm*=xm;
	}
	res=res/(xp-xm);

	free(Xc);
  }
  //////////////////////////////////////////////////////////
  // variable sur maille primale
  if (ind==1){
  	// Contruction de la matrice de Cramer A
	// sur le carré de taille n dans A
	for(int i=0; i<=n-1; i++){
	  for(int j=0; j<=n-1; j++){
	    A[i][j]=1.0;
	    // A=x^?
	    for(int k=1; k<=n-1-j; k++){ // quand j=n-1, on entre pas dans la boucle A[i][n-1]=1
	      A[i][j]*=X[(-o)+i]; 
	    } 
      } 
	}


    // Construction du second menbre 
    // sur la (n+1)ième colonne de A
    for(int j=0; j<=n-1; j++){
	  A[j][n]=Z[(-o)+j];
	}


	// Résolution du sytème linèaire
	// ref fonction d'internet
	/* algorithme de Gauss-Jordan */
	for(int k=0;k<n;k++){
	  akk=A[k][k];
	  for(int j=k;j<n+1;j++){
	    A[k][j]=A[k][j]/akk;
	  }
	  for(int i=0;i<n;i++){
	    aik=A[i][k]; 
	    if(i!=k){
	      for(int j=k;j<n+1;j++){
	        A[i][j]=A[i][j]-aik*A[k][j];
	      }
	    }
	  }
	}


	// Calcul de la valeur moyenne
	double xp, xm;
	xp=fxc(so, &X[0], rk);
	xm=fxc(so, &X[-1], rk);
	double monome_xp=xp, monome_xm=xm;
	// intégrale
	for (int i=0; i<=n-1; i++){
	  res+=(monome_xp-monome_xm)*A[n-1-i][n]/(i+1);
	  monome_xp*=xp;
	  monome_xm*=xm;
	}
    res=res/(xp-xm);
  }
  
  return res;

  freetab(A);

}


// Pour les variables centrées aux mailles :   ind=0
//   delta i+1/2 = Z[i]-Z[i-1]
// Pour les variables frontières aux mailles :   ind=1
//   delta i = Z[i+1]-Z[i]
double delta(int so,int ind, double* Z, double** dk){
  double dz=0.0;
  int k=so-2;
  if (ind==0){
  	for(int l=0; l<=k; l++){
	  //dz+=dk[k][l]*(Z[l+1] - Z[l]);
	  dz+=dk[k][l]*(Z[l] - Z[-l-1]);
	}
  }
  else if (ind==1){
  	for(int l=0; l<=k; l++){
	  //dz+=dk[k][l]*(Z[l+1] - Z[l]);
	  dz+=dk[k][l]*(Z[l+1] - Z[-l]);
	}
  }
  return dz;
}


// Pour rho.eps  , retourne AV(PI.\delta u)
// Retourne la moyenne de YdZ  
// soit : . Y=PI et Z=u  pour EPS
//        . Y=u  et Z=PI pour EKIN
// int    cas==0 : EPS    cas==1 : EKIN
/*
double YdZold(int so, int cas, double* P, double* Q, double* u, double* X, double** dk, double** rk){
  double du, res=0.0;
  int o=so-2;
  int ind;
  double* YdZtemp, *YdZ;
  YdZtemp=malloc(4*so*sizeof(double));
  YdZ=&YdZtemp[2*so];
  if(cas==0){ // EPS+=(dt/dx).PIdu
  	ind=0; // variable sur maille primale
    for(int j=-o-2; j<=o+2; j++){
      // calcul de du
      du=delta(so,1, &u[j], dk);
      // PI.du
      YdZ[j]=(P[j]+Q[j])*du;
    }
  }
  else if(cas==1){ // EKIN+=(dt/dx).udPI
  	ind=1; // variable sur maille primale
    for(int j=-o-2; j<=o+2; j++){ // on augmente le nombre de points de calcul pour être sur d'avoir ceux dont on a besoin
      // calcul de dPI
      du=0.0;
      for(int l=0; l<=o; l++){
        du+=dk[o][l]*( (P[j+l]+Q[j+l])-(P[j-l-1]+Q[j-l-1]) );
      }
      // u.dPI
      YdZ[j]=u[j]*du;
    }
  }
  res=phibar(so,ind, YdZ, X, rk);
  return res;
  
  free(YdZtemp);
}
*/

// Pour rho.eps  , retourne AV(PI.\delta u)
// Retourne la moyenne de YdZ = Y.Z' (intégrale de Y . la dérive de Z /t à x) 
// soit : . Y=PI et Z=u  pour EPS
//        . Y=u  et Z=PI pour EKIN
// int    cas==0 : EPS    cas==1 : EKIN
double YdZ(int so, int cas, double* P, double* Q, double* u, double* X, double** rk){
  double res=0.0;
  int o=so-2;
  int n=2*(so-2)+2; // taille du sytème pour le calcul de la dérivé spatial
  double* temp=malloc((n-1)*sizeof(double));
  double monome_x, x;
  double akk, aik;

  // La matrice et le second membre sont stockée sur le même tableau
  double** A=alloctabd(n,n+1);

  if(cas==0){ // EPS+=(dt/dx).PIdu
  	// Calcul de la dérivé spatial de u  (primal)

	// Contruction de la matrice de Cramer A
	// sur le carré de taille n dans A
	for(int i=0; i<=n-1; i++){
	  for(int j=0; j<=n-1; j++){
	    A[i][j]=1.0;
	    // A=x^?
	    for(int k=1; k<=n-1-j; k++){ // quand j=n-1, on entre pas dans la boucle A[i][n-1]=1
	      A[i][j]*=X[-(o)+i]; 
	    } 
      } 
	}

    // Construction du second menbre 
    // sur la (n+1)ième colonne de A
    for(int j=0; j<=n-1; j++){
	  A[j][n]=u[-(o)+j];
	}

	// Résolution du sytème linèaire
	// ref fonction d'internet
	/* algorithme de Gauss-Jordan */
	for(int k=0;k<n;k++){
	  akk=A[k][k];
	  for(int j=k;j<n+1;j++){
	    A[k][j]=A[k][j]/akk;
	  }
	  for(int i=0;i<n;i++){
	    aik=A[i][k]; 
	    if(i!=k){
	      for(int j=k;j<n+1;j++){
	        A[i][j]=A[i][j]-aik*A[k][j];
	      }
	    }
	  }
	}

	// Construction des points (P+Q).dx(u) dans le tableau temp
	for(int i=0; i<=n-2; i++){
	  temp[i]=0;
      x=fxc(so,&X[-(o)+i],rk);
      monome_x=1.;
      for(int j=0; j<=n-2; j++){
      	temp[i]+=(j+1)*monome_x*A[n-2-j][n];
      	monome_x*=x;
      }
	  temp[i]*=P[-(o)+i]+Q[-(o)+i];  
	}

	// Calcul de la valeur moyenne
	res=phibar(so, 0, &temp[o], X, rk);
  	
  } ////////////////////////////////////////////////////////////////////////
  else if(cas==1){ // EKIN+=(dt/dx).udPI
  	// Calcul de la dérivé spatial de PI  (primal)
  	double *Xc=malloc(n*sizeof(double));
  	for(int i=0; i<=n-1; i++){
  	  Xc[i]=fxc(so,&X[-(o+1)+i], rk);
  	}

	// Contruction de la matrice de Cramer A
	// sur le carré de taille n dans A
	for(int i=0; i<=n-1; i++){
	  for(int j=0; j<=n-1; j++){
	    A[i][j]=1.0;
	    // A=x^?
	    for(int k=1; k<=n-1-j; k++){ // quand j=n-1, on entre pas dans la boucle A[i][n-1]=1
	      A[i][j]*=Xc[i]; 
	    } 
      }
	}
	
    // Construction du second menbre 
    // sur la (n+1)ième colonne de A
    for(int j=0; j<=n-1; j++){
	  A[j][n]=P[-(o+1)+j]+Q[-(o+1)+j];
	}

	// Résolution du sytème linèaire
	// ref fonction d'internet
	/* algorithme de Gauss-Jordan */
	for(int k=0;k<n;k++){
	  akk=A[k][k];
	  for(int j=k;j<n+1;j++){
	    A[k][j]=A[k][j]/akk;
	  }
	  for(int i=0;i<n;i++){
	    aik=A[i][k]; 
	    if(i!=k){
	      for(int j=k;j<n+1;j++){
	        A[i][j]=A[i][j]-aik*A[k][j];
	      }
	    }
	  }
	}

	// Construction des points u.dx(P+Q) dans le tableau temp
	for(int i=0; i<=n-2; i++){
	  temp[i]=0;
      x=X[-o+i];
      monome_x=1.;
      for(int j=0; j<=n-2; j++){
      	temp[i]+=(j+1)*monome_x*A[n-2-j][n];
      	monome_x*=x;
      }
	  temp[i]*=u[-o+i];  
	}

	// Calcul de la valeur moyenne
	res=phibar(so, 1, &temp[o], X, rk);
    
  }
  return res;

  free(temp);

}



// Fonction uniquement pour le kinetic energy fix  
// sum Qbar.Delta K
// k : spatial order
// Z = &Z[i]
// Valeur moyenne d'une variable duale à partir des ses valeurs ponctuelles primales.
double phiQK(int so, double* DK, double* X, double** rk){
  double res=0.0;
  int o=so-2;
  int n=2*(so-2)+2; // taille du sytème
  double akk, aik;

  // La matrice et le second membre sont stockée sur le même tableau
  double** A=alloctabd(n,n+1);

  // Contruction de la matrice de Cramer A
  // sur le carré de taille n dans A
  for(int i=0; i<=n-1; i++){
	for(int j=0; j<=n-1; j++){
	  A[i][j]=1.0;
	  // A=x^?
      for(int k=1; k<=n-1-j; k++){ // quand j=n-1, on entre pas dans la boucle A[i][n-1]=1
        A[i][j]*=X[(-o)+i]; 
      } 
    } 
  }


  // Construction du second menbre 
  // sur la (n+1)ième colonne de A
  for(int j=0; j<=n-1; j++){
    A[j][n]=DK[(-o)+j];
  }


  // Résolution du sytème linèaire
  // ref fonction d'internet
  /* algorithme de Gauss-Jordan */
  for(int k=0;k<n;k++){
	akk=A[k][k];
	for(int j=k;j<n+1;j++){
	  A[k][j]=A[k][j]/akk;
	}
	for(int i=0;i<n;i++){
	  aik=A[i][k]; 
	  if(i!=k){
	    for(int j=k;j<n+1;j++){
	      A[i][j]=A[i][j]-aik*A[k][j];
	    }
	  }
	}
  }


  // Calcul de la valeur moyenne
  res=0;
  double xp, xm;
  xm=X[0];
  xp=X[1];
  double monome_xm=xm, monome_xp=xp;
  // intégrale sur [xm, xp]
  for (int i=0; i<=n-1; i++){ 
    res+=(monome_xp-monome_xm)*A[n-1-i][n]/(i+1);
    monome_xm*=xm;
    monome_xp*=xp;
  }
  
  return res/(xp-xm);

  freetab(A);
}


// Met a jour les condition aux limites pour la variables Z
/* Si Z est centré aux maille  jfin= ifin
      Z est frontières aux mailles jfin=ifin+1

  int iu = 0 si cond sym
      iu = 1 si cond anti-sym
*/
void condlim(int tst, int jdeb,int jfin, int nbghosts, double* Z, int ind){
  if(tst==0 || tst==1 || tst==3){
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
  }
  // Onde Acoustique 
  else if(tst==2){
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
        Z[jdeb-1-j]=Z[jdeb+j];
        // à droite
        Z[jfin+1+j]=Z[jfin-j];
      }
    }
  }
}


// The Minmod limiter
double MinMod(double r){
  return fmax(0.,fmin(1.,r));
}

// The Superbee limiter
double Superbee(double r){
  return fmax(0., fmax(fmin(1.,2.*r),fmin(2.,r)) );
}

// The Christensen limiter
double Christensen(double r, double arg){
  return fmax(0., fmin(2.,fmin(2*r,arg)) );
}


double qvis(int q, double Cq, double Cl, double* U, double tau, double c, double* DM){
  int firstnumber=floor(q/10);
  int secondnumber=q-10*firstnumber;
  double du=U[1]-U[0];
  // ordre 1
  //   firstnumber==0
  // Ordre 2
  if (firstnumber==1){
    // delta vitesse
    double dUm=U[0]-U[-1];
    double dUp=U[2]-U[1];
    // delta masse
    double dM=DM[0];
    double dMm=DM[-1];
    double dMp=DM[1];

    double rplus=(dUp/dMp)*(dM/du);
    double rmoins=(dUm/dMm)*(dM/du);
    
    double phiplus=MinMod(rplus);
    double phimoins=MinMod(rmoins);
    
    /*
    double phiplus=Superbee(rplus);
    double phimoins=Superbee(rmoins);
    */

    /*
    double argplus=(dM*rplus+dMp)/(dM+dMp);
    double argmoins=(dM*rmoins+dMm)/(dM+dMm);
    double phiplus=Christensen(rplus,argplus);
    double phimoins=Christensen(rmoins, argmoins);
    */

    du*=(1-0.5*(phiplus + phimoins));
  }

  
  // von Neumann-Ritchmyer
  if (secondnumber == 0){
    return  -Cq*(du)*fabs(du)/tau;
  }
  // Rosenbluth
  else if(secondnumber == 1){
  	if (du<0){
  	  return -Cq*(du)*fabs(du)/tau;
  	}
  	else{ return 0;}
  }
  // Landschoff
  else if(secondnumber == 2){
  	if (du<0){
  	  return -Cq*(du)*fabs(du)/tau -Cl*c*du/tau;
  	}
  	else{ return 0;}
  }
  // "Magical" q combination
  else if(secondnumber == 3){
  	if (du<0){
  	  return -Cq*(du)*fabs(du)/tau -Cl*c*du/tau;
  	}
  	else{ return -Cl*c*du/tau;}
  }
  else if(secondnumber == 4){
    return 0.;
  }
  else{
    printf("Erreur dans le choix de pseudo q\n");
  }
}




void initCoef(double** Ck, double** Ckbar, double** dk, double** Qbar, double** rk){
  // Ck
  Ck[0][0]=1.0;               Ck[0][1]=0.0;               Ck[0][2]=0.0;             Ck[0][3]=0.0;           Ck[0][4]=0.0; 
  Ck[1][0]=13./12;            Ck[1][1]=-1./24;            Ck[1][2]=0.0;             Ck[1][3]=0.0;           Ck[1][4]=0.0; 
  Ck[2][0]=1067./960;         Ck[2][1]=-29./480;          Ck[2][2]=3./640.;         Ck[2][3]=0.0;           Ck[2][4]=0.0; 
  Ck[3][0]=30251./26880;      Ck[3][1]=-7621./107520;     Ck[3][2]=159./17920;      Ck[3][3]=-5./7168;      Ck[3][4]=0.0; 
  Ck[4][0]=5851067./5160960;  Ck[4][1]=-100027./1290240;  Ck[4][2]=31471./2580480;  Ck[4][3]=-425./258048;  Ck[4][4]=35./294912; 
  // Ckbar
  Ckbar[0][0]=1.0;                 Ckbar[0][1]=0.0;                Ckbar[0][2]=0.0;                 Ckbar[0][3]=0.0;              Ckbar[0][4]=0.0; 
  Ckbar[1][0]=11./12;              Ckbar[1][1]=1./24;              Ckbar[1][2]=0.0;                 Ckbar[1][3]=0.0;              Ckbar[1][4]=0.0; 
  Ckbar[2][0]=863./960;            Ckbar[2][1]=77./1440;           Ckbar[2][2]=-17./5760;           Ckbar[2][3]=0.0;              Ckbar[2][4]=0.0;
  Ckbar[3][0]=215641./241920;      Ckbar[3][1]=6361./107520;       Ckbar[3][2]=-281./53760;         Ckbar[3][3]=367./967680;      Ckbar[3][4]=0.0; 
  Ckbar[4][0]=41208059./46448640;  Ckbar[4][1]=3629953./58060800;  Ckbar[4][2]=-801973./116121600;  Ckbar[4][3]=49879./58060800;  Ckbar[4][4]=-27859./464486400;  
  // dk
  dk[0][0]=1.0;           dk[0][1]=0.0;         dk[0][2]=0.0;         dk[0][3]=0.0;           dk[0][4]=0.0; 
  dk[1][0]=9./8;          dk[1][1]=-1./24;      dk[1][2]=0.0;         dk[1][3]=0.0;           dk[1][4]=0.0; 
  dk[2][0]=75./64;        dk[2][1]=-25./384;    dk[2][2]=3./640;      dk[2][3]=0.0;           dk[2][4]=0.0; 
  dk[3][0]=1225./1024;    dk[3][1]=-245./3072;  dk[3][2]=49./5120;    dk[3][3]=-5./7168;      dk[3][4]=0.0; 
  dk[4][0]=19845./16384;  dk[4][1]=-735./8192;  dk[4][2]=567./40960;  dk[4][3]=-405./229376;  dk[4][4]=35./294912; 
  // Qbar
  Qbar[0][0]=1./2;              Qbar[0][1]=0.0;             Qbar[0][2]=0.0;           Qbar[0][3]=0.0;              Qbar[0][4]=0.0; 
  Qbar[1][0]=13./24;            Qbar[1][1]=-1./24;          Qbar[1][2]=0.0;           Qbar[1][3]=0.0;              Qbar[1][4]=0.0; 
  Qbar[2][0]=401./720;          Qbar[2][1]=-31./480;        Qbar[2][2]=11./1440;      Qbar[2][3]=0.0;              Qbar[2][4]=0.0; 
  Qbar[3][0]=68323./120960;     Qbar[3][1]=-353./4480;      Qbar[3][2]=1879./120960;  Qbar[3][3]=-191./120960;     Qbar[3][4]=0.0; 
  Qbar[4][0]=2067169./3628800;  Qbar[4][1]=-40111./453600;  Qbar[4][2]=581./25920;    Qbar[4][3]=-28939./7257600;  Qbar[4][4]=2497./7257600; 
  // rk
  rk[0][0]=1./2;          rk[0][1]=0.0;           rk[0][2]=0.0;         rk[0][3]=0.0;          rk[0][4]=0.0; 
  rk[1][0]=9./16;         rk[1][1]=-1./16;        rk[1][2]=0.0;         rk[1][3]=0.0;          rk[1][4]=0.0; 
  rk[2][0]=75./128;       rk[2][1]=-25./256;      rk[2][2]=3./256.;     rk[2][3]=0.0;          rk[2][4]=0.0; 
  rk[3][0]=1225./2048;    rk[3][1]=-245./2048;    rk[3][2]=49./2048;    rk[3][3]=-5./2048;     rk[3][4]=0.0; 
  rk[4][0]=19845./32768;  rk[4][1]=-2205./16384;  rk[4][2]=567./16384;  rk[4][3]=-405./65536;  rk[4][4]=35./65536; 
}


void u0(int tst, double a, double b,double x, double* W){ 
  double eps=1e-8; 
  double k=2*M_PI; //////////////////////////
  double rho0=1.0; double p0=5./7;
  if (tst==0){ // Sod
    if(x<=a+(b-a)/2.0){
      W[0]=1.0; // rho
      W[1]=0.0; // u
      W[2]=1.0; // p
    }
    else{
      W[0]=0.125; // rho
      W[1]=0.0; // u
      W[2]=0.1; // p
      //W[2]=1.0;
    }
  }
  else if (tst==1){ // Bizarrium
    if(x<=a+(b-a)/2.0){
      W[0]=1./0.7e-4; // rho
      W[1]=0.0;       // u
      W[2]=1e11;      // p
    }
    else{
      W[0]=1e4; // rho
      W[1]=250.0; // u
      W[2]=0.0; // p
    }
  }
  else if (tst==2){ // Onde acoustique
    W[0]=rho0 + eps*sin(k*x);  // rho
    W[1]=eps*sin(k*x);  // u
    W[2]=p0 + eps*sin(k*x);  // p
  }
  else if (tst==3){ // 3 etats
    if(x<=a+2.0*(b-a)/5){
      W[0]=0.125; // rho
      W[1]=0.0;       // u
      W[2]=0.1;      // p
    }
    else if(a+3.0*(b-a)/5<=x){
      W[0]=0.125; // rho
      W[1]=0.0;  // u
      W[2]=0.1; // p
    }
    else{  // if(2.0*(a+b)/5<x || x<3.0*(a+b)/5){
      W[0]=1.0; // rho
      W[1]=0.0; // u
      W[2]=1.0; // p
    }
    
  }
  else{ printf("Erreur dans le choix de tst\n");}
}

void initButcher(int sch, double** A, double* THETA, double* ALPHA){
  // RK1  (Matsuno forward-backward)
  if (sch==0){
    A[0][0]=1.0;

    THETA[0]=0.0; THETA[1]=1.0;

    ALPHA[0]=1.;
  }
  // RK2 (SSP) HEUN
  else if (sch==1){
    A[0][0]=1.0;

    THETA[0]=1.0/2.0; THETA[1]=1.0/2.0;

    ALPHA[0]=1.;
  }
  // RK3 (SSP)
  else if (sch==2){
    A[0][0]=1.0; 
    A[1][0]=1.0/4.0; A[1][1]=1.0/4.0;

    THETA[0]=1.0/6.0; THETA[1]=1.0/6.0; THETA[2]=2.0/3.0;

    ALPHA[0]=1.; ALPHA[1]=1./2;
  }
  // KUTTA ORDRE 4 (Classic RK4)
  else if (sch==3){
    A[0][0]=1.0/2.0; 
    A[1][0]=0.0;     A[1][1]=1.0/2.0;
    A[2][0]=0.0;     A[2][1]=0.0;      A[2][2]=1.0;

    THETA[0]=1.0/6.0; THETA[1]=1.0/3.0; THETA[2]=1.0/3.0; THETA[3]=1.0/6.0;

    ALPHA[0]=1./2; ALPHA[1]=1./2; ALPHA[2]=1.;
  }
  // Cash-Karp (ordre 5)
  else if (sch==4){
    A[0][0]=1.0/5.0; 
    A[1][0]=3./40;        A[1][1]=9./40;
    A[2][0]=3./10;        A[2][1]=-9./10;    A[2][2]=6./5;
    A[3][0]=-11./54;      A[3][1]=5./2;      A[3][2]=-70./27;     A[3][3]=-35./27;
    A[4][0]=1631./55296;  A[4][1]=175./512;  A[4][2]=575./13824;  A[4][3]=44275./110592;  A[4][4]=253./4096;

    THETA[0]=37./378; THETA[1]=0.; THETA[2]=250./621; THETA[3]=125./594; THETA[4]=0.; THETA[5]=512./1771;  

    ALPHA[0]=1./5; ALPHA[1]=3./10; ALPHA[2]=3./5; ALPHA[3]=1.; ALPHA[4]=7./8;
  }
  // Dormand-Prince
  else if (sch==5){
    A[0][0]=1.0/5.0; 
    A[1][0]=3./40;         A[1][1]=9./40;
    A[2][0]=44./45;        A[2][1]=-56./15;        A[2][2]=32./9;
    A[3][0]=19372./6561;   A[3][1]=-25360./2187;   A[3][2]=64448./6561;   A[3][3]=-212./729;
    A[4][0]=9017./3168;    A[4][1]=-355./33;       A[4][2]=46732./5247;  A[4][3]=49./176;    A[4][4]=-5103./18656;
    A[5][0]=35./384;       A[5][1]=0.;             A[5][2]=500./1113;    A[5][3]=125./192;   A[5][4]=-2187./6784;  A[5][5]=11./84;

    THETA[0]=35./384; THETA[1]=0.; THETA[2]=500./1113; THETA[3]=125./192; THETA[4]=-2187./6784; THETA[5]=11./84; THETA[6]=0.;

    ALPHA[0]=1./5; ALPHA[1]=3./10; ALPHA[2]=4./5; ALPHA[3]=8./9; ALPHA[4]=1.; ALPHA[5]=1.;
  }
  // HEUN  ordre 3
  else if (sch==6){
    A[0][0]=1.0/3; 
    A[1][0]=0; A[1][1]=2.0/3;

    THETA[0]=1.0/4; THETA[1]=0; THETA[2]=3./4;

    ALPHA[0]=1./3; ALPHA[1]=2./3;
  }
  // RALSTON  ordre 3
  else if (sch==7){
    A[0][0]=1.0/2; 
    A[1][0]=0; A[1][1]=3.0/4;

    THETA[0]=2.0/9; THETA[1]=1./3; THETA[2]=4./9;

    ALPHA[0]=1./2; ALPHA[1]=3./4;
  }
  // Bogacki-Shampine ordre 3 
  else if (sch==8){
    A[0][0]=1.0/2.0; 
    A[1][0]=0.0;     A[1][1]=3.0/4;
    A[2][0]=2./9;    A[2][1]=1./3;      A[2][2]=4./9;

    THETA[0]=2./9; THETA[1]=1.0/3; THETA[2]=4./9; THETA[3]=0;

    ALPHA[0]=1./2; ALPHA[1]=3./4; ALPHA[2]=1.; 
  }
  // Bogacki-Shampine ordre 2 
  else if (sch==9){
    A[0][0]=1.0/2.0; 
    A[1][0]=0.0;     A[1][1]=3.0/4;
    A[2][0]=2./9;    A[2][1]=1./3;      A[2][2]=4./9;

    THETA[0]=7./24; THETA[1]=1.0/4; THETA[2]=1./3; THETA[3]=1./8;

    ALPHA[0]=1./2; ALPHA[1]=3./4; ALPHA[2]=1.; 
  }
  // KUTTA odre 3
  else if (sch==10){
    A[0][0]=1./2; 
    A[1][0]=-1.; A[1][1]=2.0;

    THETA[0]=1./6; THETA[1]=2./3; THETA[2]=1./6;

    ALPHA[0]=1./2; ALPHA[1]=1.;
  }
  // KUTTA ordre 2 (Explicit Midpoint Method) 
  else if (sch==11){
    A[0][0]=1.0/2;

    THETA[0]=0.0; THETA[1]=1.0;

    ALPHA[0]=1./2;
  }
  // KUTTA ordre 2 Ralston
  else if (sch==12){
    A[0][0]=2.0/3;

    THETA[0]=1.0/4; THETA[1]=3./4;

    ALPHA[0]=2./3;
  }
  // EULER forward
  else if (sch==13){
    A[0][0]=0;

    THETA[0]=0; THETA[1]=1.0;

    ALPHA[0]=0.0;
  }
  // SSPRK4(5,4) Spiteri-Ruuth
  else if (sch==14){
    A[0][0]=0.39175222700392; 
    A[1][0]=0.21766909633821;  A[1][1]=0.36841059262959;
    A[2][0]=0.08269208670950;  A[2][1]=0.13995850206999;  A[2][2]=0.25189177424738;
    A[3][0]=0.06796628370320;  A[3][1]=0.11503469844438;  A[3][2]=0.20703489864929;  A[3][3]=0.54497475021237;

    THETA[0]=0.14681187618661; THETA[1]=0.24848290924556; THETA[2]=0.10425883036650; THETA[3]=0.27443890091960; THETA[4]=0.22600748319395; 

    ALPHA[0]=0.39175222700392; ALPHA[1]=0.58607968896779; ALPHA[2]=0.47454236302687; ALPHA[3]=0.93501063100924; 
  }
  // SSPRK3(4,3) Spiteri-Ruuth
  else if (sch==15){
    A[0][0]=0.5; 
    A[1][0]=0.5;   A[1][1]=0.5;
    A[2][0]=1./6;  A[2][1]=1./6;  A[2][2]=1./6;

    THETA[0]=1./6; THETA[1]=1./6; THETA[2]=1./6; THETA[3]=0.5; 

    ALPHA[0]=0.5; ALPHA[1]=1.0; ALPHA[2]=0.5; 
  }
  // SSPRK3(5,3) Spiteri-Ruuth
  else if (sch==16){
    A[0][0]=0.37726891511710; 
    A[1][0]=0.75453783023419;  A[1][1]=0.37726891511710;
    A[2][0]=0.49056882269314;  A[2][1]=0.16352294089771;  A[2][2]=0.16352294089771;
    A[3][0]=0.78784303014311;  A[3][1]=0.14831273384724;  A[3][2]=0.14831273384724;  A[3][3]=0.34217696850008;

    THETA[0]=0.19707596384481; THETA[1]=0.11780316509765; THETA[2]=0.11709725193772; THETA[3]=0.27015874934251; THETA[4]=0.29786487010104; 

    ALPHA[0]=0.37726891511710; ALPHA[1]=0.75453783023419; ALPHA[2]=0.49056882269314; ALPHA[3]=0.78784303014311; 
  }
  // SPPRK1(3,1) Spiteri-Ruuth
  else if (sch==17){
    A[0][0]=1./3; 
    A[1][0]=1./3; A[1][1]=1./3;

    THETA[0]=1./3; THETA[1]=1./3; THETA[2]=1./3;

    ALPHA[0]=1./3; ALPHA[1]=2./3;
  }
  // SPPRK2(3,2) Spiteri-Ruuth
  else if (sch==18){
    A[0][0]=1./2; 
    A[1][0]=1./2; A[1][1]=1./2;

    THETA[0]=1./3; THETA[1]=1./3; THETA[2]=1./3;

    ALPHA[0]=1./2; ALPHA[1]=1.;
  }
  // SPPRK2(4,2) Spiteri-Ruuth
  else if (sch==19){
    A[0][0]=1./3; 
    A[1][0]=1./3; A[1][1]=1./3;
    A[2][0]=1./3; A[2][1]=1./3;  A[2][2]=1./3;

    THETA[0]=1./4; THETA[1]=1./4; THETA[2]=1./4; THETA[3]=1./4;

    ALPHA[0]=1./3; ALPHA[1]=2./3; ALPHA[2]=1.;
  }
  // SPPRK1(2,1) Spiteri-Ruuth
  else if (sch==20){
    A[0][0]=1./2; 

    THETA[0]=1./2; THETA[1]=1./2; 

    ALPHA[0]=1./2; 
  }
  // SPPRK3(4,3) Spiteri-Ruuth
  else if (sch==21){
    A[0][0]=1./2; 
    A[1][0]=1./2; A[1][1]=1./2;
    A[2][0]=1./6; A[2][1]=1./6;  A[2][2]=1./6;

    THETA[0]=1./6; THETA[1]=1./6; THETA[2]=1./6; THETA[3]=1./2;

    ALPHA[0]=1./2; ALPHA[1]=1.; ALPHA[2]=1./2;
  }
  else{ printf("Erreur choix sch (initButcher)\n");}
}


