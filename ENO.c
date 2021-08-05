#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "head2.h"


//  Fonction necessaire pour la procedure ENO



// Fonction stencil
/*  renvoie le stencil qui ne comporte pas de discocntuité
    l'alloctation de du tableau S se fait en dehors
    double* S=malloc(so*sizeof(double))
*/
void stencil(int so, double* X, double* Z, double** S){
  int ind_left=0;
  int nb=2;  // taille du stencil
  double resL, resR;
  for(int i=3; i<=so; i++){
  	resL=diff_div(nb,ind_left-1, X, Z);
  	resR=diff_div(nb,ind_left, X, Z);
    if(resL<resR){
      ind_left--;
    }
    nb++;
  }
  
  // valeur de sortie
  for(int i=0; i<=so-1 ; i++){
  	S[0][i]=X[ind_left+i];
  	S[1][i]=Z[ind_left+i];
  }

}


// Fonction  différence divisées
/*  renvoie la différence divisé en [X[1],X[0]] 
    à l'ordre nb
*/
double diff_div(int nb, int left, double* X, double* Z){
  double** tab=alloctabd(nb,nb);

  for(int i=0; i<=nb-1 ; i++){
  	tab[0][i]=Z[left+i];
  }

  // Routine de différence divisées
  for(int i=1; i<=nb-1; i++){
  	for(int j=i; j<=nb-1; j++){
  	  tab[i][j]=(tab[i-1][j]-tab[i-1][j-1])/(X[j]-X[j-i]);
  	}
  }
  
  return tab[nb-1][nb-1];

  freetab(tab);
}




// Z :  variable AVERAGE dont on veut la valeur ponctuelle (&Z[i])
// X :  maillage primal   (&X[i])
// ind : 0-centré  1-décentré
double phipw(int so, int ind, double* Z, double* X, double** rk){
  // entier utile :
  int o=so-2; //ordre de la méthode
  int n=2*(so-2)+2; // taille du sytème
  double akk, aik;
  double res=0.0;

  // Calcul du stencil de calcul
  double** S=alloctabd(2,n);
  if (ind==0){
  	// Xc=fxc (X) ...
    stencil(so, Xc, Z, S);
  }
  else if(tst==1){
    stencil(so, X, Z, S);	
  }

  // IDEE  :  renvoyer ind_left et remplacer le -(o) par -(ind_left) dans le code dessous
  //            plutot que de renvoyer la matrice S
  // X[-(o+1)+i]

  

  
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
