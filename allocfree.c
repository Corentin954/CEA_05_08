#include <stdlib.h>
#include "head.h"

/* 
--------------------------------------------------------------------------------
  Cette fonction libere la memoire allouee par alloctab.
--------------------------------------------------------------------------------
*/
void freetab(void *ptr) {
  void **ptrT=ptr;
  free(ptrT[0]);
  free(ptr);
}

// Version avec des int
int **alloctab(int dim1, int dim2) {
  int **ptr;
  ptr = malloc(dim1*sizeof(int *));
  if (ptr != NULL) {
    int i, taille_ligne = dim2*sizeof(int);
    int *tmp = malloc(dim1*taille_ligne);
    if (tmp != NULL) {
      for (i=0; i<dim1; i++) {
  	    ptr[i] = tmp;
  	    tmp += dim2;
  	  }
    }
    else
      ptr = NULL;
    }
  return(ptr);
}

// Version avec des double
double **alloctabd(int dim1, int dim2) {
  double **ptr;
  ptr = malloc(dim1*sizeof(double *));
  if (ptr != NULL) {
    int i, taille_ligne = dim2*sizeof(double);
    double *tmp = malloc(dim1*taille_ligne);
    if (tmp != NULL) {
      for (i=0; i<dim1; i++) {
  	    ptr[i] = tmp;
  	    tmp += dim2;
  	  }
    } else
      ptr = NULL;
    }
  return(ptr);
}
