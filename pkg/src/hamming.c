/* Hamming distance function
 */

#define USE_RINTERNALS
#include <R.h>
#include <Rdefines.h>
#include "utils.h"

int hamming(unsigned int *a, unsigned int *b, int n, int maxDistance){
   int i, h=0;
   for(i=0; i<n; ++i){
      h += (a[i] == b[i]) ? 0 : 1;
      if ( maxDistance > 0 && maxDistance < h ){
         return -1;
      }
   }
   return h;
}


// -- R interface

SEXP R_hm(SEXP a, SEXP b, SEXP maxDistance){
  PROTECT(a);
  PROTECT(b);
  PROTECT(maxDistance);

  int na = length(a);
  int nb = length(b);
  int nt = ( na > nb) ? na : nb;
  SEXP yy;
  PROTECT(yy = allocVector(INTSXP,nt));
  int *y = INTEGER(yy);
  int i=0, j=0, k=0, nchar;
  int maxDist = INTEGER(maxDistance)[0];

  for ( k=0; k<nt; ++k){
    if ( INTEGER(VECTOR_ELT(a,i))[0] == NA_INTEGER || INTEGER(VECTOR_ELT(b,j))[0] == NA_INTEGER ){
      y[k] = NA_INTEGER;
      continue;         
    }
    nchar = length(VECTOR_ELT(a,i));
    if ( nchar != length(VECTOR_ELT(b,j)) ){
      y[k] = -1;
      continue;
    }
    y[k] = hamming(
      (unsigned int *) INTEGER(VECTOR_ELT(a,i)),
      (unsigned int *) INTEGER(VECTOR_ELT(b,j)),
      nchar,
      maxDist
    );
    i = RECYCLE(i+1,na);
    j = RECYCLE(j+1,nb);
  }
  UNPROTECT(4);
  return yy;
}
/*
#include <stdio.h>
void main(){
   int h = hamming("aa","ab",0);
   printf("h = %d\n",h);
}

*/
