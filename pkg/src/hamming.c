/* Hamming distance function
 */

#define USE_RINTERNALS
#include <R.h>
#include <Rdefines.h>
#include "utils.h"

double hamming(unsigned int *a, unsigned int *b, int n, int maxDistance){
   double h=0.0;
   for(int i=0; i<n; ++i){
      if (a[i] != b[i]) h++;
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
  PROTECT(yy = allocVector(REALSXP,nt));
  double *y = REAL(yy);
  int i=0, j=0, k=0, nchar;
  int maxDist = INTEGER(maxDistance)[0];

  for ( k=0; k<nt; ++k){
    if ( INTEGER(VECTOR_ELT(a,i))[0] == NA_INTEGER || INTEGER(VECTOR_ELT(b,j))[0] == NA_INTEGER ){
      y[k] = NA_REAL;
      continue;         
    }
    nchar = length(VECTOR_ELT(a,i));
    if ( nchar != length(VECTOR_ELT(b,j)) ){
      y[k] = R_PosInf;
      continue;
    }
    y[k] = hamming(
      (unsigned int *) INTEGER(VECTOR_ELT(a,i)),
      (unsigned int *) INTEGER(VECTOR_ELT(b,j)),
      nchar,
      maxDist
    );
    if (y[k] < 0) y[k] = R_PosInf;
    i = RECYCLE(i+1,na);
    j = RECYCLE(j+1,nb);
  }
  UNPROTECT(4);
  return yy;
}



