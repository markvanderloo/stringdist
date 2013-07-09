
#define USE_RINTERNALS
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include "utils.h"

/* Longest common substring
 * - basically edit distance, only allowing insertions and substitutions, at the cost of 1.
 */
static double lcs(unsigned int *a, int na, unsigned int *b, int nb, int maxDistance, int *scores){
  if (na == 0){
    if ( maxDistance > 0 && maxDistance < nb ){
      return -1.0;
    } else {
      return (double) nb;
    }
  }
  if (nb == 0){
    if (maxDistance > 0 && maxDistance < na){
      return -1.0;
    } else {
      return (double) na;
    }
  }

  int i, j;
  int I = na+1, J = nb+1;
  int mincol;
  for ( i = 0; i < I; ++i ){
    scores[i] = i;
  }
  for ( j = 1; j < J; ++j ){
    scores[I*j] = j;
  }

  for ( i = 1; i <= na; ++i ){
    mincol = na + nb + 1;
    for ( j = 1; j <= nb; ++j ){
      if ( a[i-1] == b[j-1] ){ // equality, copy previous score
        scores[i + I*j] = scores[i-1 + I*(j-1)];
      } else {
        scores[i + I*j] = min2(
          scores[i-1 + I*j    ] + 1 ,     // deletion
          scores[i   + I*(j-1)] + 1       // insertion
        );
      }
      mincol = min2(mincol,scores[i+I*j]);
    }
    if ( maxDistance > 0 && mincol > maxDistance ){
      return -1;
    }
  }
  double dist = (double) scores[I*J-1];
  return dist;
}

//-- interface with R


SEXP R_lcs(SEXP a, SEXP b, SEXP maxDistance){
   PROTECT(a);
   PROTECT(b);
   PROTECT(maxDistance);

   int na = length(a), nb = length(b);
   int *scores; 
   int maxDist = INTEGER(maxDistance)[0];

   scores = (int *) malloc( (max_length(a) + 1) * (max_length(b) + 1) * sizeof(int)); 
   if ( scores == NULL ){
      error("%s\n","unable to allocate enough memory for workspace");
   }

   // output vector
   int nt = (na > nb) ? na : nb;   
   int i=0, j=0;
   SEXP yy;
   PROTECT(yy = allocVector(REALSXP, nt));
   double *y = REAL(yy);   
   
   for ( int k=0; k < nt; ++k ){
      if (INTEGER(VECTOR_ELT(a,i))[0] == NA_INTEGER || INTEGER(VECTOR_ELT(b,j))[0] == NA_INTEGER){
         y[k] = NA_REAL;
         continue;
      }
      y[k] = lcs(
        (unsigned int *) INTEGER(VECTOR_ELT(a,i)), 
         length(VECTOR_ELT(a,i)), 
        (unsigned int *) INTEGER(VECTOR_ELT(b,j)), 
         length(VECTOR_ELT(b,j)), 
         maxDist,
         scores
      );
      if (y[k] < 0 ) y[k] = R_PosInf;
      i = RECYCLE(i+1,na);
      j = RECYCLE(j+1,nb);
   }
   
   free(scores);
   UNPROTECT(4);
   return(yy);
}

