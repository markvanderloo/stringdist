
#define USE_RINTERNALS
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include "utils.h"

/* Longest common substring
 * - basically edit distance, only allowing insertions and substitutions, at the cost of 1.
 */
static int lcs(unsigned int *a, int na, unsigned int *b, int nb, int maxDistance, int *scores){
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
  return scores[I*J - 1];
}

//-- interface with R


SEXP R_lcs(SEXP a, SEXP b, SEXP maxDistance){
  PROTECT(a);
  PROTECT(b);
  PROTECT(maxDistance);

  int na = length(a), nb = length(b);
  int *scores; 
  int maxDist = INTEGER(maxDistance)[0];

  // space for the workfunction
  scores = (int *) malloc( (max_length(a) + 1) * (max_length(b) + 1) * sizeof(int)); 
  if ( scores == NULL ){
    UNPROTECT(3);
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



//-- Match function interface with R

SEXP R_match_lcs(SEXP x, SEXP table, SEXP nomatch, SEXP matchNA, SEXP maxDistance){
  PROTECT(x);
  PROTECT(table);
  PROTECT(nomatch);
  PROTECT(matchNA);
  PROTECT(maxDistance);

  int nx = length(x), ntable = length(table);
  int no_match = INTEGER(nomatch)[0];
  int match_na = INTEGER(matchNA)[0];
  int max_dist = INTEGER(maxDistance)[0];

  // space for the workfunction
  int *scores = (int *) malloc( (max_length(x) + 1) * (max_length(table) + 1) * sizeof(int)); 
  if ( scores == NULL ){
    UNPROTECT(3);
    error("%s\n","unable to allocate enough memory for workspace");
  }

  // output vector
  SEXP yy;
  PROTECT(yy = allocVector(INTSXP, nx));
  int *y = INTEGER(yy);

  int *X, *T;
  double d = R_PosInf, d1 = R_PosInf;
  int index, xNA, tNA;

  for ( int i=0; i<nx; i++){
    index = no_match;

    X = INTEGER(VECTOR_ELT(x,i));
    xNA = (X[0] == NA_INTEGER);

    for ( int j=0; j<ntable; j++){

      T = INTEGER(VECTOR_ELT(table,j));
      tNA = (T[0] == NA_INTEGER);

      if ( !xNA && !tNA ){        // both are char (usual case)
        d = (double) lcs(
          (unsigned int *) X, 
          length(VECTOR_ELT(x,i)), 
          (unsigned int *) T,
          length(VECTOR_ELT(table,j)), 
          max_dist,
          scores
        );
        if ( d > -1 && d < d1){ 
          index = j + 1;
          if ( d == 0.0  ) break;
          d1 = d;
        }
      } else if ( xNA && tNA ) {  // both are NA
        index = match_na ? j + 1 : no_match;
        break;
      }
    }
    
    y[i] = index;
  }
  free(scores);  
  UNPROTECT(6);
  return(yy);
}
