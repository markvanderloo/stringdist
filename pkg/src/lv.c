
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include "utils.h"

/* Levenshtein distance
 * Computes Levenshtein distance
 * - Simplified from restricted DL pseudocode at http://en.wikipedia.org/wiki/Damerau%E2%80%93Levenshtein_distance
 * - Extended with custom weights and maxDistance
 */
static double osa(unsigned int *a, int na, unsigned int *b, int nb, double *weight, double maxDistance, double *scores){
   if (na == 0) return(nb);
   if (nb == 0) return(na);
   int i, j;
   int I = na+1, J = nb+1;
   double sub;

   for ( i = 0; i < I; ++i ){
      scores[i] = i;
   }
   for ( j = 1; j < J; ++j ){
      scores[I*j] = j;
   }

   for ( i = 1; i <= na; ++i ){
      for ( j = 1; j <= nb; ++j ){
         sub = (a[i-1] == b[j-1]) ? 0 : weight[2];
         scores[i + I*j] = min3( 
            scores[i-1 + I*j    ] + weight[0],     // deletion
            scores[i   + I*(j-1)] + weight[1],     // insertion
            scores[i-1 + I*(j-1)] + sub            // substitution
         );
         if ( maxDistance > 0 && scores[i + I*j] > maxDistance ){
            return -1;
         }
      }
   }
   return(scores[I*J-1]);
}

//-- interface with R


SEXP R_lv(SEXP a, SEXP b, SEXP weight, SEXP maxDistance){
   PROTECT(a);
   PROTECT(b);
   PROTECT(weight);
   PROTECT(maxDistance);

   int na = length(a), nb = length(b);
   double *scores; 
   double *w = REAL(weight);
   double maxDist = REAL(maxDistance)[0];

   scores = calloc(get_dp_matrix_size(a,b), sizeof(double)); 
   if ( scores == NULL ){
      error("%s\n","unable to allocate enough memory for workspace");
   }

   // output vector
   int nt = (na > nb) ? na : nb;   
   int i,j,k;
   SEXP yy;
   PROTECT(yy = allocVector(REALSXP, nt));
   double *y = REAL(yy);   
   
   for ( int k=0; k < nt; ++k ){
      i = k % na;
      j = k % nb;      
      if (INTEGER(VECTOR_ELT(a,i))[0] == NA_INTEGER || INTEGER(VECTOR_ELT(b,j))[0] == NA_INTEGER){
         y[k] = NA_REAL;
         continue;
      }
      y[k] = osa(
         INTEGER(VECTOR_ELT(a,i)), 
         length(VECTOR_ELT(a,i)), 
         INTEGER(VECTOR_ELT(b,j)), 
         length(VECTOR_ELT(b,j)), 
         w,
         maxDist,
         scores
      );
   }
   
   free(scores);
   UNPROTECT(5);
   return(yy);
}



