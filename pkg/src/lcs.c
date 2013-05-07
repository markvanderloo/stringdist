
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include "utils.h"

/* Longest common substring
 * - basically edit distance, only allowing insertions and substitutions, at the cost of 1.
 */
static int lcs(unsigned int *a, int na, unsigned int *b, int nb, int maxDistance, int *scores){
   int i, j;
   int I = na+1, J = nb+1;

   for ( i = 0; i < I; ++i ){
      scores[i] = i;
   }
   for ( j = 1; j < J; ++j ){
      scores[I*j] = j;
   }

   for ( i = 1; i <= na; ++i ){
      for ( j = 1; j <= nb; ++j ){

         if ( a[i-1] == b[j-1] ){ // equality, copy previous score
            scores[i + I*j] = scores[i-1 + I*(j-1)];
         } else {
          scores[i + I*j] = min2(
            scores[i-1 + I*j    ] + 1 ,     // deletion
            scores[i   + I*(j-1)] + 1       // insertion
          );
         }
         if ( maxDistance > 0 && scores[i + I*j] > maxDistance ){
            return -1;
         }
      }
   }
   return(scores[I*J-1]);
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
   int i,j,k;
   SEXP yy;
   PROTECT(yy = allocVector(INTSXP, nt));
   int *y = INTEGER(yy);   
   
   for ( int k=0; k < nt; ++k ){
      i = k % na;
      j = k % nb;
      if (INTEGER(VECTOR_ELT(a,i))[0] == NA_INTEGER || INTEGER(VECTOR_ELT(b,j))[0] == NA_INTEGER){
         y[k] = NA_REAL;
         continue;
      }
      y[k] = lcs(
         INTEGER(VECTOR_ELT(a,i)), 
         length(VECTOR_ELT(a,i)), 
         INTEGER(VECTOR_ELT(b,j)), 
         length(VECTOR_ELT(b,j)), 
         maxDist,
         scores
      );
   }
   
   free(scores);
   UNPROTECT(4);
   return(yy);
}

