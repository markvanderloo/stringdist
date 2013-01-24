
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include "utils.h"

/* Optimal string alignment algorithm. 
 * Computes Damerau-Levenshtein distance, restricted to single transpositions.
 * - See pseudocode at http://en.wikipedia.org/wiki/Damerau%E2%80%93Levenshtein_distance
 * - Extended with custom weights and maxDistance
 */
static double osa(const char *a, int na, const char *b, int nb, double *weight, double maxDistance, double *scores){
   int i, j;
   int I = na+1, J = nb+1;
   double sub, tran;

   for ( i = 0; i < I; ++i ){
      scores[i] = i;
   }
   for ( j = 1; j < J; ++j ){
      scores[I*j] = j;
   }

   for ( i = 1; i <= na; ++i ){
      for ( j = 1; j <= nb; ++j ){
         if (a[i-1] == b[j-1]){
            sub = 0;
            tran= 0;
         } else {
            sub = weight[2];
            tran= weight[3];
         }
         
         scores[i + I*j] = min3( 
            scores[i-1 + I*j    ] + weight[0],     // deletion
            scores[i   + I*(j-1)] + weight[1],     // insertion
            scores[i-1 + I*(j-1)] + sub            // substitution
         );
         if ( i>1 && j>1 && a[i-1] == b[j-2] && a[i-2] == b[j-1] ){
            scores[i + I*j] = min2(scores[i + I*j], scores[i-2+I*(j-2)]) + tran; // transposition
         }
         if ( maxDistance > 0 && scores[i + I*j] > maxDistance ){
            return -1;
         }
      }
   }
   return(scores[I*J-1]);
}

//-- interface with R

static int vmax(int *x, int n){
   double m = x[0];
   for ( int i = 1; i < n; ++i ){
      if ( x[i] > m ){ 
         m = x[i];
      }
   }
   return(m);
}

SEXP R_osa(SEXP a, SEXP b, SEXP ncharA, SEXP ncharB, SEXP weight, SEXP maxDistance){
   PROTECT(a);
   PROTECT(b);
   PROTECT(ncharA);
   PROTECT(ncharB);
   PROTECT(weight);
   PROTECT(maxDistance);

   int dsize, na = length(a), nb = length(b);
   double *scores; 
   double *w = REAL(weight);
   double maxDist = REAL(maxDistance)[0];

   // determine workspace size 
   dsize = (vmax(INTEGER(ncharA),na)+1) * (vmax(INTEGER(ncharB),nb)+1);
   scores = calloc(dsize, sizeof(double)); 
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
      if (STRING_ELT(a,i) == NA_STRING || STRING_ELT(b,j) == NA_STRING){
         y[k] = NA_REAL;
         continue;
      }
      y[k] = osa(
         CHAR(STRING_ELT(a,i)), 
         INTEGER(ncharA)[i], 
         CHAR(STRING_ELT(b,j)), 
         INTEGER(ncharB)[j], 
         w,
         maxDist,
         scores
      );
   }
   
   free(scores);
   UNPROTECT(7);
   return(yy);
}



