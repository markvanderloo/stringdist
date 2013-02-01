/* Hamming distance function
 */

#include <R.h>
#include <Rdefines.h>

int hamming(int *a, int *b, int n, int maxDistance){
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
   int i,j,k,nchar;
   int maxDist = INTEGER(maxDistance)[0];

   for ( k=0; k<nt; ++k){
      i = k % na;
      j = k % nb;
      if ( INTEGER(VECTOR_ELT(a,i))[0] == NA_INTEGER || INTEGER(VECTOR_ELT(b,j))[0] == NA_INTEGER ){
         y[k] = NA_INTEGER;
         continue;         
      }
      nchar = length(VECTOR_ELT(a,i));
      if ( nchar != length(VECTOR_ELT(b,j)) ){
         error("Characters strings a[%d] and b[%d] have different number of characters", i, j);
      }
      y[k] = hamming(
         INTEGER(VECTOR_ELT(a,i)),
         INTEGER(VECTOR_ELT(b,j)),
         nchar,
         maxDist
      );
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
