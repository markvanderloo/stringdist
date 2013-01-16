
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>

static double min3(double x, double y, double z){
   if ( x <= y && x <= z ){ 
      return(x);
   } else if ( y <= x && y <= z ) {
      return(y);
   } else {
      return(z);
   }
}

static double min2(double x, double y){
   if ( x <= y ){
      return(x);
   } else {
      return(y);
   }
}

/* Optimal string alignment algorithm. 
 * Computes Damerau-Levenshtein distance, restricted to single transpositions.
 * Algorithm taken from http://en.wikipedia.org/wiki/Damerau%E2%80%93Levenshtein_distance
 * and extended with custom weights.
 */
int osa(const char *a, int na, const char *b, int nb, double *w, double *d){
   int i, j, sub=0, tran=0;
   int I = na+1, J = nb+1;

   for ( i = 0; i < I; ++i ){
      d[i] = i;
   }
   for ( j = 1; j < J; ++j ){
      d[I*j] = j;
   }

   for ( i = 1; i <= na; ++i ){
      for ( j = 1; j <= nb; ++j ){
         if (a[i-1] == b[j-1]){
            sub = 0;
            tran= 0;
         } else {
            sub = w[2];
            tran= w[3];
         }
         
         d[i + I*j] = min3( 
            d[i-1 + I*j    ] + w[0],     // deletion
            d[i   + I*(j-1)] + w[1],     // insertion
            d[i-1 + I*(j-1)] + sub       // substitution
         );
         if ( i>1 && j>1 && a[i-1] == b[j-2] && a[i-2] == b[j-1] ){
            d[i + I*j] = min2(d[i + I*j],d[i-2+I*(j-2)]) + tran; // transposition
         }
      }
   }
   return(d[I*J-1]);
}

//----- R interface ---------------------------------

int vmax(int *x, int n){
   double m = x[0];
   for ( int i = 1; i < n; ++i ){
      if ( x[i] > m ){ 
         m = x[i];
      }
   }
   return(m);
}

SEXP R_osa(SEXP A, SEXP B, SEXP ncharA, SEXP ncharB, SEXP w){
   PROTECT(A);
   PROTECT(B);
   PROTECT(ncharA);
   PROTECT(ncharB);
   PROTECT(w);
   int t, dsize, NA = length(A), NB = length(B);
   double *workspace; 

   // determine workspace size 
   dsize = (vmax(INTEGER(ncharA),NA)+1) * (vmax(INTEGER(ncharB),NB)+1);
   workspace = calloc(dsize, sizeof(double)); 
   if ( workspace == NULL ){
      error("%s\n","unable to allocate enough memory for workspace");
   }

   // output vector
   t = (NA > NB) ? NA : NB;   
   SEXP out;
   PROTECT(out = allocVector(INTSXP,t));
   
   int k,l;
   for ( int i=0; i < t; ++i ){
      k = i % NA;
      l = i % NB; 
      INTEGER(out)[i] = osa(
         CHAR(STRING_ELT(A,k)), 
         INTEGER(ncharA)[k], 
         CHAR(STRING_ELT(B,l)), 
         INTEGER(ncharB)[l], 
         REAL(w), 
         workspace
      );
   }
   
   free(workspace);
   UNPROTECT(6);
   return(out);
}



/* Simple test
void main(){
   int na=3, nb=2;
   char a[3] = "abc";
   char b[2] = "b0";
   int *d;
   d = (int *) calloc((na+1)*(nb+1), sizeof(int));
   for ( int i=0; i < (na+1)*(nb+1); printf("%d, ",d[i++]));
   printf("\n");
   int x = osa(a,na,b,nb,d);
   printf("Distance: %d\n",x);
   free(d);
}

*/

