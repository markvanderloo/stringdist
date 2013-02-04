#include <R.h>
#include <Rdefines.h>


int get_dp_matrix_size(SEXP a, SEXP b){
   int max_a=0, max_b=0, t;
   for ( int i=0; i<length(a); ++i){
      t = length(VECTOR_ELT(a,i));
      if ( max_a < t ) max_a = t;
   }
   for (int i=0; i<length(b); ++i){
      t = length(VECTOR_ELT(b,i));
      if ( max_b < t ) max_b = t;
   }
   return ((max_a + 1)*(max_b + 1));
}

unsigned int max_length(SEXP x){
  unsigned int t=0, m;
  for (int i=0; i<length(x); ++i){
    m = length(VECTOR_ELT(x,i));
    if (t < m) t = m;
  }
  return t;
}

double min3(double x, double y, double z){
   if ( x <= y && x <= z ){ 
      return(x);
   } else if ( y <= x && y <= z ) {
      return(y);
   } else {
      return(z);
   }
}

double min2(double x, double y){
   if ( x <= y ){
      return(x);
   } else {
      return(y);
   }
}

