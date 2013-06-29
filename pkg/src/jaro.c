
#include <R.h>
#include <Rdefines.h>
#include "utils.h"
#include <string.h>

static inline int max(int x, int y){
  int m = (x < y) ? y : x;
  return m;
}

static inline int min(int x, int y){
  int m = (x < y) ? x : y;
  return m; 
}


/* First match of a in b[] 
 * Returns -1 if no match is found
 * Parameter 'guard; indicates which elements of b have been matched before to avoid
 * matching two instances of the same character to the same position in b (which we treat read-only).
 */
static int match_int(unsigned int a, unsigned int *b, int *guard, int width){

  int i = 0;
  while ( 
      b[i] && 
      ( i < width ) && 
      ( b[i] != a || (b[i] == a && guard[i])) 
  ){
    ++i;
  }
  if ( b[i] == a ){
    guard[i] = 1;
    return i;
  } 
  return -1;
}


/* jaro distance (see http://en.wikipedia.org/wiki/Jaro%E2%80%93Winkler_distance).
 *
 * a    : string (in uint rep)
 * b    : string (in uint rep)
 * x    : length of a (in uints)
 * y    : length of b (in uints)
 * work : workspace, minimally of length max(x,y)
 *
 */
static double jaro(
             unsigned int *a, 
             unsigned int *b,
             int x,
             int y,
             int *work
        ){

  // edge case
  if ( x == 0 && y == 0 ) return 0;

  // swap arguments if necessary, so we always loop over the shortest string
  if ( x > y ){
    unsigned int *c = b;
    int z = y;
    b = a;
    a = c;
    y = x;
    x = z;
  }

  // max transposition distance
  int M = max(max(x,y)/2 - 1,0);
  // transposition counter
  double t = 0.0;
  // number of matches 
  double m = 0.0;
  
  int left, right, J, jmax=0;
  
  for ( int i=0; i < x; ++i ){
    left  = max(0, i-M);

    if ( left >= y ){
      J = -1;
    } else {
      right = min(y, i+M);
      J =  match_int(a[i], b + left, work + left, right - left);
    }
    if ( J >= 0 ){
      ++m;
      t += (J + left < jmax ) ? 1 : 0; 
      jmax = max(jmax, J + left);
    }
  }
  double d;
  if ( m < 1 ){
    d = 1.0;
  } else {
    d = 1.0 - (1.0/3.0)*(m/x + m/y + (m-t)/m);
  }
  memset(work,0,sizeof(int) * y);
  return d;
}

// Winkler's l-factor (nr of matching characters at beginning of the string).
static double get_l(unsigned int *a, unsigned int *b, int n){
  int i;
  double l;
  while ( a[i] == b[i] && i < n ){ 
    i++;
  }
  l = (double) i;
  return l;
}


/*----------- R interface ------------------------------------------------*/

SEXP R_jaro_winkler(SEXP a, SEXP b, SEXP p){
  PROTECT(a);
  PROTECT(b);

  // find the length of longest strings
  int max_a = max_length(a);
  int max_b = max_length(b);
  int max_char = max(max_a,max_b);
  // recycling parameters
  int na = length(a);
  int nb = length(b);
  int nt = max(na,nb);
  
  int *s, *t;
  int length_s, length_t;

  // workspace for worker function
  int *work = (int *) calloc( sizeof(int), max_char );

  // output variable
  SEXP yy;
  PROTECT(yy = allocVector(REALSXP,nt));
  double *y = REAL(yy);

  // compute distances, skipping NA's
  int i=0,j=0,l,n;
  double pp = REAL(p)[0];
  for ( int k=0; k < nt; ++k ){
    length_s = length(VECTOR_ELT(a,i));
    length_t = length(VECTOR_ELT(b,j));
    s = (unsigned int *) INTEGER(VECTOR_ELT(a,i));
    t = (unsigned int *) INTEGER(VECTOR_ELT(b,j));
    if ( s[0] == NA_INTEGER || t[0] == NA_INTEGER){
      y[k] = NA_REAL;
      continue;
    } else { // jaro distance
      y[k] = jaro(s, t, length_s, length_t, work);
    } 
    // Winkler's penalty factor
    if ( pp > 0 && y[k] != NA_REAL && y[k] > 0 ){
      n = min(min(length_s,length_t),4);
      y[k] =  y[k] - get_l(s,t,n)*pp*y[k]; 

    }
    i = RECYCLE(i+1,na);
    j = RECYCLE(j+1,nb);
  }
    

  UNPROTECT(3);
  free(work);
  return yy;
}




