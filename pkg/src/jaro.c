
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
static int match_int(int a, int *b, int *guard, int width){

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

// Winkler's l-factor (nr of matching characters at beginning of the string).
static double get_l(int *a, int *b, int n){
  int i=0;
  double l;
  while ( a[i] == b[i] && i < n ){ 
    i++;
  }
  l = (double) i;
  return l;
}


/* jaro distance (see http://en.wikipedia.org/wiki/Jaro%E2%80%93Winkler_distance).
 *
 * a    : string (in int rep)
 * b    : string (in int rep)
 * x    : length of a (in uints)
 * y    : length of b (in uints)
 * p    : Winkler's p-factor in (0,0.25)
 * work : workspace, minimally of length max(x,y)
 *
 */
static double jaro_winkler(
             int *a, 
             int *b,
             int x,
             int y,
             double p,
             int *work
        ){

  // edge case
  if ( x == 0 && y == 0 ) return 0;

  // swap arguments if necessary, so we always loop over the shortest string
  if ( x > y ){
    int *c = b;
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

  // Winkler's penalty factor
  if ( p > 0 && d > 0 ){
    int n = min(min(x,y),4);
    d =  d - get_l(a,b,n)*p*d; 
  }

  return d;
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
  int i=0,j=0,n;
  double pp = REAL(p)[0];
  for ( int k=0; k < nt; ++k ){
    length_s = length(VECTOR_ELT(a,i));
    length_t = length(VECTOR_ELT(b,j));
    s = INTEGER(VECTOR_ELT(a,i));
    t = INTEGER(VECTOR_ELT(b,j));
    if ( s[0] == NA_INTEGER || t[0] == NA_INTEGER){
      y[k] = NA_REAL;
      continue;
    } else { // jaro-winkler distance
      y[k] = jaro_winkler(s, t, length_s, length_t, pp, work);
    } 
    i = RECYCLE(i+1,na);
    j = RECYCLE(j+1,nb);
  }
    
  UNPROTECT(3);
  free(work);
  return yy;
}


//-- Match function interface with R
/*
SEXP R_match_jaro_winkler(SEXP x, SEXP table, SEXP nomatch, SEXP matchNA, SEXP p){
  PROTECT(x);
  PROTECT(table);
  PROTECT(nomatch);
  PROTECT(matchNA);
  PROTECT(weight);
  PROTECT(maxDistance);

  int nx = length(x), ntable = length(table);
  int no_match = INTEGER(nomatch)[0];
  int match_na = INTEGER(matchNA)[0];
  double *w = REAL(weight);
  double maxDist = REAL(maxDistance)[0];
  
  // claim space for workhorse 
  int max_x = max_length(x);
  int max_table = max_length(table);
  double *scores = (double *) malloc( (max_x + 3) * (max_table + 2) * sizeof(double) );
  if ( scores == NULL ){
     error("%s\n","unable to allocate enough memory");
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
        d = osa(
          (unsigned int *) X, 
          length(VECTOR_ELT(x,i)), 
          (unsigned int *) T, 
          length(VECTOR_ELT(table,j)), 
          w,
          maxDist,
          scores
        );
        if ( d > -1 && d < d1){ 
          index = j + 1;
          d1 = d;
        }
      } else if ( xNA && tNA ) {  // both are NA
        index = match_na ? j + 1 : no_match;
        break;
      }
    }
    
    y[i] = index;
  }  
  UNPROTECT(7);
  free(scores);
  return(yy);
}
*/

