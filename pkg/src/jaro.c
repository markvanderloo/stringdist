/*  stringdist - a C library of string distance algorithms with an interface to R.
 *  Copyright (C) 2013  Mark van der Loo
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 *
 *  You can contact the author at: mark _dot_ vanderloo _at_ gmail _dot_ com
 */

#define USE_RINTERNALS
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

// Winkler's l-factor (nr of matching characters at beginning of the string).
static double get_l(unsigned int *a, unsigned int *b, int n){
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
             unsigned int *a, 
             unsigned int *b,
             int x,
             int y,
             double p,
             int *work
        ){

  // edge case
  if ( x == 0 && y == 0 ) return 0;

  // swap arguments if necessary, so we always loop over the shortest string
  if ( x > y ){
    unsigned int *c = b;
    unsigned int z = y;
    b = a;
    a = c;
    y = x;
    x = z;
  }

  memset(work,0,sizeof(int) * y);

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

  // Winkler's penalty factor
  if ( p > 0 && d > 0 ){
    int n = min(min(x,y),4);
    d =  d - get_l(a,b,n)*p*d; 
  }

  return d;
}


/*----------- R interface ------------------------------------------------*/

SEXP R_jw(SEXP a, SEXP b, SEXP p){
  PROTECT(a);
  PROTECT(b);
  PROTECT(p);

  // find the length of longest strings
  int ml_a = max_length(a)
    , ml_b = max_length(b)
    , na = length(a)
    , nb = length(b)
    , nt = max(na,nb)
    , bytes = IS_CHARACTER(a);
  
  double pp = REAL(p)[0];

  // workspace for worker function
  int *work = (int *) malloc( sizeof(int) * max(ml_a,ml_b) );
  unsigned int *s = NULL, *t = NULL;
  if (bytes){
    s = (unsigned int *) malloc((ml_a + ml_b) * sizeof(int));
    t = s + ml_a;
  }
  if ( (work == NULL) | (bytes && s == NULL) ){
     UNPROTECT(3); free(s); free(work);
     error("Unable to allocate enough memory");
  }

  // output variable
  SEXP yy;
  PROTECT(yy = allocVector(REALSXP,nt));
  double *y = REAL(yy);

  // compute distances, skipping NA's
  int i=0,j=0, len_s, len_t, isna_s, isna_t;

  for ( int k=0; k < nt; 
      ++k 
    , i = RECYCLE(i+1,na)
    , j = RECYCLE(j+1,nb) ){

    s = get_elem(a, i, bytes, &len_s, &isna_s, s);
    t = get_elem(b, j, bytes, &len_t, &isna_t, t);
    if ( isna_s || isna_t ){
      y[k] = NA_REAL;
      continue;
    } else { // jaro-winkler distance
      y[k] = jaro_winkler(s, t, len_s, len_t, pp, work);
    } 
  }
    
  free(work);
  if (bytes) free(s);
  UNPROTECT(4);
  return yy;
}


//-- Match function interface with R

SEXP R_match_jw(SEXP x, SEXP table, SEXP nomatch, SEXP matchNA, SEXP p, SEXP maxDist){
  PROTECT(x);
  PROTECT(table);
  PROTECT(nomatch);
  PROTECT(matchNA);
  PROTECT(p);
  PROTECT(maxDist);

  int nx = length(x)
    , ntable = length(table)
    , no_match = INTEGER(nomatch)[0]
    , match_na = INTEGER(matchNA)[0]
    , bytes = IS_CHARACTER(x)
    , ml_x = max_length(x)
    , ml_t = max_length(table);

  double pp = REAL(p)[0];
  double max_dist = REAL(maxDist)[0] == 0.0 ? R_PosInf : REAL(maxDist)[0];
  
  // workspace for worker function
  int *work = (int *) malloc( sizeof(int) * max(ml_x, ml_t) );
  unsigned int *X = NULL, *T = NULL;
  if (bytes){
    X = (unsigned int *) malloc( (ml_x + ml_t) * sizeof(int));
    T = X + ml_x;
  }
  if ( (work == NULL) | (bytes && X == NULL) ){
    UNPROTECT(6); free(work); free(X);
    error ("Unable to allocate enough memory\n");
  }

  // output vector
  SEXP yy;
  PROTECT(yy = allocVector(INTSXP, nx));
  int *y = INTEGER(yy);


  double d = R_PosInf, d1 = R_PosInf;
  int index, isna_X, isna_T, len_X,len_T;


  for ( int i=0; i<nx; i++){
    index = no_match;
    X = get_elem(x, i, bytes, &len_X, &isna_X, X);
    d1 = R_PosInf;
    for ( int j=0; j<ntable; j++){
      T = get_elem(table, j, bytes, &len_T, &isna_T, T);

      if ( !isna_X && !isna_T ){        // both are char (usual case)
        d = jaro_winkler(X, T, len_X, len_T, pp, work);
        if ( d > max_dist ){
          continue;
        } else if ( d > -1 && d < d1){ 
          index = j + 1;
          if ( abs(d) < 1e-14 ) break;
          d1 = d;
        }
      } else if ( isna_X && isna_T ) {  // both are NA
        index = match_na ? j + 1 : no_match;
        break;
      }
    }
    
    y[i] = index;
  }   

  if (bytes) free(X); 
  free(work);
  UNPROTECT(7);
  return(yy);
}


