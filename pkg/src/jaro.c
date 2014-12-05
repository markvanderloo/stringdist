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
#ifdef _OPENMP
#include <omp.h>
#endif



/* First match of a in b[] 
 * Returns -1 if no match is found
 * Parameter 'guard; indicates which elements of b have been matched before to avoid
 * matching two instances of the same character to the same position in b (which we treat read-only).
 */
static int match_int(unsigned int a, unsigned int *b, int *guard, int width, int m){
  int i = 0;
  while ( 
      ( i < width ) && 
      ( b[i] != a || (b[i] == a && guard[i])) 
  ){
    ++i;
  }
  // ugly edge case workaround
  if ( !(m && i==width) && b[i] == a ){
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
             double *w,
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
  memset(work,0, sizeof(int) * y);

  // max transposition distance
  int M = MAX(MAX(x,y)/2 - 1,0);
  // transposition counter
  double t = 0.0;
  // number of matches 
  double m = 0.0;
  int max_reached; 
  int left, right, J, jmax=0;
  
  for ( int i=0; i < x; ++i ){
    left  = MAX(0, i-M);
    if ( left >= y ){
      J = -1;
    } else {
      right = MIN(y, i+M);
      // ugly workaround: I should rewrite match_int.
      max_reached = (right == y) ? 1 : 0;
      J =  match_int(a[i], b + left, work + left, right - left, max_reached);
    }

    if ( J >= 0 ){
      ++m;
      t += (J + left < jmax ) ? 1 : 0; 
      jmax = MAX(jmax, J + left);
    }
  }
  double d;
  if ( m < 1 ){
    d = 1.0;
  } else {
    d = 1.0 - (1.0/3.0)*(w[0]*m/x + w[1]*m/y + w[2]*(m-t)/m);
  }

  // Winkler's penalty factor
  if ( p > 0 && d > 0 ){
    int n = MIN(MIN(x,y),4);
    d =  d - get_l(a,b,n)*p*d; 
  }

  return d;
}


/*----------- R interface ------------------------------------------------*/

SEXP R_jw(SEXP a, SEXP b, SEXP p, SEXP weight, SEXP useBytes, SEXP nthrd){
  PROTECT(a);
  PROTECT(b);
  PROTECT(p);
  PROTECT(weight);
  PROTECT(useBytes);
  PROTECT(nthrd);

  // find the length of longest strings
  int ml_a = max_length(a)
    , ml_b = max_length(b)
    , na = length(a)
    , nb = length(b)
    , nt = MAX(na,nb)
    , bytes = INTEGER(useBytes)[0];
 
  // output variable
  SEXP yy;
  PROTECT(yy = allocVector(REALSXP,nt));
  double *y = REAL(yy)
    , pp = REAL(p)[0]
    , *w = REAL(weight);


  #ifdef _OPENMP 
  int  nthreads = INTEGER(nthrd)[0];
  #pragma omp parallel num_threads(nthreads) default(none) \
      shared(y, w, pp, R_PosInf, NA_REAL, bytes, na, nb, ml_a, ml_b, nt, a, b)
  #endif
  {
    // workspace for worker function
    int *work = (int *) malloc( sizeof(int) * MAX(ml_a,ml_b) );
    unsigned int *s = NULL, *t = NULL;
    s = (unsigned int *) malloc((2L + ml_a + ml_b) * sizeof(int));
    t = s + ml_a + 1L;
    if ( (work == NULL) | (bytes && s == NULL) ) nt = -1;


    // compute distances, skipping NA's
    int len_s, len_t, isna_s, isna_t
      , i = 0, j = 0, ID = 0, num_threads = 1;

    #ifdef _OPENMP
    ID = omp_get_thread_num();
    num_threads = omp_get_num_threads();
    i = recycle(ID-num_threads, num_threads, na);
    j = recycle(ID-num_threads, num_threads, nb);
    #endif
    for ( int k = ID; k < nt; k += num_threads){
      get_elem1(a, i, bytes, &len_s, &isna_s, s);
      get_elem1(b, j, bytes, &len_t, &isna_t, t);
      if ( isna_s || isna_t ){
        y[k] = NA_REAL;
        continue;
      } else { // jaro-winkler distance
        y[k] = jaro_winkler(s, t, len_s, len_t, pp, w, work);
      }
      i = recycle(i, num_threads, na);
      j = recycle(j, num_threads, nb); 
    }
      
    free(work);
    if (bytes) free(s);
  }
  UNPROTECT(6);
  if ( nt < 0 ) error("Unable to allocate enough memory");
  return yy;
}


//-- Match function interface with R

SEXP R_match_jw(SEXP x, SEXP table, SEXP nomatch, SEXP matchNA, SEXP p
    , SEXP weight, SEXP maxDist, SEXP nthrd){
  PROTECT(x);
  PROTECT(table);
  PROTECT(nomatch);
  PROTECT(matchNA);
  PROTECT(p);
  PROTECT(weight);
  PROTECT(maxDist);
  PROTECT(nthrd);

  int nx = length(x)
    , ntable = length(table)
    , no_match = INTEGER(nomatch)[0]
    , match_na = INTEGER(matchNA)[0]
    , bytes = IS_CHARACTER(x)
    , ml_x = max_length(x)
    , ml_t = max_length(table)
    , nthreads = INTEGER(nthrd)[0];

  double pp = REAL(p)[0]
    , *w = REAL(weight)
    , max_dist = REAL(maxDist)[0] == 0.0 ? R_PosInf : REAL(maxDist)[0];
  
  // output vector
  SEXP yy;
  PROTECT(yy = allocVector(INTSXP, nx));
  int *y = INTEGER(yy);


  #ifdef _OPENMP
  #pragma omp parallel num_threads(nthreads) default(none) \
    shared(x,table, y, R_PosInf, nx, ntable, no_match, match_na, bytes,ml_x,ml_t,pp,w, max_dist)
  #endif
  {
    // workspace for worker function
    int *work = (int *) malloc( sizeof(int) * MAX(ml_x, ml_t) );
    unsigned int *X = NULL, *T = NULL;
    if (bytes){
      X = (unsigned int *) malloc( (ml_x + ml_t) * sizeof(int));
      T = X + ml_x;
    }
    if ( (work == NULL) | (bytes && X == NULL) ) nx = -1;

    double d = R_PosInf, d1 = R_PosInf;
    int index, isna_X, isna_T, len_X,len_T;


    #ifdef _OPENMP
    #pragma omp for
    #endif
    for ( int i=0; i<nx; i++){
      index = no_match;
      X = get_elem(x, i, bytes, &len_X, &isna_X, X);
      d1 = R_PosInf;
      for ( int j=0; j<ntable; j++){
        T = get_elem(table, j, bytes, &len_T, &isna_T, T);

        if ( !isna_X && !isna_T ){        // both are char (usual case)
          d = jaro_winkler(X, T, len_X, len_T, pp, w, work);
          if ( d > max_dist ){
            continue;
          } else if ( d > -1.0 && d < d1){ 
            index = j + 1;
            if ( ABS(d) < 1e-14 ) break; // exact match
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
  } // end of parallel region
  UNPROTECT(8);
  if ( nx < 0 ) error ("Unable to allocate enough memory");
  return(yy);
}


