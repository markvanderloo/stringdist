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
 *
 *
 * This code is gratefully based on Nick Logan's github repository
 * https://github.com/ugexe/Text--Levenshtein--Damerau--XS/blob/master/damerau-int.c
 *
 */ 



#define USE_RINTERNALS
#include <R.h>
#include <Rdefines.h>
#include "utils.h"
#ifdef _OPENMP
#include <omp.h>
#endif

static double hamming_dist(unsigned int *a, int len_a, unsigned int *b, int len_b){
  double h=0;
    if (len_a != len_b) return 1.0/0.0;
    for(int i=0; i < len_a; ++i){
     if (a[i] != b[i]) h++;
    }
  return h;
}

// -- R interface

SEXP R_hm(SEXP a, SEXP b, SEXP useBytes, SEXP nthrd){
  PROTECT(a);
  PROTECT(b);
  PROTECT(useBytes);
  PROTECT(nthrd);

  int na = length(a)
    , nb = length(b)
    , nt = ( na > nb) ? na : nb
    , bytes = INTEGER(useBytes)[0]
    , ml_a = max_length(a)
    , ml_b = max_length(b);

  // create answer vector.
  SEXP yy;
  PROTECT(yy = allocVector(REALSXP,nt));
  double *y = REAL(yy);

  /* formally, the pragma statement need not be included in the #ifdef, but 
   * at least in gcc it generates warning (hence CRAN-trouble) when not 
   * recognized.
   */
  #ifdef _OPENMP 
  int  nthreads = INTEGER(nthrd)[0];
  #pragma omp parallel num_threads(nthreads) default(none) \
      shared(y, R_PosInf, NA_REAL, bytes, na, nb, ml_a, ml_b, nt, a, b)
  #endif
  {
    unsigned int *s = NULL, *t = NULL;
    s = (unsigned int *) malloc( (2L + ml_a + ml_b) * sizeof(int));
    if ( s == NULL ) error("Unable to allocate enough memory");
    t = s + ml_a + 1L;
    
    int k, len_s, len_t, isna_s, isna_t
      , i = 0, j = 0, ID = 0, num_threads = 1;

    #ifdef _OPENMP
    ID = omp_get_thread_num();
    num_threads = omp_get_num_threads();
    i = recycle(ID-num_threads, num_threads, na);
    j = recycle(ID-num_threads, num_threads, nb);
    #endif
    for ( k = ID; k < nt; k += num_threads ){
      get_elem1(a, i, bytes, &len_s, &isna_s, s);
      get_elem1(b, j, bytes, &len_t, &isna_t, t);
      if ( isna_s || isna_t ){
        y[k] = NA_REAL;
      } else if ( len_s != len_t ){
          y[k] = R_PosInf;
      } else {
        y[k] = hamming_dist(s, len_s, t, len_t);
      }
      i = recycle(i, num_threads, na);
      j = recycle(j, num_threads, nb);
    }
    free(s);
  } // end of parallel region
  UNPROTECT(5);
  return yy;
}


//-- Match function interface with R

SEXP R_match_hm(SEXP x, SEXP table, SEXP nomatch, SEXP matchNA
    , SEXP maxDistance, SEXP useBytes, SEXP nthrd){
  PROTECT(x);
  PROTECT(table);
  PROTECT(nomatch);
  PROTECT(matchNA);
  PROTECT(useBytes);
  PROTECT(nthrd);

  int nx = length(x)
    , ntable = length(table)
    , no_match = INTEGER(nomatch)[0]
    , match_na = INTEGER(matchNA)[0]
    , max_dist = INTEGER(maxDistance)[0]
    , bytes = INTEGER(useBytes)[0];

  // convert to integer. 
  Stringset *X = new_stringset(x, bytes);
  Stringset *T = new_stringset(table, bytes);

  // output vector
  SEXP yy;
  PROTECT(yy = allocVector(INTSXP, nx));
  int *y = INTEGER(yy);

  #ifdef _OPENMP
  int nthreads = INTEGER(nthrd)[0];
  #pragma omp parallel num_threads(nthreads) default(none) \
    shared(X, T, y, R_PosInf, NA_INTEGER, nx, ntable, no_match, match_na, bytes, max_dist)
  #endif
  {
    double d = R_PosInf, d1 = R_PosInf;
    int index, len_X, len_T;
    unsigned int *str, **tab;

    #ifdef _OPENMP
    #pragma omp for
    #endif
    for ( int i=0; i<nx; i++){
      index = no_match;
      len_X = X->str_len[i];
      d1 = R_PosInf;
      str = X->string[i];
      tab = T->string;
      for ( int j=0; j<ntable; j++, tab++){
        len_T = T->str_len[j];
        if ( len_X != len_T ) continue;
        if ( len_X != NA_INTEGER && len_T != NA_INTEGER ){        // both are char (usual case)
          d = hamming_dist( str, len_X, *tab, len_T );
          if ( d <= max_dist && d < d1){ 
            index = j + 1;
            if ( d == 0.0 ) break;
            d1 = d;
          }
        } else if ( len_X == NA_INTEGER && len_T == NA_INTEGER ) {  // both are NA
          index = match_na ? j + 1 : no_match;
          break;
        }
      }
      
      y[i] = index;
    }  
  } // end of parallel region
  free_stringset(X);
  free_stringset(T);
  UNPROTECT(7);
  if (nx < 0 ) error("Unable to allocate enough memory");
  return(yy);
}

