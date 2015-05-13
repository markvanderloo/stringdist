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
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include "utils.h"
#include "stringdist.h"
#ifdef _OPENMP
#include <omp.h>
#endif
// TODO: catch error and report.
static Stringdist *R_open_stringdist(Distance d, int max_len_a, int max_len_b, SEXP weight, SEXP p, SEXP q){

  Stringdist *sd = NULL;
  if (d == osa || d == lv || d == dl || d == hamming || d == lcs){
    sd = open_stringdist(d, max_len_a, max_len_b, REAL(weight));
  } else if ( d == qgram || d == cosine || d == jaccard ){
    sd = open_stringdist(d, max_len_a, max_len_b, (unsigned int) INTEGER(q)[0]);
  } else if ( d == jw ){
    sd = open_stringdist(d, max_len_a, max_len_b, REAL(weight), REAL(p)[0]);
  } else if (d == soundex) {
    sd = open_stringdist(d, max_len_a, max_len_b);
  }

  return sd;
}



SEXP R_stringdist(SEXP a, SEXP b, SEXP method
  , SEXP weight, SEXP p, SEXP q
  , SEXP useBytes, SEXP nthrd){
  PROTECT(a);
  PROTECT(b);
  PROTECT(method);
  PROTECT(weight);
  PROTECT(p);
  PROTECT(q);
  PROTECT(useBytes);
  PROTECT(nthrd);

  int na = length(a)
    , nb = length(b)
    , bytes = INTEGER(useBytes)[0]
    , ml_a = max_length(a)
    , ml_b = max_length(b)
    , nt = (na > nb) ? na : nb;
  
  // output vector
  SEXP yy;
  PROTECT(yy = allocVector(REALSXP, nt));
  double *y = REAL(yy);


  #ifdef _OPENMP 
  int  nthreads = INTEGER(nthrd)[0];
  #pragma omp parallel num_threads(nthreads) default(none) \
      shared(y,na,nb, R_PosInf, NA_REAL, bytes, method, weight, p, q, ml_a, ml_b, nt, a, b)
  #endif
  {

    Stringdist *sd = R_open_stringdist( (Distance) INTEGER(method)[0]
        , ml_a, ml_b
        , weight
        , p
        , q
    );

    unsigned int *s = NULL, *t = NULL;
    s = (unsigned int *) malloc(( 2L + ml_a + ml_b) * sizeof(int));

    if ( (sd==NULL) | (bytes && s == NULL) ) nt = -1;
    t = s + ml_a + 1L;
      
    int len_s, len_t, isna_s, isna_t
      , i = 0, j = 0, ID = 0, num_threads = 1;

    #ifdef _OPENMP
    ID = omp_get_thread_num();
    num_threads = omp_get_num_threads();
    i = recycle(ID-num_threads, num_threads, na);
    j = recycle(ID-num_threads, num_threads, nb);
    #endif
    for ( int k=ID; k < nt; k += num_threads ){
      get_elem1(a, i, bytes, &len_s, &isna_s, s);
      get_elem1(b, j, bytes, &len_t, &isna_t, t);
      if (isna_s || isna_t){
        y[k] = NA_REAL;
      } else {
        y[k] = stringdist(sd, s, len_s, t, len_t);
        if ( y[k] < 0 ) y[k] = R_PosInf;
      }
      i = recycle(i, num_threads, na);
      j = recycle(j, num_threads, nb);
    }
    

    close_stringdist(sd);

    free(s);
  } // end of parallel region

  UNPROTECT(9);
  if (nt < 0 ) error("Unable to allocate enough memory");
  return(yy);
}

/* amatch
 *
 */
SEXP R_amatch(SEXP x, SEXP table, SEXP method 
  , SEXP nomatch, SEXP matchNA
  , SEXP weight, SEXP p, SEXP q
  , SEXP maxDistance, SEXP useBytes
  , SEXP nthrd){

  PROTECT(x);
  PROTECT(table);
  PROTECT(method);
  PROTECT(nomatch);
  PROTECT(matchNA);
  PROTECT(weight);
  PROTECT(p);
  PROTECT(q);
  PROTECT(maxDistance);
  PROTECT(useBytes);
  PROTECT(nthrd);

  int nx = length(x)
    , ntable = length(table)
    , no_match = INTEGER(nomatch)[0]
    , match_na = INTEGER(matchNA)[0]
    , bytes = INTEGER(x)[0]
    , ml_x = max_length(x)
    , ml_t = max_length(table);


  double maxDist = REAL(maxDistance)[0];

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
    shared(X, T, y, R_PosInf, NA_INTEGER, nx, ntable, no_match, match_na, bytes, ml_x, ml_t, method, weight, p, q, maxDist)
  #endif
  {
    /* claim space for workhorse */

    Stringdist *sd = R_open_stringdist( (Distance) INTEGER(method)[0]
        , ml_x, ml_t
        , weight
        , p
        , q
    );
    double d = R_PosInf, d1 = R_PosInf;
    int index, len_X, len_T;
    unsigned int *str;
    unsigned int **tab;

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
        if (len_X != NA_INTEGER && len_T != NA_INTEGER ){        // both are char (usual case)
          d = stringdist(sd, str, len_X, *tab, len_T);
          if ( d <= maxDist && d < d1){ 
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
    close_stringdist(sd);
  } // end of parallel region
  free_stringset(X);
  free_stringset(T);
  UNPROTECT(12);
  if (nx < 0 ) error("Unable to allocate enough memory");
  return(yy);
} // end R_amatch






