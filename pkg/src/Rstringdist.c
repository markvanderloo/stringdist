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
#include <stdint.h>
#include <R.h>
#include <Rdefines.h>
#include <math.h>
#include "utils.h"
#include "stringdist.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))


static Stringdist *R_open_stringdist(Distance d, int max_len_a, int max_len_b, SEXP weight, SEXP p, SEXP bt, SEXP q){

  Stringdist *sd = NULL;
  if (d == osa || d == lv || d == dl || d == hamming || d == lcs){
    sd = open_stringdist(d, max_len_a, max_len_b, REAL(weight));
  } else if (d == qgram || d == cosine || d == jaccard || d == running_cosine){
    sd = open_stringdist(d, max_len_a, max_len_b, (unsigned int) INTEGER(q)[0]);
  } else if ( d == jw ){
    sd = open_stringdist(d, max_len_a, max_len_b, REAL(weight), REAL(p)[0], REAL(bt)[0]);
  } else if (d == soundex) {
    sd = open_stringdist(d, max_len_a, max_len_b);
  }
  if ( sd == NULL ){
    error("Could not allocate enough memory");
  }
  return sd;
}



SEXP R_stringdist(SEXP a, SEXP b, SEXP method
  , SEXP weight, SEXP p, SEXP bt, SEXP q
  , SEXP useBytes, SEXP nthrd){

  int na = length(a)
    , nb = length(b)
    , bytes = INTEGER(useBytes)[0]
    , ml_a = max_length(a)
    , ml_b = max_length(b)
    , nt = (na > nb) ? na : nb
    , intdist = TYPEOF(a) == VECSXP ? 1 : 0; // expect lists of integers? 
 
  // output vector
  SEXP yy;
  PROTECT(yy = allocVector(REALSXP, nt));
  double *y = REAL(yy);

  #ifdef _OPENMP 
  int  nthreads = MIN(INTEGER(nthrd)[0],MAX(na,nb));
  #pragma omp parallel num_threads(nthreads) default(none) \
      shared(y,na,nb, R_PosInf, NA_REAL, bytes, intdist, method, weight, p, bt, q, ml_a, ml_b, nt, a, b)
  #endif
  {

    Stringdist *sd = R_open_stringdist( (Distance) INTEGER(method)[0]
        , ml_a, ml_b
        , weight
        , p
        , bt
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
      get_elem(a, i, bytes, intdist, &len_s, &isna_s, s);
      get_elem(b, j, bytes, intdist, &len_t, &isna_t, t);
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

  UNPROTECT(1);
  if (nt < 0 ) error("Unable to allocate enough memory");
  return(yy);
}

/* amatch
 *
 */
SEXP R_amatch(SEXP x, SEXP table, SEXP method 
  , SEXP nomatch, SEXP matchNA
  , SEXP weight, SEXP p, SEXP bt, SEXP q
  , SEXP maxDistance, SEXP useBytes
  , SEXP nthrd){


  int nx = length(x)
    , ntable = length(table)
    , no_match = INTEGER(nomatch)[0]
    , match_na = INTEGER(matchNA)[0]
    , bytes = INTEGER(x)[0]
    , ml_x = max_length(x)
    , ml_t = max_length(table)
    , intdist = TYPEOF(x) == VECSXP ? 1 : 0; // list of integers?


  double maxDist = REAL(maxDistance)[0];

  // convert to integer. 
  Stringset *X = new_stringset(x, bytes, intdist);
  Stringset *T = new_stringset(table, bytes, intdist);

  // output vector
  SEXP yy;
  PROTECT(yy = allocVector(INTSXP, nx));
  int *y = INTEGER(yy);
  
  #ifdef _OPENMP
  int nthreads = MAX(MIN(INTEGER(nthrd)[0],nx),0);
  #pragma omp parallel num_threads(nthreads) default(none) \
    shared(X, T, y, R_PosInf, NA_INTEGER, nx, ntable, no_match, match_na, ml_x, ml_t, method, weight, p, bt, q, maxDist)
  #endif
  {
    /* claim space for workhorse */

    Stringdist *sd = R_open_stringdist( (Distance) INTEGER(method)[0]
        , ml_x, ml_t
        , weight
        , p
        , bt
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
            if ( fabs(d) < 1e-14 ){ 
              break; // exact match
            }
            d1 = d;
          }
        } else if ( len_X == NA_INTEGER && len_T == NA_INTEGER ) {  // both are NA
          index = match_na ? j + 1 : no_match;
          break;
        }
      }
      y[i] = index;
    }
    str=NULL;
    tab=NULL;
    close_stringdist(sd);
  } // end of parallel region
  free_stringset(X);
  free_stringset(T);
  UNPROTECT(1);

  return(yy);
} // end R_amatch


// Lower tridiagonal distance matrix for a single vector argument.

static int get_j(R_xlen_t k, int n){
  double nd = (double) n;
  double kd = (double) k;
  double u = ceil( (2.*nd - 3.)/2. - sqrt(pow(nd-.5,2.) - 2.*(kd + 1.)) );

  return (int) u;
}

/* max n for objects of length n(n-1).
 * 
*/
#ifdef LONG_VECTOR_SUPPORT
#define MAXN ( (R_xlen_t) (0.5 + 1.5 * sqrt((double) R_XLEN_T_MAX)) )
#else
#define MAXN ( (R_xlen_t) (0.5 + 1.5 * sqrt((double) R_LEN_T_MAX)) )
#endif

SEXP R_lower_tri(SEXP a, SEXP method
  , SEXP weight, SEXP p,  SEXP bt, SEXP q
  , SEXP useBytes, SEXP nthrd){

  int bytes = INTEGER(useBytes)[0]
    , ml = max_length(a)
    , intdist = TYPEOF(a) == VECSXP ? 1 : 0; // expect list of integer vectors? 

  // Long vectors on platforms where LONG_VECTOR_SUPPORT is defined.
  R_xlen_t n = xlength(a)
    , N = n*(n-1)/2;

  if ( n > MAXN ){
    error("Length of input vector (%d) exceeds maximum allowed for this platform (%d)",n,MAXN);
  }  
 
 
  // output vector
  SEXP yy;
  PROTECT(yy = allocVector(REALSXP, N));
  // nothing to do if n=1 
  if (n == 1L) goto end ;
  double *y = REAL(yy);


  #ifdef _OPENMP 
  int  nthreads = MIN(INTEGER(nthrd)[0],N);
  nthreads = MIN(nthreads, n);
  #pragma omp parallel num_threads(nthreads) default(none) \
      shared(y,n,N, R_PosInf, NA_REAL, bytes, intdist, method, weight, p, bt, q, ml, a)
  #endif
  {

    Stringdist *sd = R_open_stringdist( (Distance) INTEGER(method)[0]
        , ml, ml
        , weight
        , p
        , bt
        , q
    );

    unsigned int *s = NULL, *t = NULL;
    s = (unsigned int *) malloc(( 2L + 2*ml) * sizeof(int));

    if ( (sd==NULL) | (bytes && s == NULL) ) n = -1;
    t = s + ml + 1L;
      
    int len_s, len_t, isna_s, isna_t
      , i = 0, j = 0
      , thread_id = 0
      , n_threads = 1
      , col_max = n-1;

    R_xlen_t pp = 0
      , k_start = 0
      , k_end = N;

    #ifdef _OPENMP
      thread_id  = omp_get_thread_num();
      n_threads = omp_get_num_threads();
    #endif
      // some administration to parallelize the loop.
      pp = N / n_threads;
      k_start = thread_id * pp;
      k_end   = (thread_id < n_threads - 1 ) ? k_start + pp : N;
      j = get_j(k_start,n);
      i = k_start + j * (j - 2*n + 3)/2;
    for ( R_xlen_t k=k_start; k < k_end; k++ ){
      i++;
      get_elem(a, i, bytes, intdist, &len_s, &isna_s, s);
      get_elem(a, j, bytes, intdist, &len_t, &isna_t, t);

      if (isna_s || isna_t){
        y[k] = NA_REAL;
      } else {
        y[k] = stringdist(sd, s, len_s, t, len_t);
        if ( y[k] < 0 ) y[k] = R_PosInf;
      }
      if ( i == col_max ){
        j++;
        i = j;
      }
    }

    free(s);
    close_stringdist(sd);
  } // end of parallel region

  end:
  UNPROTECT(1);
  if (n < 0 ) error("Unable to allocate enough memory");
  return(yy);
}

// afind
// For each string in 'a', return the starting position of
// the best match with 'pattern'.
SEXP R_afind(SEXP a, SEXP pattern, SEXP width
  , SEXP method, SEXP weight, SEXP p, SEXP bt
  , SEXP q, SEXP useBytes, SEXP nthrd)
{
  
  int na = length(a)              // nr of  texts to search
    , npat = length(pattern)      // nr of patterns
    , ml_a   = max_length(a)      // max length of searched string
    , ml_b = max_length(pattern)  // max length of the pattern.
    , intdist = 0                 // no distances between integer sequences (yet)
    , bytes = INTEGER(useBytes)[0];


  int *window = INTEGER(width);   // access the window widths.

  // output list
  SEXP out_list;
  PROTECT(out_list = allocVector(VECSXP, 2));

  // output location
  SEXP out_loc;
  PROTECT(out_loc = allocMatrix(INTSXP, na, npat));
  VECTOR_ELT(out_list,0) = out_loc;
  int *yloc = INTEGER(out_loc);

  // output distance
  SEXP out_dist;
  PROTECT(out_dist = allocMatrix(REALSXP, na, npat));
  VECTOR_ELT(out_list,1) = out_dist;
  double *ydist = REAL(out_dist);
  // Setup stringdist structure.
  // find maximum window length
  int max_window = 0;
  for ( int i=0; i<npat; i++){
    if (max_window < window[i]){
      max_window = window[i];
    }
  }

 
  #ifdef _OPENMP 
  int  nthreads = MIN(INTEGER(nthrd)[0],na);
  #pragma omp parallel num_threads(nthreads) default(none) \
      shared(yloc,ydist, na, npat, R_PosInf, NA_REAL, NA_INTEGER, bytes, intdist, \
      method, weight, p, bt, q, ml_a, ml_b, window, max_window, a, pattern)
  #endif
  {  // start parallel region


    Stringdist *sd = R_open_stringdist( (Distance) INTEGER(method)[0]
        , max_window, ml_b
        , weight
        , p
        , bt
        , q
    );

    // allocate memory to store the strings
    unsigned int *s = NULL, *t = NULL;
    s = (unsigned int *) malloc(( 2L + ml_a + ml_b) * sizeof(int));
    
    // t is the location of the pattern
    t = s + ml_a + 1L;
    
    int len_s, len_t, isna_s, isna_t, max_k, k_min, current_window, offset;
    int ID=0, num_threads=1;
    
    double d, d_min;

    #ifdef _OPENMP
    ID = omp_get_thread_num();
    num_threads = omp_get_num_threads();
    #endif
    for ( int i = ID; i < na; i += num_threads ){
      // get text to search
      get_elem(a, i, bytes, intdist, &len_s, &isna_s, s);
      for( int j = 0; j < npat; j++){
        // get pattern
        get_elem(pattern, j, bytes, intdist, &len_t, &isna_t, t);
        current_window = window[j];
        offset = j*na;
        if (isna_s || isna_t){ // something to search in, or find?
          yloc[offset + i]  = NA_INTEGER;
          ydist[offset + i] = NA_REAL;
        } else if ( current_window >= len_s ){ // is the text shorter than the window?
          yloc[offset + i]  = 1L;
          ydist[offset + i] = stringdist(sd, s, len_s, t, len_t); 
        } else { // slide window over text and compute distances
          max_k = len_s - current_window;
          d_min = R_PosInf;
          k_min = 0;
          for (int k = 0; k <= max_k; k++){
            d = stringdist(sd, s + k, current_window, t, len_t);
            if ( d < d_min ){
              d_min = d;
              k_min = k;
            }
          } // end loop over windows
          yloc[offset + i]  = k_min + 1;
          ydist[offset + i] = d_min;
          reset_stringdist(sd);
        }
      } // end loop over patterns
    } // end loop over strings
    close_stringdist(sd);
  } // end parallel region
  UNPROTECT(3);
  return(out_list);

}

// helper function to determine  whether all is INTSXP

SEXP R_all_int(SEXP X){
  PROTECT(X);
  SEXP all_int;
  all_int = PROTECT(allocVector(LGLSXP,1L));

  int n = length(X);
  LOGICAL(all_int)[0] =  1L;
  for (int i=0; i<n; i++){
    if (TYPEOF(VECTOR_ELT(X,i)) != INTSXP){
      LOGICAL(all_int)[0] = 0L;
      break;
    }
  }

  UNPROTECT(2);
  return all_int;

}

// helper function determining the lengths of all elements of a list.

SEXP R_lengths(SEXP X){
  PROTECT(X);
  int n = length(X);
  SEXP out;
  out = PROTECT(allocVector(INTSXP, n));
  
  int *outp = INTEGER(out);  

  for ( int i=0; i<n; i++, outp++ ) {
    (*outp) = length(VECTOR_ELT(X,i));
  }
  UNPROTECT(2);
  return out;
}















