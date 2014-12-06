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
#ifdef _OPENMP
#include <omp.h>
#endif
#include "utils.h"

/* Longest common substring
 * - basically edit distance, only allowing insertions and deletions, at the cost of 1.
 */
static int lcs(unsigned int *a, int na, unsigned int *b, int nb, int *scores){
  if (!na){
    return (double) nb;
  }
  if (!nb){
    return (double) na;
  }

  int i, j;
  int M, I = na+1, L = na+1, J = nb+1;
  
  for ( i = 0; i < I; ++i ){
    scores[i] = i;
  }
  for ( j = 1; j < J; ++j, L += I ){
    scores[L] = j;
  }

  for ( i = 1; i <= na; ++i ){
    M = 0; L = I;
    for ( j = 1; j <= nb; ++j, L += I, M += I ){
      if ( a[i-1] == b[j-1] ){ // equality, copy previous score
        scores[i + L] = scores[i-1 + M];
      } else {
        scores[i + L] = MIN(
          scores[i-1 + L] + 1 ,     // deletion
          scores[i   + M] + 1       // insertion
        );
      }
      
    }
  }
  int score = scores[I*J - 1];
  return score;
}

//-- interface with R


SEXP R_lcs(SEXP a, SEXP b, SEXP useBytes, SEXP nthrd){
  PROTECT(a);
  PROTECT(b);
  PROTECT(useBytes);
  PROTECT(nthrd);

  int na = length(a)
    , nb = length(b)
    , ml_a = max_length(a)
    , ml_b = max_length(b)
    , bytes = INTEGER(useBytes)[0]
    , nt = (na > nb) ? na : nb;

  // output vector
  SEXP yy;
  PROTECT(yy = allocVector(REALSXP, nt));
  double *y = REAL(yy);   
  
  #ifdef _OPENMP 
  int  nthreads = INTEGER(nthrd)[0];
  #pragma omp parallel num_threads(nthreads) default(none) \
      shared(y, R_PosInf, NA_REAL, bytes, na, nb, ml_a, ml_b, nt, a, b)
  #endif
  {
    // space for the workfunction
    int *scores; 
    scores = (int *) malloc( (ml_a + 1) * (ml_b + 1) * sizeof(int)); 

    unsigned int *s = NULL, *t = NULL;
    s = (unsigned int *) malloc( (2L + ml_a + ml_b) * sizeof(int));
    t = s + ml_a + 1L; 
    // no memory = no joy = no looping.
    if ( (scores == NULL) | (bytes && s == NULL) ) nt = -1;

    int len_s, len_t, isna_s, isna_t
      , i = 0, j = 0, ID = 0, num_threads=1;

    #ifdef _OPENMP
    ID = omp_get_thread_num();
    num_threads = omp_get_num_threads();
    i = recycle(ID-num_threads, num_threads, na);
    j = recycle(ID-num_threads, num_threads, nb);
    #endif
    
    for ( int k=ID; k < nt; k += num_threads ){
      get_elem1(a, i, bytes, &len_s, &isna_s, s);
      get_elem1(b, j, bytes, &len_t, &isna_t, t);
      if ( isna_s || isna_t ){
        y[k] = NA_REAL;
        continue;
      }
      y[k] = lcs(s, len_s, t, len_t, scores );
      if (y[k] < 0 ) y[k] = R_PosInf;
      i = recycle(i, num_threads, na);
      j = recycle(j, num_threads, nb);
    } 
    free(scores);
    free(s);
  } // end parallel region

  UNPROTECT(5);
  if (nt < 0)  error("Unable to allocate enough memory");
  return(yy);
}



//-- Match function interface with R

SEXP R_match_lcs(SEXP x, SEXP table, SEXP nomatch, SEXP matchNA
    , SEXP maxDistance, SEXP useBytes, SEXP nthrd){
  PROTECT(x);
  PROTECT(table);
  PROTECT(nomatch);
  PROTECT(matchNA);
  PROTECT(maxDistance);
  PROTECT(useBytes);
  PROTECT(nthrd);

  int nx = length(x)
    , ntable = length(table)
    , no_match = INTEGER(nomatch)[0]
    , match_na = INTEGER(matchNA)[0]
    , max_dist = INTEGER(maxDistance)[0]
    , bytes = INTEGER(useBytes)[0]
    , ml_x = max_length(x)
    , ml_t = max_length(table);


  // output vector
  SEXP yy;
  PROTECT(yy = allocVector(INTSXP, nx));
  int *y = INTEGER(yy);

  #ifdef _OPENMP
  int nthreads = INTEGER(nthrd)[0];
  #pragma omp parallel num_threads(nthreads) default(none) \
    shared(x,table, y, R_PosInf, nx, ntable, no_match, match_na, bytes, ml_x, ml_t, max_dist)
  #endif
  {
    // space for the workfunction
    int *work = (int *) malloc( (max_length(x) + 1) * (max_length(table) + 1) * sizeof(int)); 
    unsigned int *X = NULL, *T = NULL;

    X = (unsigned int *) malloc((2L + ml_x + ml_t)*sizeof(int));
    T = X + ml_x + 1L;
    if ( ( work == NULL) | (bytes && X == NULL) ) nx = -1;


    double d = R_PosInf, d1 = R_PosInf;

    int index, len_X, isna_X, len_T, isna_T;
    
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for ( int i=0; i<nx; i++){
      index = no_match;


      get_elem1(x, i, bytes, &len_X, &isna_X, X);
      d1 = R_PosInf;
      for ( int j=0; j<ntable; j++){

        get_elem1(table, j, bytes, &len_T, &isna_T, T);
        if ( !isna_X && !isna_T ){        // both are char (usual case)
          d = (double) lcs(
            X, len_X, T, len_T, work
          );
          if ( d <= max_dist && d < d1){ 
            index = j + 1;
            if ( d == 0.0  ) break;
            d1 = d;
          }
        } else if ( isna_X && isna_T ) {  // both are NA
          index = match_na ? j + 1 : no_match;
          break;
        }
      }
      
      y[i] = index;
    }
    free(X);
    free(work);
  } // end of parallel region
  UNPROTECT(8);
  if (nx < 0) error("Unable to allocate enough memory");
  return(yy);
}
