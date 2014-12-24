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
#ifdef _OPENMP
#include <omp.h>
#endif

/* Optimal string alignment algorithm. 
 * Computes Damerau-Levenshtein distance, restricted to single transpositions.
 * - See pseudocode at http://en.wikipedia.org/wiki/Damerau%E2%80%93Levenshtein_distance
 * - Extended with custom weights and maxDistance
 */
static double osa(unsigned int *a, int na, unsigned int *b, int nb, double *weight, double *scores){

  if (!na){
    return (double) nb;
  }
  if (!nb){
    return (double) na;
  }

  int i, j;
  int M, I = na+1, L=na+1, J = nb+1;
  double sub, tran;

   for ( i = 0; i < I; ++i ){
      scores[i] = i;
   }
   for ( j = 1; j < J; ++j, L += I ){
      scores[L] = j;
   }

for ( i = 1; i <= na; ++i ){
  L = I; M = 0;
  for ( j = 1; j <= nb; ++j, L += I, M += I ){
    if (a[i-1] == b[j-1]){
      sub = 0;
      tran= 0;
    } else {
      sub = weight[2];
      tran= weight[3];
    }

    scores[i + L] = MIN(MIN( 
      scores[i-1 + L] + weight[0],     // deletion
      scores[i   + M] + weight[1]),    // insertion
      scores[i-1 + M] + sub            // substitution
    );
    if ( i>1 && j>1 && a[i-1] == b[j-2] && a[i-2] == b[j-1] ){
      scores[i + L] = MIN(scores[i + L], scores[i-2 + M-I] + tran); // transposition
    }
  }
}
  double score = scores[I*J-1];
  return score;
}

//-- Distance function interface with R


SEXP R_osa(SEXP a, SEXP b, SEXP weight, SEXP useBytes, SEXP nthrd){
  PROTECT(a);
  PROTECT(b);
  PROTECT(weight);
  PROTECT(useBytes);
  PROTECT(nthrd);

  int na = length(a)
    , nb = length(b)
    , bytes = INTEGER(useBytes)[0]
    , ml_a = max_length(a)
    , ml_b = max_length(b)
    , nt = (na > nb) ? na : nb;
  
  double *w = REAL(weight);

  // output vector
  SEXP yy;
  PROTECT(yy = allocVector(REALSXP, nt));
  double *y = REAL(yy);


  #ifdef _OPENMP 
  int  nthreads = INTEGER(nthrd)[0];
  #pragma omp parallel num_threads(nthreads) default(none) \
      shared(y, w, R_PosInf, NA_REAL, bytes, na, nb, ml_a, ml_b, nt, a, b)
  #endif
  {
    double *scores = (double *) malloc( (ml_a + 1) * (ml_b + 1) * sizeof(double)); 

    unsigned int *s = NULL, *t = NULL;
    s = (unsigned int *) malloc(( 2L + ml_a + ml_b) * sizeof(int));

    if ( (scores == NULL) | (bytes && s == NULL) ) nt = -1;
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
        y[k] = osa(
           s, len_s 
         , t, len_t
         , w, scores
        );
        if ( y[k] < 0 ) y[k] = R_PosInf;
      }
      i = recycle(i, num_threads, na);
      j = recycle(j, num_threads, nb);
    }
     
    free(scores);
    if (bytes) free(s);
  } // end of parallel region

  UNPROTECT(6);
  if (nt < 0 ) error("Unable to allocate enough memory");
  return(yy);
}



//-- Match function interface with R

SEXP R_match_osa(SEXP x, SEXP table, SEXP nomatch, SEXP matchNA, SEXP weight
  , SEXP maxDistance, SEXP useBytes, SEXP nthrd){
  PROTECT(x);
  PROTECT(table);
  PROTECT(nomatch);
  PROTECT(matchNA);
  PROTECT(weight);
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

  double *w = REAL(weight);
  double maxDist = REAL(maxDistance)[0];
  // output vector
  SEXP yy;
  PROTECT(yy = allocVector(INTSXP, nx));
  int *y = INTEGER(yy);
  
  #ifdef _OPENMP
  int nthreads = INTEGER(nthrd)[0];
  #pragma omp parallel num_threads(nthreads) default(none) \
    shared(x,table, y, R_PosInf, nx, ntable, no_match, match_na, bytes, ml_x, ml_t, w, maxDist)
  #endif
  {
    /* claim space for workhorse */
    double *work = (double *) malloc( (ml_x + 3) * (ml_t + 2) * sizeof(double) );

    unsigned int *X = NULL, *T = NULL;
    X = (unsigned int *) malloc( (2L + ml_x + ml_t) * sizeof(int));
    T = X + ml_x + 1L;

    if ( (work == NULL) | (bytes && X == NULL) ) nx = -1;

    double d = R_PosInf, d1 = R_PosInf;
    int index, isna_X, isna_T, len_X, len_T;

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

          d = osa(
            X, len_X, T, len_T, w, work
          );
          if ( d <= maxDist && d < d1){ 
            index = j + 1;
            if ( d == 0.0 ) break;
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
  UNPROTECT(9);
  if (nx < 0 ) error("Unable to allocate enough memory");
  return(yy);
}


