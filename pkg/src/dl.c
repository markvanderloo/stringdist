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
 *
 * Changes/additions wrt original code:
 * - Added R.h, Rdefines.h inclusion
 * - Added R interface function
 * - Added edit weights (function is now of type double)
 * - Added corner cases for length-zero strings.
 * - Replaced linked list dictionary with fixed-size struct for loop 
 *    externalization of memory allocation.
 * - Externalized allocation of dynamic programming matrix.
 * 
 * 
 */


#define USE_RINTERNALS
#include <R.h>
#include <Rdefines.h>
#include "utils.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/* Our unsorted dictionary  */
/* Note we use character ints, not chars. */


typedef struct {
  unsigned int *key;
  unsigned int *value;
  unsigned int length;
} dictionary;

static void reset_dictionary(dictionary *d){
  int nbytes = sizeof(unsigned int)*d->length;
  memset(d->key  , 0, nbytes);
  memset(d->value, 0, nbytes);
}

static dictionary *new_dictionary(unsigned int length){
  dictionary *d = (dictionary *) malloc(sizeof(dictionary));
  if ( d == NULL ){
    return NULL;
  }
  d->key   = (unsigned int *) malloc(length*sizeof(int));
  d->value = (unsigned int *) malloc(length*sizeof(int));
  if ( d->key == NULL || d->value == NULL){
    free(d->key);
    free(d->value);
    free(d);
    return NULL;
  }
  d->length = length;
  reset_dictionary(d);
  return d;
}

static void free_dictionary(dictionary *d){
  free(d->key);
  free(d->value);
  free(d);
}

static void uniquePush(dictionary *d, unsigned int key){
  int i=0;
  while (d->key[i] && d->key[i] != key){
    ++i;
  }
  d->key[i] = key;
}

static unsigned int which(dictionary *d, unsigned int key){
  int i=0;
  while( d->key[i] != key ){
     ++i;
  }
  return i;
}

/* End of Dictionary Stuff */


/* All calculations/work are done here */
// note: src (tgt) will be indexed to their x + 1 (y+1).
static double distance(
      unsigned int *src,
      unsigned int *tgt,
      unsigned int x,
      unsigned int y,
      double *weight,
      dictionary *dict,
      double *scores
    ){

  if (!x){
    return (double) y;
  }
  if (!y){
    return (double) x;
  }

  unsigned int swapCount, targetCharCount,i,j;
  double delScore, insScore, subScore, swapScore;
  unsigned int score_ceil = x + y;
  
  /* intialize matrix start values */
  scores[0] = score_ceil;  
  scores[1 * (y + 2) + 0] = score_ceil;
  scores[0 * (y + 2) + 1] = score_ceil;
  scores[1 * (y + 2) + 1] = 0;

  uniquePush(dict,src[0]);
  uniquePush(dict,tgt[0]);

  /* work loops    */
  /* i = src index */
  /* j = tgt index */
  for(i=1;i<=x;i++){ 
    uniquePush(dict,src[i]);
    scores[(i+1) * (y + 2) + 1] = i;
    scores[(i+1) * (y + 2) + 0] = score_ceil;
    swapCount = 0;
    
    for(j=1;j<=y;j++){
      if(i == 1) {
        uniquePush(dict,tgt[j]);
        scores[1 * (y + 2) + (j + 1)] = j;
        scores[0 * (y + 2) + (j + 1)] = score_ceil;
      }
      targetCharCount = dict->value[which(dict, tgt[j-1])];
      swapScore = scores[targetCharCount * (y + 2) + swapCount] + (i - targetCharCount - 1 + j - swapCount) *  weight[3];

      if(src[i-1] != tgt[j-1]){
        subScore = scores[i * (y + 2) + j] + weight[2];
        insScore = scores[(i+1) * (y + 2) + j] + weight[1];
        delScore = scores[i * (y + 2) + (j + 1)] + weight[0];
        scores[(i+1) * (y + 2) + (j + 1)] = MIN(swapScore, MIN(delScore, MIN(insScore, subScore) ));
      } else {
        swapCount = j;
        scores[(i+1) * (y + 2) + (j + 1)] = MIN(scores[i * (y + 2) + j], swapScore);
      }
    }
    
   dict->value[which(dict,src[i-1])] = i;    
  }

  double score = scores[(x+1) * (y + 2) + (y + 1)];
  reset_dictionary(dict);
  return score;
}

/* End of workhorse */

// -- interface with R 


SEXP R_dl(SEXP a, SEXP b, SEXP weight, SEXP nthrd){
  PROTECT(a);
  PROTECT(b);
  PROTECT(weight);
  PROTECT(nthrd);
   
  int na = length(a)
    , nb = length(b)
    , nt = (na > nb) ? na : nb
    , bytes = IS_CHARACTER(a)
    , ml_a = max_length(a)
    , ml_b = max_length(b);
  
  double *w = REAL(weight);

  // output
  SEXP yy; 
  PROTECT(yy = allocVector(REALSXP, nt));
  double *y = REAL(yy);



  #ifdef _OPENMP 
  int  nthreads = INTEGER(nthrd)[0];
  #pragma omp parallel num_threads(nthreads) default(none) \
      shared(y, w, R_PosInf, NA_REAL, bytes, na, nb, nt, ml_a, ml_b, a, b)
  #endif
  {
    /* claim space for workhorse */
    unsigned int *s=NULL, *t=NULL;
    dictionary *dict = new_dictionary( ml_a + ml_b + 1 );

    double *scores = (double *) malloc( (ml_a + 3) * (ml_b + 2) * sizeof(double) );

    int slen = (ml_a + ml_b + 2) * sizeof(int);
    s = (unsigned int *) malloc(slen);

    if ( (scores == NULL) | ( s == NULL ) ){
      nt = -1;
    } else {
      t = s + ml_a + 1;
      memset(s, 0, slen);
    }
    /* start working */
    
    unsigned int *s1, *t1;
    int k, len_s, len_t, isna_s, isna_t
      , i = 0, j = 0, ID = 0, num_threads = 1;

    #ifdef _OPENMP
    ID = omp_get_thread_num();
    num_threads = omp_get_num_threads();
    i = recycle(ID-num_threads, num_threads, na);
    j = recycle(ID-num_threads, num_threads, nb);
    #endif

    for ( k=ID; k < nt; k += num_threads ){
      if (bytes){
        s = get_elem(a, i, bytes, &len_s, &isna_s, s);
        t = get_elem(b, j, bytes, &len_t, &isna_t, t);
      } else { // make sure there's an extra 0 at the end of the string.
        s1 = get_elem(a, i, bytes, &len_s, &isna_s, s);
        t1 = get_elem(b, j, bytes, &len_t, &isna_t, t);
        memcpy(s,s1,len_s*sizeof(int));
        memcpy(t,t1,len_t*sizeof(int));
      }
      if ( isna_s || isna_t ){
        y[k] = NA_REAL;
        continue;
      }

      y[k] = distance(
       s, t, len_s, len_t,
       w, dict, scores
      );
      if (y[k] < 0 ) y[k] = R_PosInf;
      i = recycle(i, num_threads, na);
      j = recycle(j, num_threads, nb);
      memset(s, 0, slen);
    }
    
    free_dictionary(dict);
    free(scores);
    free(s);
  }
  UNPROTECT(5);
  if (nt < 0)  error("Unable to allocate enough memory");
 
  return yy;
} 

//-- Match function interface with R

SEXP R_match_dl(SEXP x, SEXP table, SEXP nomatch, SEXP matchNA
    , SEXP weight, SEXP maxDistance, SEXP nthrd){
  PROTECT(x);
  PROTECT(table);
  PROTECT(nomatch);
  PROTECT(matchNA);
  PROTECT(weight);
  PROTECT(maxDistance);
  PROTECT(nthrd);

  int nx = length(x)
    , ntable = length(table)
    , no_match = INTEGER(nomatch)[0]
    , match_na = INTEGER(matchNA)[0]
    , bytes = IS_CHARACTER(x)
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
    shared(x,table, y, w, R_PosInf, nx, ntable, no_match, match_na, bytes, ml_x, ml_t, maxDist)
  #endif
  {
    /* claim space for workhorse */
    dictionary *dict = new_dictionary( ml_x + ml_t + 1 );
    double *scores = (double *) malloc( (ml_x + 3) * (ml_t + 2) * sizeof(double) );

    unsigned int *X = NULL, *T = NULL;

    X = (unsigned int *) malloc( (ml_x + ml_t + 2) * sizeof(int) );

    if ( (scores == NULL) ||  (X == NULL) ) nx = -1;

    T = X + ml_x + 1;
    memset(X, 0, (ml_x + ml_t + 2)*sizeof(int));

    double d = R_PosInf, d1 = R_PosInf;
    int index, len_X, len_T, isna_X, isna_T;
    unsigned int *X1, *T1;
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for ( int i=0; i<nx; i++){
      index = no_match;
      if ( bytes ){
        X = get_elem(x, i , bytes, &len_X, &isna_X, X);
      } else {
        X1 = get_elem(x, i , bytes, &len_X, &isna_X, X);
        memcpy(X, X1, len_X*sizeof(int));
      }
      d1 = R_PosInf;

      for ( int j=0; j<ntable; j++){
        if ( bytes ){
          T = get_elem(table, j, bytes, &len_T, &isna_T, T);
        } else {
          T1 = get_elem(table, j, bytes, &len_T, &isna_T, T);
          memcpy(T, T1, len_T * sizeof(int));
        }
        if ( !isna_X && !isna_T ){        // both are char (usual case)
          d = distance(
            X, T, len_X, len_T, w, dict, scores
          );
          memset(T,0, (ml_t+1)*sizeof(int));
          if ( d <= maxDist && d < d1){ 
            index = j + 1;
            if ( fabs(d) < 1e-14 ) break;
            d1 = d;
          }
        } else if ( isna_X && isna_T ) {  // both are NA
          index = match_na ? j + 1 : no_match;
          break;
        }
      }
      
      y[i] = index;
      memset(X,0,(ml_x + 1)*sizeof(int));
    }  
    free(X);
    free_dictionary(dict);
    free(scores);
  } // end of parallel region
  UNPROTECT(8);
  if (nx < 0) error("Unable to allocate enough memory");
  return(yy);
}
