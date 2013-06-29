/* This code is gratefully based on Nick Logan's github repository
 * https://github.com/ugexe/Text--Levenshtein--Damerau--XS/blob/master/damerau-int.c
 * 
 *
 * Changes/additions wrt original code:
 * - Added R.h, Rdefines.h inclusion
 * - Added R interface function
 * - Added edit weights (function is now of type double)
 * - Added corner cases for length-zero strings.
 * - Replaced linked list dictionary with fixed-size struct for loop externalization.
 * - Externalized allocation of dynamic programming matrix.
 * 
 * mark.vanderloo@gmail.com
 */

/* ugexe@cpan.org (Nick Logan)    */

#define USE_RINTERNALS
#include <R.h>
#include <Rdefines.h>
#include "utils.h"


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

static double distance(
      unsigned int *src,
      unsigned int *tgt,
      unsigned int x,
      unsigned int y,
      double *weight,
      double maxDistance,
      dictionary *dict,
      double *scores
    ){

  if (x == 0){
    if ( maxDistance > 0 && maxDistance < y ){
      return -1;
    } else {
      return (double) y;
    }
  }
  if (y == 0){
    if (maxDistance > 0 && maxDistance < x){
      return -1;
    } else {
      return (double) x;
    }
  }

  unsigned int swapCount, targetCharCount,i,j;
  double delScore, insScore, subScore, swapScore;
  unsigned int score_ceil = x + y;
  double colmin;
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
    colmin = (double) x + y + 1;
    for(j=1;j<=y;j++){
      if(i == 1) {
        uniquePush(dict,tgt[j]);
        scores[1 * (y + 2) + (j + 1)] = j;
        scores[0 * (y + 2) + (j + 1)] = score_ceil;
      }
      targetCharCount = dict->value[which(dict, tgt[j-1])];
      swapScore = scores[targetCharCount * (y + 2) + swapCount] + i - targetCharCount - 2 + j - swapCount + weight[3];

      if(src[i-1] != tgt[j-1]){
        subScore = scores[i * (y + 2) + j] + weight[2];
        insScore = scores[(i+1) * (y + 2) + j] + weight[1];
        delScore = scores[i * (y + 2) + (j + 1)] + weight[0];
        scores[(i+1) * (y + 2) + (j + 1)] = min2(swapScore, min3(delScore, insScore, subScore));
      } else {
        swapCount = j;
        scores[(i+1) * (y + 2) + (j + 1)] = min2(scores[i * (y + 2) + j], swapScore);
      } 
      colmin = min2(colmin,scores[(i+1)*(y+2) + (j+1)]);
    }
    /* We will return -1 here if the */
    /* current minimum > maxDistance   */
    if(maxDistance > 0 && maxDistance < colmin) {

      reset_dictionary(dict);
      return -1;
    }

      dict->value[which(dict,src[i-1])] = i;    
  }

  {
  double score = scores[(x+1) * (y + 2) + (y + 1)];
  reset_dictionary(dict);
  return score;
  }
}

/* End of workhorse */

// -- interface with R 


SEXP R_dl(SEXP a, SEXP b, SEXP weight, SEXP maxDistance){
  PROTECT(a);
  PROTECT(b);
  PROTECT(maxDistance);
  PROTECT(weight);
   
  int i=0, j=0, k=0;
  int na = length(a);
  int nb = length(b);
  int nt = (na > nb) ? na : nb;
  double maxDist = REAL(maxDistance)[0];
  double *w = REAL(weight);

  SEXP yy; 
  PROTECT(yy = allocVector(REALSXP, nt));
  double *y = REAL(yy);

  /* claim space for workhorse */
  int max_a = max_length(a);
  int max_b = max_length(b);
  dictionary *dict = new_dictionary( max_length(a) + max_length(b) + 1 );
  double *scores = (double *) malloc( (max_a + 3) * (max_b + 2) * sizeof(double) );

  for ( k=0; k < nt; ++k ){
    if ( INTEGER(VECTOR_ELT(a,i))[0] == NA_INTEGER || INTEGER(VECTOR_ELT(b,j))[0] == NA_INTEGER){
      y[k] = NA_REAL;
      continue;
    }
    y[k] = distance(
     (unsigned int *) INTEGER(VECTOR_ELT(a,i)),
     (unsigned int *) INTEGER(VECTOR_ELT(b,j)),
      length(VECTOR_ELT(a,i)),
      length(VECTOR_ELT(b,j)),
      w,
      maxDist,
      dict,
      scores
    );
    i = RECYCLE(i+1,na);
    j = RECYCLE(j+1,nb);
  }
  free_dictionary(dict);
  free(scores);
  UNPROTECT(5);
  return yy;
} 

