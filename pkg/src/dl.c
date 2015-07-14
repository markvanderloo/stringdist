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


#include "utils.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "dictionary.h"

/*
static void print_dict(dictionary *d){
  for ( int i=0; i<d->length; i++){
    Rprintf("d[%d] = %d; ", i, d->key[i]);
  }
  Rprintf("\n");
}
*/

static void reset_dictionary(dictionary *d){
  int nbytes = sizeof(unsigned int)*(d->length);
  memset(d->key  , 0, nbytes);
  memset(d->value, 0, nbytes);
}

dictionary *new_dictionary(unsigned int length){
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

void free_dictionary(dictionary *d){
  if ( d != NULL ){
    free(d->key);
    free(d->value);
    free(d);
  }
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


// note: src (tgt) will be indexed to their x + 1 (y+1).
double dl_dist(
      unsigned int *src,
      int x,
      unsigned int *tgt,
      int y,
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
    scores[(i+1) * (y + 2) + 1] = i * weight[0];
    scores[(i+1) * (y + 2) + 0] = score_ceil;
    swapCount = 0;
    
    for(j=1;j<=y;j++){
      if(i == 1) {
        uniquePush(dict,tgt[j]);
        scores[1 * (y + 2) + (j + 1)] = j * weight[0];
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


