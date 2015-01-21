
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

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "dist.h"
#include "stringdist.h"

#define MAX(X,Y) ((X) < (Y) ? (X) : (Y))


/* 
 *
 *
 *
 */
Stringdist *open_stringdist(Distance d, int str_len_a, int str_len_b, ...){
  va_list args;
  va_start(args, str_len_b);

  Stringdist *S = (Stringdist *) malloc(sizeof(Stringdist)); 
  (*S) = (Stringdist) {d, NULL, NULL, NULL, NULL, 0L, 0.0};
  
  switch (d){
    case osa :
      S->work = (double *) malloc( (str_len_a + 1) * (str_len_b + 1) * sizeof(double)); 
      S->weight = (double *) malloc(4*sizeof(double));
      memcpy(S->weight, va_arg(args, double *), 4*sizeof(double));
    case lv :
      S->work = (double *) malloc( (str_len_a + 1) * (str_len_b + 1) *sizeof(double));
      S->weight = (double *) malloc(3 * sizeof(double));
      memcpy(S->weight, va_arg(args, double *), 3*sizeof(double));
    case dl :
      S->dict = new_dictionary( str_len_a + str_len_b + 1);
      S->work = (double *) malloc( (str_len_a + 3) * (str_len_b + 3) * sizeof(double)); 
      S->weight = (double *) malloc(4*sizeof(double));
      memcpy(S->weight, va_arg(args, double *), 4*sizeof(double));
    case hamming :
      break;
    case lcs :
      S->work = (double *) malloc( (str_len_a + 1) * (str_len_b + 1) *sizeof(double));
    case qgram :
      S->tree = new_qtree(va_arg(args, unsigned int), 2L); 
    case cosine :
      S->tree = new_qtree(va_arg(args, unsigned int), 2L); 
    case jaccard :
      S->tree = new_qtree(va_arg(args, unsigned int), 2L); 
    case jw :
      S->work = (double *) malloc( sizeof(double) * MAX(str_len_a,str_len_b));
      S->weight = (double *) malloc(3*sizeof(double));
      memcpy(S->weight, va_arg(args, double *), 3*sizeof(double));
      S->p = va_arg(args, double);
    case soundex :
      break;
    default :
      break;
      // set errno, return NULL
  };

  va_end(args);
  return S; 
  
}

void close_stringdist(Stringdist *S){
  free(S->work);
  free(S->weight);
  free_dictionary(S->dict);
  free_qtree(S->tree);
  free(S);
}



double stringdist(Stringdist *S, unsigned int *str_a, int len_a, unsigned int *str_b, int len_b){
  double d = -1.0;
  unsigned int ifail;
  switch(S->distance){
    case osa :
      d = osa_dist(str_a, len_a, str_b, len_b, S->weight, S->work);
    case lv :
      d = lv_dist( str_a, len_a, str_b, len_b, S->weight, S->work);
    case dl :
      d = dl_dist(str_a, len_a, str_b, len_b, S->weight, S->dict, S->work);
    case hamming :
      d = hamming_dist(str_a, len_a, str_b, len_b);
    case lcs :
      d = lcs_dist(str_a, len_a, str_b, len_b, S->work); 
    case qgram :
      d = qgram_dist(str_a, len_a, str_b, len_b, S->q, S->tree, 0L);
    case cosine :
      d = qgram_dist(str_a, len_a, str_b, len_b, S->q, S->tree, 1L);
    case jaccard :
      d = qgram_dist(str_a, len_a, str_b, len_b, S->q, S->tree, 2L);
    case jw :
      d = jaro_winkler_dist(str_a, len_a, str_b, len_b, S->p, S->weight, S->work);
    case soundex :
      d = soundex_dist(str_a,len_a,str_b, len_b, &ifail);
      break;
    default :
      break;
      // set errno, return -1
  }

  return d;
}






