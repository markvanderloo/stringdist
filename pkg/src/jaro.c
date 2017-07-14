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

#include "utils.h"
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif


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
 * bt   : Winkler's (undocumented?) boost threshold
 * work : workspace, minimally of length x + y
 *
 */
double jaro_winkler_dist(
     unsigned int *a 
     , int x
     , unsigned int *b
     , int y
     , double p
     , double bt
     , double *w
     , double *work
  ){


  // edge case
  if ( x == 0 && y == 0 ) return 0;

  for (int k=0; k < x + y; k++) work[k] = 0;

  // we need space for integers (or do a lot of conversions)
  unsigned int *wrk = (unsigned int*) work;
  unsigned int *matcha = wrk
             , *matchb = wrk + x;  
  unsigned int left, right;

  // number of matches
  int m = 0;
  // max transposition distance
  int M = MAX(MAX(x,y)/2 - 1,0);

  // store the match indices. Indices are stored as i+1 because 0 is used as 'no match'
  for ( int i = 0; i < x; ++i){
    left = MAX(0,i-M);
    right = MIN(y,i+M);
    for ( int j = left; j <= right; j++){
      if ((a[i] == b[j]) && (matchb[j]==0)){
         matcha[i] = i+1;
         matchb[j] = j+1;
         m += 1;
         break;
      }
    }
  }

  // copy matches so they're easy to compare for transposition counting
  int j = 0;
  for (int i=0; i < x; ++i){
    if (matcha[i]){ 
      matcha[j] = a[matcha[i]-1];
      ++j;
    }
  }
  j = 0;
  for (int i=0; i < y; ++i){
    if (matchb[i]){ 
      matchb[j] = b[matchb[i]-1];
      ++j;
    }
  }

  // count 'transpositions', the Jaro way.
  double t = 0.0;
  for ( int k=0; k<m; ++k){
    if (matcha[k] != matchb[k]) t += 0.5;
  }

  double d;
  if ( m < 1 ){
    d = 1.0;
  } else {
    d = 1.0 - (1.0/3.0)*(w[0]*m/((double) x) + w[1]*m/((double) y) + w[2]*(m-t)/m);
  }
  // Winkler's penalty factor
  if ( p > 0 && d > bt ){
    int n = MIN(MIN(x,y),4);
    d =  d - get_l(a,b,n)*p*d; 
  }

  return d;
}





