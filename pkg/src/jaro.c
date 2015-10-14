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



/* First match of a in b[] 
 * Returns -1 if no match is found
 * Parameter 'guard; indicates which elements of b have been matched before to avoid
 * matching two instances of the same character to the same position in b (which we treat read-only).
 */
static int match_int(unsigned int a, unsigned int *b, double *guard, int width, int m){
  int i = 0;
  while ( 
      ( i < width ) && 
      ( b[i] != a || (b[i] == a && guard[i] > 0)) 
  ){
    ++i;
  }
  // ugly edge case workaround
  if ( !(m && i==width) && b[i] == a ){
    guard[i] = 1.0;
    return i;
  } 
  return -1;
}

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
 * work : workspace, minimally of length max(x,y)
 *
 */
double jaro_winkler_dist(
             unsigned int *a, 
             int x,
             unsigned int *b,
             int y,
             double p,
             double *w,
             double *work
        ){

  // edge case
  if ( x == 0 && y == 0 ) return 0;
  // swap arguments if necessary, so we always loop over the shortest string
  if ( x > y ){
    unsigned int *c = b;
    unsigned int z = y;
    b = a;
    a = c;
    y = x;
    x = z;
  }

  for (int k=0; k<MAX(x,y); k++) work[k] = 0.0;

  // max transposition distance
  int M = MAX(MAX(x,y)/2 - 1,0);
  // transposition counter
  double t = 0.0;
  // number of matches 
  double m = 0.0;
  int max_reached; 
  int left, right, J, jmax=0;
  
  for ( int i=0; i < x; ++i ){
    left  = MAX(0, i-M);
    if ( left >= y ){
      J = -1;
    } else {
      right = MIN(y, i+M);
      // ugly workaround: I should rewrite match_int.
      max_reached = (right == y) ? 1 : 0;
      J =  match_int(a[i], b + left, work + left, right - left, max_reached);
    }

    if ( J >= 0 ){
      ++m;
      t += (J + left < jmax ) ? 1 : 0; 
      jmax = MAX(jmax, J + left);
    }
  }
  double d;
  if ( m < 1 ){
    d = 1.0;
  } else {
    double tot_weight = w[0] + w[1] + w[2];
    d = 1.0 - (1.0/tot_weight)*(w[0]*m/((double) x) + w[1]*m/((double) y) + w[2]*(m-t)/m);
  }

  // Winkler's penalty factor
  if ( p > 0 && d > 0 ){
    int n = MIN(MIN(x,y),4);
    d =  d - get_l(a,b,n)*p*d; 
  }
  return d;
}



