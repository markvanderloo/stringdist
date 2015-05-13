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

#ifdef _OPENMP
#include <omp.h>
#endif
#include "utils.h"

/* Longest common substring
 * - basically edit distance, only allowing insertions and deletions, at the cost of 1.
 */
double lcs_dist(unsigned int *a, int na, unsigned int *b, int nb, double *scores){
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
   
  return scores[I*J - 1];
}

