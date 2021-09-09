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

/* Optimal string alignment algorithm. 
 * Computes Damerau-Levenshtein distance, restricted to single transpositions.
 * - See pseudocode at http://en.wikipedia.org/wiki/Damerau%E2%80%93Levenshtein_distance
 * - Extended with custom weights and maxDistance
 */
double osa_dist(unsigned int *a, int na, unsigned int *b, int nb, double *weight, double *scores){

  if (!na){
    return (double) nb * weight[1]; // ins weight
  }
  if (!nb){
    return (double) na * weight[0]; // del weight
  }

  int i, j;
  int M, I = na+1, L=na+1, J = nb+1;
  double sub, tran;

   for ( i = 0; i < I; ++i ){
      scores[i] = i * weight[0];
   }
   for ( j = 1; j < J; ++j, L += I ){
      scores[L] = j * weight[1];
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

