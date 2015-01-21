
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


#ifndef SD_STRINGDIST_H
#define SD_STRINGDIST_H

#include "dictionary.h"
#include "qtree.h"

typedef enum Distance { osa, lv, dl, hamming, lcs, qgram, cosine, jaccard, jw, soundex} Distance;
typedef struct {
  Distance distance;
  // workspace
  double *work;
  // [optional] weight vector
  double *weight;
  // dictionary object for dl-distance
  dictionary *dict;
  // tree object to store q-grams
  qtree *tree;
  // the q in qgrams
  unsigned int q;
  // Winkler's penalty factor
  double p;
} Stringdist;

Stringdist *open_stringdist(Distance, int, int, ...);

double stringdist(Stringdist *, unsigned int *, int, unsigned int *, int);

void close_stringdist(Stringdist *S);



#endif



