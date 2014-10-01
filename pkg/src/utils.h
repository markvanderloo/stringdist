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

#ifndef sd_utils_h
#define sd_utils_h

/* integer recycling macro  */
#define RECYCLE(X,Y) ( (X) == (Y) ? 0 : (X) )

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define ABS(X) ((X)<0 ? -1*(X) : (X))

unsigned int max_length(SEXP);

/* Get element from SEXP list and determine some parameters.
 *
 * Input:
 * x: A list of integer vectors or a character vector
 * i: Index in x: what element to extract.
 * bytes: (boolean) if (bytes) then x is assumed to be a character vector.
 * 
 * Output:
 * len  : the length of the i'th object in x.
 * isna : wether the i'th object represents an NA
 * c    : if (bytes) then c will contain the values of the i'th element in x, coerced to integers.
 *
 * Return value: 
 * A pointer to the integer representation of the i'th object in x.
 *
 */
unsigned int *get_elem(SEXP x, int i, int bytes, int *len, int *isna, unsigned int *c);


/* (mutlithreaded) recycling.
 *
 * This avoids having to compute i % ni at every iteration while
 * recycling over a vector.
 *
 * Input:
 * i : integer, current index.
 * nthreads : number of threads working on the vector.
 * ni : vector length
 *
 *
 */
static inline int recycle(int i, int nthreads, int ni){
  i += nthreads;
  if ( i >= ni )
    i = (nthreads < ni) ? (i - ni) : (i % ni);
  return i;
}




#endif
