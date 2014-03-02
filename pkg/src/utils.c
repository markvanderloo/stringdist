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
 */


#define USE_RINTERNALS
#include <R.h>
#include <Rdefines.h>


unsigned int max_length(SEXP x){
  unsigned int t=0, m;
  for (int i=0; i<length(x); ++i){
    m = length(VECTOR_ELT(x,i));
    if (t < m) t = m;
  }
  return t;
}

unsigned int *get_elem(SEXP x, int i, int bytes, int *len, int *isna, unsigned int *c){
  unsigned int *out;
  if (bytes){
    *len  = length(STRING_ELT(x,i));
    *isna = ( STRING_ELT(x,i) == NA_STRING );
    for (int j=0; j < *len; j++ )
      c[j] =  CHAR(STRING_ELT(x,i))[j];
    out = c;
  } else {
    *len  = length(VECTOR_ELT(x,i));
    *isna = (INTEGER(VECTOR_ELT(x,i))[0] == NA_INTEGER);
    out = (unsigned int *) INTEGER(VECTOR_ELT(x,i));
  }
  return out;
}

