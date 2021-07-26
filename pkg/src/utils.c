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

  if (TYPEOF(x) == VECSXP){
    for (int i=0; i<length(x); ++i){
      m = length(VECTOR_ELT(x,i));
      if (t < m) t = m;
    }
  } else {
    for (int i=0; i<length(x); ++i){
      m = length(STRING_ELT(x,i));
      if (t < m) t = m;
    }

  }
  return t;
}


