
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

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>


#ifdef __cplusplus
extern "C" {
#endif
  
/*
 FUNCTIONS
 */
SEXP R_stringdist(SEXP a, SEXP b, SEXP method
                    , SEXP weight, SEXP p, SEXP bt, SEXP q
                    , SEXP useBytes, SEXP nthrd);
SEXP R_amatch(SEXP x, SEXP table, SEXP method 
                , SEXP nomatch, SEXP matchNA
                , SEXP weight, SEXP p, SEXP bt, SEXP q
                , SEXP maxDistance, SEXP useBytes
                , SEXP nthrd);
SEXP R_lower_tri(SEXP a, SEXP method
                   , SEXP weight, SEXP p,  SEXP bt, SEXP q
                   , SEXP useBytes, SEXP nthrd);
SEXP R_soundex(SEXP x, SEXP useBytes);
SEXP R_get_qgrams(SEXP a, SEXP qq);
SEXP R_lengths(SEXP X);
SEXP R_all_int(SEXP X);


//#endif
  
#ifdef __cplusplus
}
#endif