
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

#ifndef _STRINGDIST_API_H
#define _STRINGDIST_API_H

#include <stringdist.h>		// also includes R.h, Rinternals.h, Rdefines.h

#include <Rconfig.h>
#include <R_ext/Rdynload.h>

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_hidden
#endif

#ifdef __cplusplus
extern "C" {
#endif


SEXP attribute_hidden sd_all_int(SEXP X)
{
  static SEXP(*fun)(SEXP) = NULL;
  if (fun == NULL) fun = (SEXP(*)(SEXP)) R_GetCCallable("stringdist","R_all_int");
  return fun(X);
}

SEXP attribute_hidden sd_amatch(SEXP x, SEXP table, SEXP method 
                                  , SEXP nomatch, SEXP matchNA
                                  , SEXP weight, SEXP p, SEXP bt, SEXP q
                                  , SEXP maxDistance, SEXP useBytes
                                  , SEXP nthrd)
{
  static SEXP(*fun)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = NULL;
  if (fun == NULL) fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("stringdist","R_amatch");
  return fun(x, table, method, nomatch, matchNA, weight, p, bt, q, maxDistance, useBytes, nthrd);
}

SEXP attribute_hidden sd_get_qgrams(SEXP a, SEXP qq)
{
  static SEXP(*fun)(SEXP, SEXP) = NULL;
  if (fun == NULL) fun = (SEXP(*)(SEXP, SEXP)) R_GetCCallable("stringdist","R_get_qgrams");
  return fun(a, qq);
}

SEXP attribute_hidden sd_lengths(SEXP X)
{
  static SEXP(*fun)(SEXP) = NULL;
  if (fun == NULL) fun = (SEXP(*)(SEXP)) R_GetCCallable("stringdist","R_lengths");
  return fun(X);
}

SEXP attribute_hidden sd_lower_tri(SEXP a, SEXP method
                                     , SEXP weight, SEXP p,  SEXP bt, SEXP q
                                     , SEXP useBytes, SEXP nthrd)
{
  static SEXP(*fun)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = NULL;
  if (fun == NULL) fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("stringdist","R_lower_tri");
  return fun(a, method, weight, p, bt, q, useBytes, nthrd);
}

SEXP attribute_hidden sd_soundex(SEXP x, SEXP useBytes)
{
  static SEXP(*fun)(SEXP, SEXP) = NULL;
  if (fun == NULL) fun = (SEXP(*)(SEXP, SEXP)) R_GetCCallable("stringdist","R_soundex");
  return fun(x, useBytes);
}

SEXP attribute_hidden sd_stringdist(SEXP a, SEXP b, SEXP method
                                      , SEXP weight, SEXP p, SEXP bt, SEXP q
                                      , SEXP useBytes, SEXP nthrd)
{
  static SEXP(*fun)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = NULL;
  if (fun == NULL) fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("stringdist","R_stringdist");
  return fun(a, b, method, weight, p, bt, q, useBytes, nthrd);
}


#ifdef __cplusplus
}
#endif

#endif