
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

#include <R_ext/Rdynload.h>

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_hidden
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @mainpage Stringdist C API
 *
 * @author Mark van der Loo, Jan van der Laan, R Core Team, Paul Hsieh, Chris Muir
 * @version `R package stringdist` version `0.9.5.0` and higher.
 *
 * @section using Using the stringdist C API
 * To call the functions described here from your package you need to:
 * 
 * 1. Make sure that `stringdist` is installed.
 * 2. Add `stringdist` to `Imports` (or `Depends`) and `LinkingTo` in the `DESCRIPTION` file.
 * 3. In your source file under the package's `/src` directory, add the line 
 * ```
 * #include <stringdist_api.h>
 * ```
 * 
 * 
 * An example of a published package using this API is
 * [refinr](https://CRAN.R-project.org/package=refinr). A minimal example can be
 * found [here](https://github.com/markvanderloo/linkstringdist).
 *
 * @section encoding Character encoding
 * All `character` vector input is expected to be in `UTF-8` (this also allows
 * `ASCII`).  Distance computations are based on UTF [code
 * points](https://en.wikipedia.org/wiki/Code_point) unless `useBytes` is
 * `TRUE`, in which case distances are computed over byte sequences. Using
 * non-UTF-8 encoded strings is untested and is highly likely to result in
 * errors.
 *
 * @section threads Thread safety
 *
 * It is not safe to call functions from `stringdist` C API from 
 * multiple concurrent threads.
 *
 *
 *
 */


/** 
 * @file stringdist_api.h
 * @brief Functions exported from the stringdist package.
 *
 */


/*
SEXP attribute_hidden sd_all_int(SEXP X)
{
  static SEXP(*fun)(SEXP) = NULL;
  if (fun == NULL) fun = (SEXP(*)(SEXP)) R_GetCCallable("stringdist","R_all_int");
  return fun(X);
}
*/

/**
 * @brief Find the location of values in `x` in `table` by approximate matching.
 * 
 * @param x `[character]` vector.
 * @param table `[character]` vector (lookup table)
 * @param method `[integer]` scalar, indicating the distance method as follows
 *   @parblock
 *    - 0: Optimal String Alignment (`"osa"`)
 *    - 1: Levenshtein (`"lv"`)
 *    - 2: Damerau-Levenshtein (`"dl"`)
 *    - 3: Hamming (`"hamming"`)
 *    - 4: Longest Common Substring (`"lcs"`)
 *    - 5: q-gram (`"qgram"`)
 *    - 6: cosine (`"cosine"`)
 *    - 7: Jaccard (`"jaccard"`)
 *    - 8: Jaro-Winkler (`"jw"`)
 *    - 9: Soundex (`"soundex"`)
 *   @endparblock
 * @param nomatch `[integer]` The value to be returned when no match is found. 
 * @param matchNA Should `NA`s be matched? Default behaviour mimics the
 *   behaviour of base `match`, meaning that `NA` matches `NA` (see also the note
 *   on `NA` handling below).
 * @param weight `[numeric]` vector. Edit penalty
 *   @parblock
 *     For `method='osa'` or`'dl'`, the penalty for
 *     deletion, insertion, substitution and transposition, in that order. When
 *     `method='lv'`, the penalty for transposition is ignored. When
 *     `method='jw'`, the weights associated with characters of `a`,
 *     characters from `b` and the transposition weight, in that order. 
 *     Weights must be positive and not exceed 1. `weight` is ignored
 *     completely for other methods
 *   @endparblock
 * @param q  `[integer]` scalar. Size of the q-gram; must be nonnegative. Only
 *   applies to `method='qgram'`, `'jaccard'` or `'cosine'`.
 * @param maxDistance `[numeric]` scalar. The maximum distance allowed for matching.
 * @param p `[numeric]` scalar. Penalty factor for Jaro-Winkler distance. The
 *   valid range for `p` is `0 <= p <= 0.25`. If `p=0` (default), the
 *   Jaro-distance is returned. Applies only to `method='jw'`.
 * @param bt `[numeric]` vector. Winkler's boost threshold. Winkler's penalty
 *   factor is only applied when the Jaro distance is larger than `bt`.  Applies
 *   only to `method='jw'` and `p>0`.
 * @param useBytes Perform byte-wise comparison (i.e. do not translate UTF-8 to
 *   integer prior to distance calculation)
 * @param nthread `[integer]` scalar. Maximum number of threads to use. 
 *
 *
 * @return
 * `[integer]` vector of `length(x)` with indices in `table`.
 */
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

/**
 * @brief Compute q-gram counts
 * 
 * @param a `[character]` vector 
 * @param qq `[integer`] scalar. 
 * 
 * @return
 * A `[numeric]` vector of `length(a)*n_qgrams`, where `n_qrams` is the number
 * of different `qgrams` observed in the elements of `a`. The output vector has
 * an attribute called `qgrams`, which is an integer vector of size
 * `q*n_qgrams` containing integer (UTF-32) labels for the q-grams
 * sequentially.
 * 
 */
SEXP attribute_hidden sd_get_qgrams(SEXP a, SEXP qq)
{
  static SEXP(*fun)(SEXP, SEXP) = NULL;
  if (fun == NULL) fun = (SEXP(*)(SEXP, SEXP)) R_GetCCallable("stringdist","R_get_qgrams");
  return fun(a, qq);
}

/*
SEXP attribute_hidden sd_lengths(SEXP X)
{
  static SEXP(*fun)(SEXP) = NULL;
  if (fun == NULL) fun = (SEXP(*)(SEXP)) R_GetCCallable("stringdist","R_lengths");
  return fun(X);
}
*/

/**
 * @brief Lower tridiagonal elements of distance matrix.
 *
 * @param a `[character]` vector
 * @param method `[integer]` scalar, indicating the distance method as follows
 *   @parblock
 *    - 0: Optimal String Alignment (`"osa"`)
 *    - 1: Levenshtein (`"lv"`)
 *    - 2: Damerau-Levenshtein (`"dl"`)
 *    - 3: Hamming (`"hamming"`)
 *    - 4: Longest Common Substring (`"lcs"`)
 *    - 5: q-gram (`"qgram"`)
 *    - 6: cosine (`"cosine"`)
 *    - 7: Jaccard (`"jaccard"`)
 *    - 8: Jaro-Winkler (`"jw"`)
 *    - 9: Soundex (`"soundex"`)
 *   @endparblock
 * @param weight `[numeric]` vector. Edit penalty
 *   @parblock
 *     For `method='osa'` or`'dl'`, the penalty for deletion, insertion,
 *     substitution and transposition, in that order. When `method='lv'`, the
 *     penalty for transposition is ignored. When `method='jw'`, the weights
 *     associated with characters of `a`, characters from `b` and the
 *     transposition weight, in that order.  Weights must be positive and not
 *     exceed 1. `weight` is ignored completely for other methods
 *   @endparblock
 * @param q  `[integer]` scalar. Size of the q-gram; must be nonnegative. Only
 *   applies to `method='qgram'`, `'jaccard'` or `'cosine'`.
 * @param p `[numeric]` scalar. Penalty factor for Jaro-Winkler distance. The
 *   valid range for `p` is `0 <= p <= 0.25`. If `p=0` (default), the
 *   Jaro-distance is returned. Applies only to `method='jw'`.
 * @param bt `[numeric]` vector. Winkler's boost threshold. Winkler's penalty
 *   factor is only applied when the Jaro distance is larger than `bt`.  Applies
 *   only to `method='jw'` and `p>0`.
 * @param useBytes Perform byte-wise comparison (i.e. do not translate UTF-8 to
 *   integer prior to distance calculation)
 * @param nthread `[integer]` scalar. Maximum number of threads to use. 
 *
 * @return
 *   A `[numeric]` vector of length `n*(n-1)/2`, where `n=length(a)`. It contains
 *   the positive values of consequtive columns of the distance matrix. Also see
 *   the R-code in  `stringdist:::lower_tri`.
 */
SEXP attribute_hidden sd_lower_tri(SEXP a, SEXP method
                                     , SEXP weight, SEXP p,  SEXP bt, SEXP q
                                     , SEXP useBytes, SEXP nthrd)
{
  static SEXP(*fun)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = NULL;
  if (fun == NULL) fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("stringdist","R_lower_tri");
  return fun(a, method, weight, p, bt, q, useBytes, nthrd);
}

/**
 * @brief Compute soundex code
 *
 * @param[in] x `[character]` vector
 * @param[in] useBytes `[logical]` scalar.
 *
 * @return
 *
 * A `list` with `length(x)` element. Each element is a length 4 integer vector 
 * representing a 4-character soundex code. The integers are ASCII code points.
 *
 */
SEXP attribute_hidden sd_soundex(SEXP x, SEXP useBytes)
{
  static SEXP(*fun)(SEXP, SEXP) = NULL;
  if (fun == NULL) fun = (SEXP(*)(SEXP, SEXP)) R_GetCCallable("stringdist","R_soundex");
  return fun(x, useBytes);
}

/**
 * @brief compute string distances
 * 
 * @param a `[character]` vector
 * @param b `[character]` vector
 * @param method `[integer]` scalar, indicating the distance method as follows
 *   @parblock
 *    - 0: Optimal String Alignment (`"osa"`)
 *    - 1: Levenshtein (`"lv"`)
 *    - 2: Damerau-Levenshtein (`"dl"`)
 *    - 3: Hamming (`"hamming"`)
 *    - 4: Longest Common Substring (`"lcs"`)
 *    - 5: q-gram (`"qgram"`)
 *    - 6: cosine (`"cosine"`)
 *    - 7: Jaccard (`"jaccard"`)
 *    - 8: Jaro-Winkler (`"jw"`)
 *    - 9: Soundex (`"soundex"`)
 *   @endparblock
 * @param weight `[numeric]` vector. Edit penalty
 *   @parblock
 *     For `method='osa'` or`'dl'`, the penalty for deletion, insertion,
 *     substitution and transposition, in that order. When `method='lv'`, the
 *     penalty for transposition is ignored. When `method='jw'`, the weights
 *     associated with characters of `a`, characters from `b` and the
 *     transposition weight, in that order.  Weights must be positive and not
 *     exceed 1. `weight` is ignored completely for other methods
 *   @endparblock
 * @param q  `[integer]` scalar. Size of the q-gram; must be nonnegative. Only
 *   applies to `method='qgram'`, `'jaccard'` or `'cosine'`.
 * @param p `[numeric]` scalar. Penalty factor for Jaro-Winkler distance. The
 *   valid range for `p` is `0 <= p <= 0.25`. If `p=0` (default), the
 *   Jaro-distance is returned. Applies only to `method='jw'`.
 * @param bt `[numeric]` vector. Winkler's boost threshold. Winkler's penalty
 *   factor is only applied when the Jaro distance is larger than `bt`.  Applies
 *   only to `method='jw'` and `p>0`.
 * @param useBytes Perform byte-wise comparison (i.e. do not translate UTF-8 to
 *   integer prior to distance calculation)
 * @param nthread `[integer]` scalar. Maximum number of threads to use. 
 *
 *
 * @return
 *   A `[numeric]` vector of length `max(length(a),length(b))` where the shortest
 *   vector is recycled over the longer (no warnings are given when the longer
 *   length is not an integer multiple of the shorter length).
 *
 *
 */
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
