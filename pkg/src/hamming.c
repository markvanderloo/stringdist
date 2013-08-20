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
 *
 * This code is gratefully based on Nick Logan's github repository
 * https://github.com/ugexe/Text--Levenshtein--Damerau--XS/blob/master/damerau-int.c
 *
 */ 



#define USE_RINTERNALS
#include <R.h>
#include <Rdefines.h>
#include "utils.h"

int hamming(unsigned int *a, unsigned int *b, int n, int maxDistance){
   int h=0;
   for(int i=0; i<n; ++i){
      if (a[i] != b[i]) h++;
      if ( maxDistance > 0 && maxDistance < h ){
         return -1;
      }
   }
   return h;
}


// -- R interface

SEXP R_hm(SEXP a, SEXP b, SEXP maxDistance){
  PROTECT(a);
  PROTECT(b);
  PROTECT(maxDistance);

  int na = length(a)
    , nb = length(b)
    , nt = ( na > nb) ? na : nb
    , bytes = IS_CHARACTER(a)
    , maxDist = INTEGER(maxDistance)[0]
    , ml_a = max_length(a)
    , ml_b = max_length(b);

  unsigned int *s, *t;
  if ( bytes ){
    s = (unsigned int *) malloc( (ml_a + ml_b) * sizeof(int));
    if ( s == NULL ) error("Unable to allocate enough memory");
    t = s + ml_a;
  }

  SEXP yy;
  PROTECT(yy = allocVector(REALSXP,nt));
  double *y = REAL(yy);

  int i=0, j=0, k=0, len_s, len_t, isna_s, isna_t;
  for ( k=0; k<nt; 
        ++k
      , i = RECYCLE(i+1,na)
      , j = RECYCLE(j+1,nb) ){

    s = get_elem(a, i, bytes, &len_s, &isna_s, s);
    t = get_elem(b, j, bytes, &len_t, &isna_t, t);
    if ( isna_s || isna_t ){
      y[k] = NA_REAL;
      continue;         
    }
    if ( len_s != len_t ){
      y[k] = R_PosInf;
      continue;
    }
    y[k] = (double) hamming(s, t, len_s, maxDist);
    if (y[k] < 0) y[k] = R_PosInf;
  }

  if (bytes) free(s);
  UNPROTECT(4);
  return yy;
}


//-- Match function interface with R

SEXP R_match_hm(SEXP x, SEXP table, SEXP nomatch, SEXP matchNA, SEXP maxDistance){
  PROTECT(x);
  PROTECT(table);
  PROTECT(nomatch);
  PROTECT(matchNA);
  PROTECT(maxDistance);

  int nx = length(x), ntable = length(table);
  int no_match = INTEGER(nomatch)[0];
  int match_na = INTEGER(matchNA)[0];
  int max_dist = INTEGER(maxDistance)[0];


  // output vector
  SEXP yy;
  PROTECT(yy = allocVector(INTSXP, nx));
  int *y = INTEGER(yy);
  int *X, *T;


  double d = R_PosInf, d1 = R_PosInf;
  int nchar, index, xNA, tNA;

  for ( int i=0; i<nx; i++){
    index = no_match;
    nchar = length(VECTOR_ELT(x,i));

    X = INTEGER(VECTOR_ELT(x,i));
    xNA = (X[0] == NA_INTEGER);

    for ( int j=0; j<ntable; j++){
      if ( nchar != length(VECTOR_ELT(table,j)) ) continue;

      T = INTEGER(VECTOR_ELT(table,j));
      tNA = (T[0] == NA_INTEGER);

      if ( !xNA && !tNA ){        // both are char (usual case)
        d = (double) hamming(
          (unsigned int *) X,
          (unsigned int *) T,
          nchar,
          max_dist
        );
        if ( d > -1 && d < d1){ 
          index = j + 1;
          if ( d == 0.0 ) break;
          d1 = d;
        }
      } else if ( xNA && tNA ) {  // both are NA
        index = match_na ? j + 1 : no_match;
        break;
      }
    }
    
    y[i] = index;
  }  
  UNPROTECT(6);
  return(yy);
}

