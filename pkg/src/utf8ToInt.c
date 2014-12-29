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

//#define USE_RINTERNALS
#include <R.h>
#include <Rdefines.h>
#include "utils.h"
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif


/* This function is gratefully copied from the R core distribution.
 * It is replicated here to facilitate multicore processing through
 * openmp (as it is not part of R's API).
 *
 * Convert a single character to integer.
 *
 * *s input buffer (must be utf8)
 * *w output buffer
 *
 * value:
 * >0 : The number of bytes in the multi-byte representation
 * 0  : End of string reached.
 * <0 : Invalid input (not interpretable as UTF-8)
 * 
 * 
 */
static int mbrtoint(unsigned int *w, const char *s)
{
    unsigned int byte;
    byte = *((unsigned char *)s);

    if (byte == 0) {
  *w = 0;
  return 0;
    } else if (byte < 0xC0) {
  *w = (int) byte;
  return 1;
    } else if (byte < 0xE0) {
  if (!s[1]) return -2;
  if ((s[1] & 0xC0) == 0x80) {
      *w = (int) (((byte & 0x1F) << 6) | (s[1] & 0x3F));
      return 2;
  } else return -1;
    } else if (byte < 0xF0) {
  if (!s[1] || !s[2]) return -2;
  if (((s[1] & 0xC0) == 0x80) && ((s[2] & 0xC0) == 0x80)) {
      *w = (int) (((byte & 0x0F) << 12)
      | ((s[1] & 0x3F) << 6) | (s[2] & 0x3F));
      byte = *w;
      if (byte >= 0xD800 && byte <= 0xDFFF) return -1; /* surrogate */
      if (byte == 0xFFFE || byte == 0xFFFF) return -1;
      return 3;
  } else return -1;
    } else if (byte < 0xF8) {
  if (!s[1] || !s[2] || !s[3]) return -2;
  if (((s[1] & 0xC0) == 0x80)
      && ((s[2] & 0xC0) == 0x80)
      && ((s[3] & 0xC0) == 0x80)) {
      *w = (int) (((byte & 0x07) << 18)
      | ((s[1] & 0x3F) << 12)
      | ((s[2] & 0x3F) << 6)
      | (s[3] & 0x3F));
      byte = *w;
      return 4;
  } else return -1;
    } else if (byte < 0xFC) {
  if (!s[1] || !s[2] || !s[3] || !s[4]) return -2;
  if (((s[1] & 0xC0) == 0x80)
      && ((s[2] & 0xC0) == 0x80)
      && ((s[3] & 0xC0) == 0x80)
      && ((s[4] & 0xC0) == 0x80)) {
      *w = (int) (((byte & 0x03) << 24)
      | ((s[1] & 0x3F) << 18)
      | ((s[2] & 0x3F) << 12)
      | ((s[3] & 0x3F) << 6)
      | (s[4] & 0x3F));
      byte = *w;
      return 5;
  } else return -1;
    } else {
  if (!s[1] || !s[2] || !s[3] || !s[4] || !s[5]) return -2;
  if (((s[1] & 0xC0) == 0x80)
      && ((s[2] & 0xC0) == 0x80)
      && ((s[3] & 0xC0) == 0x80)
      && ((s[4] & 0xC0) == 0x80)
      && ((s[5] & 0xC0) == 0x80)) {
      *w = (int) (((byte & 0x01) << 30)
      | ((s[1] & 0x3F) << 24)
      | ((s[2] & 0x3F) << 18)
      | ((s[3] & 0x3F) << 12)
      | ((s[5] & 0x3F) << 6)
      | (s[5] & 0x3F));
      byte = *w;
      return 6;
  } else return -1;
    }
    /* return -2; not reached */
}


/* Translate a UTF-8 string to integers.
 *
 * Input
 *
 * str   : pointer to input string.
 * outbuf: pointer to output buffer.
 *
 * Returns:
 * The number of logical characters converted.
 *
 */
static int utf8_to_int(const char *str, unsigned int *outbuf){

  unsigned int *p = outbuf;
  char *s = (char *) str;
  int nbytes
    , str_len = 0L;

  for (int i=0;; i++){
    nbytes = mbrtoint(p, s);
    if (nbytes > 0){
      p += 1L;
      str_len += 1L;
      s += nbytes;
    } else if (nbytes == -1L) { // non-utf-8 sequence encountered.
      return -1; 
    } else if (nbytes == 0L ){
      return str_len;
    }
  }
}


/* todo: catch *len=-1 */
unsigned int *get_elem1(SEXP x, int i, int bytes, int *len, int *isna, unsigned int *c){

  *isna = ( STRING_ELT(x,i) == NA_STRING );
  if (bytes){
    (*len)  = length(STRING_ELT(x,i));
    for (int j=0; j < *len; j++ )
      c[j] =  CHAR(STRING_ELT(x,i))[j];
  } else {
    (*len)  = utf8_to_int( CHAR(STRING_ELT(x,i)), c);
  }
  if ( *len < 0 ){
    error("Encountered byte sequence not representing an utf-8 character.\n");
  }
  return  c;
}

/* byte-by-byte char to int translation 
 *
 *
 * Return value: the number of bytes converted.
 *
 *
 */
static int char_to_int(const char *str, unsigned int *outbuf){
  unsigned int *p = outbuf;
  char *s = (char *) str;
  int str_len = 0L;
  while (*s){
    (*p) =  (unsigned int) *s;
    p++;
    s++;
    str_len++;
  }
  return str_len;
}

Stringset *new_stringset(SEXP str, int bytes){
  size_t nstr = length(str);
  Stringset *s;
  s = (Stringset *) malloc(sizeof(Stringset));

  // get and set string lengths.
  s->str_len = (int *) malloc(nstr * sizeof(int));

  size_t nbytes = 0L;
  for (size_t i=0; i<nstr; i++){
    nbytes += length(STRING_ELT(str,i));
  }

  s->string = (unsigned int **) malloc(nstr * sizeof(int *));
  // room for int rep of strings, including a trailing zero (needed by e.g. by full dl-distance)
  // this is enough room for byte-by-byte translation, so for UTF-8 it will be too much.
  s->data = (unsigned int *) malloc( (nstr + nbytes) * sizeof(int));

  int *t = s->str_len;
  unsigned int *d = s->data;
  for (size_t i=0L; i < nstr; i++, t++){
    if ( STRING_ELT(str,i) == NA_STRING ){
      (*t) = NA_INTEGER; 
    } else {
      if (bytes){
        (*t) = char_to_int(CHAR(STRING_ELT(str,i)), d);
      } else {
        (*t) = utf8_to_int(CHAR(STRING_ELT(str,i)), d); 
      }
      s->string[i] = d;
      (*(d + (*t))) = 0L; // append a zero.
      d += (*t) + 1L;
    }
  }
  return s;
}

void free_stringset(Stringset *s){
  free(s->string);
  free(s->data);
  free(s->str_len);
  free(s);  
}




/*
SEXP stringset(SEXP str, SEXP useBytes){
  int bytes = INTEGER(useBytes)[0];

  Stringset *s = new_stringset(str, bytes);
  for ( size_t i=0; i < length(str); i++){
    Rprintf("string of length %d, ints:", s->str_len[i]);
    for ( int j=0; j < s->str_len[i]; j++){
      Rprintf("%d ", *(s->string[i] + j));
    }
    Rprintf("\n");
  }
  for (int i=0; i<10; ++i){
    Rprintf("%d, ",s->data[i]);
  }; Rprintf("\n");
 
  free_stringset(s);
  return R_NilValue;
}
*/
















