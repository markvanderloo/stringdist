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

/* This code is kindly contributed by Jan van der Laan (August 2014)
 */

#define USE_RINTERNALS
#include <R.h>
#include <Rdefines.h>
#include "utils.h"
#include <ctype.h>


// Translate similar sounding consonants to numeric codes; vowels are all 
// translated to 'a' and voiceless characters (and other characters) are 
// translated to 'h'.
// Upper and lower case ASCII characters are treated as separate cases,
// avoiding the use of 'tolower' whose effect depends on locale.
unsigned int translate_soundex(unsigned int c) {
  switch ( c ) {
    case 'b':
    case 'f':
    case 'p':
    case 'v':
    case 'B':
    case 'F':
    case 'P':
    case 'V':
      return '1';
    case 'c':
    case 'g':
    case 'j':
    case 'k':
    case 'q':
    case 's':
    case 'x':
    case 'z':
    case 'C':
    case 'G':
    case 'J':
    case 'K':
    case 'Q':
    case 'S':
    case 'X':
    case 'Z':
      return '2';
    case 'd':
    case 't':
    case 'D':
    case 'T':
      return '3';
    case 'l':
    case 'L':
      return '4';
    case 'm':
    case 'n':
    case 'M':
    case 'N':
      return '5';
    case 'r':
    case 'R':
      return '6';
    case 'h':
    case 'w':
    case 'H':
    case 'W':
      return 'h';
    case 'a':
    case 'e':
    case 'i':
    case 'o':
    case 'u':
    case 'y':
    case 'A':
    case 'E':
    case 'I':
    case 'O':
    case 'U':
    case 'Y':
      return 'a'; // use 'a' to encode vowels
    case '!': // we will allow all printable ASCII characters.
    case '"':
    case '#':
    case '$':
    case '%':
    case '&':
    case '\'':
    case '(':
    case ')':
    case '*':
    case '+':
    case ',':
    case '-':
    case '.':
    case '/':
    case ':':
    case ';':
    case '<':
    case '=':
    case '>':
    case '?':
    case '@':
    case '[':
    case '\\':
    case ']':
    case '^':
    case '_':
    case '`':
    case '{':
    case '|':
    case '}':
    case '~':
    case '0':
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
      return 'h'; // ignored characters; voiceless symbols.
    default:
      return 'w'; // other characters are ignored with a warning
  }
}

// Translate a string to a soundex phonetic code
//
// str: the input string
// str_len: the length of the input string
// result: the character vector in which the soundex code is copied. This 
//    should be a vector of a least length 4.
// output: the number of non-ascii or non-printable ascii characters
// encountered during translation.
unsigned int soundex(const unsigned int* str, unsigned int str_len, char* result) {
  if (!str || !result) return 0;
  if (str_len == 0) {
    unsigned int j;
    for (j = 0; j < 4; ++j) result[j] = '0';
    return 0;
  }
  unsigned int i = 0, j = 0, nfail = 0;
  char cj = translate_soundex(str[j]);
  // the first character is copied directly and not translated to a numerical
  // code
  result[0] = toupper(str[0]);
  for (i = 1; i < str_len && j < 4; ++i) {
    char ci = translate_soundex(str[i]);
    if (ci == 'a') {
      // vowels are not added to the result; but we do set the previous
      // character to the vowel because two consonants with a vowel in between
      // are not merged
      cj = ci;
    } else if (ci != 'h') {
      // a consonant that is not equal to the previous consonant is added to 
      // the result
      if (ci != cj) {
        result[++j] = ci;
        cj = ci;
      }
    }
    if ( ci == 'w' ){
      // the translated character is non-printable ASCII or non-ASCII.
      ++nfail;
    }
  }
  // pad with zeros
  for (++j ; j < 4; ++j) result[j] = '0';
  return nfail;
}

double soundex_dist(unsigned int *a, unsigned int *b, unsigned int a_len, 
    unsigned int b_len, unsigned int *nfail) {
  char sa[4];
  char sb[4];
  (*nfail) += soundex(a, a_len, sa);
  (*nfail) += soundex(b, b_len, sb);
  int c = strncmp(sa, sb, 4);
  return (c != 0)*1.0;
}

// ================================ R INTERFACE ===============================

void charwarning(unsigned int nfail){
  warning("soundex encountered %d non-printable ASCII or non-ASCII"
  "\n  characters. Results may be unreliable, see ?printable_ascii",nfail);
}

SEXP R_soundex(SEXP x) {
  int n = length(x);
  int bytes = IS_CHARACTER(x);

  // when a and b are character vectors; create unsigned int vectors in which
  // the elements of and b will be copied
  unsigned int *s = NULL;
  if (bytes) {
    int ml = max_length(x);
    s = (unsigned int *) malloc(ml*sizeof(unsigned int));
    if (s == NULL) {
       free(s);
       error("Unable to allocate enough memory");
    }
  }

  // create output variable
  SEXP y = allocVector(STRSXP, n);
  PROTECT(y);

  // compute distances, skipping NA's
  int len_s, isna_s;
  unsigned int nfail = 0;
  char sndx[5];
  for (int i = 0; i < n; ++i) {
    s = get_elem(x, i, bytes, &len_s, &isna_s, s);
    if (isna_s) {
      SET_STRING_ELT(y, i, R_NaString);
    } else { 
      nfail += soundex(s, len_s, sndx);
      sndx[4] = 0;
      SET_STRING_ELT(y, i, mkChar(sndx));
    } 
  }
  if ( nfail > 0 ) charwarning(nfail);
  // cleanup and return
  if (bytes) free(s);
  UNPROTECT(1);
  return y;
}


SEXP R_soundex_dist(SEXP a, SEXP b) {
  int na = length(a);
  int nb = length(b);
  int nt = MAX(na,nb);
  int bytes = IS_CHARACTER(a);

  // when a and b are character vectors; create unsigned int vectors in which
  // the elements of and b will be copied
  unsigned int *s = NULL, *t = NULL;
  if (bytes) {
    int ml_a = max_length(a);
    int ml_b = max_length(b);
    s = (unsigned int *) malloc((ml_a + ml_b) * sizeof(unsigned int));
    t = s + ml_a;
    if (s == NULL) {
       free(s);
       error("Unable to allocate enough memory");
    }
  }

  // create output variable
  SEXP yy = allocVector(REALSXP, nt);
  PROTECT(yy);
  double *y = REAL(yy);

  // compute distances, skipping NA's
  int len_s, len_t, isna_s, isna_t;
  unsigned int nfail = 0;
  for (int k=0; k < nt; ++k, ++y) {
    s = get_elem(a, k % na, bytes, &len_s, &isna_s, s);
    t = get_elem(b, k % nb, bytes, &len_t, &isna_t, t);
    if (isna_s || isna_t) {
      (*y) = NA_REAL;
    } else { 
      (*y) = soundex_dist(s, t, len_s, len_t, &nfail);
    } 
  }

  if ( nfail > 0 ) charwarning(nfail);

  // cleanup and return
  if (bytes) free(s);
  UNPROTECT(1);
  return yy;
}

SEXP R_match_soundex(SEXP x, SEXP table, SEXP nomatch, SEXP matchNA) {

  int nx = length(x);
  int ntable = length(table);
  int no_match = INTEGER(nomatch)[0];
  int match_na = INTEGER(matchNA)[0];
  int bytes = IS_CHARACTER(x);

  // when a and b are character vectors; create unsigned int vectors in which
  // the elements of and b will be copied
  unsigned int *s = NULL, *t = NULL;
  if (bytes) {
    int ml_x = max_length(x);
    int ml_t = max_length(table);
    s = (unsigned int *) malloc((ml_x + ml_t) * sizeof(unsigned int));
    t = s + ml_x;
    if (s == NULL) {
       free(s);
       error("Unable to allocate enough memory");
    }
  }

  // output vector
  SEXP yy = allocVector(INTSXP, nx);
  PROTECT(yy);
  int* y = INTEGER(yy);

  int index, isna_s, isna_t, len_s, len_t;
  unsigned int nfail = 0;
  double d;
  for (int i=0; i<nx; ++i) {
    index = no_match;
    s = get_elem(x, i, bytes, &len_s, &isna_s, s);

    for (int j=0; j<ntable; ++j) {
      t = get_elem(table, j, bytes, &len_t, &isna_t, t);

      if (!isna_s && !isna_t) {        // both are char (usual case)
        d = soundex_dist(s, t, len_s, len_t, &nfail);
        if (d < 0.5) { // exact match as d can only take on values 0 and 1
          index = j + 1;
          break;
        } 
      } else if (isna_s && isna_t) {  // both are NA
        index = match_na ? j + 1 : no_match;
        break;
      }
    }
    y[i] = index;
  }   

  if ( nfail > 0 ) charwarning(nfail);
  if (bytes) free(s); 
  UNPROTECT(1);
  return(yy);
}


