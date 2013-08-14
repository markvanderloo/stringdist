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


#define USE_RINTERNALS
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rdefines.h>
#include "utils.h"
/* binary tree; dictionary of qgrams */

typedef struct qnode {
  // the payload q-gram
  unsigned int *qgram;
  // vector of q-gram counts 
  // (double to prevent int overload 
  // when calulating distances)
  double *n; 
  struct qnode *left;
  struct qnode *right;
} qtree;

static void free_qtree(qtree *Q){
  if ( Q == NULL ) return;
  free_qtree(Q->left);
  free_qtree(Q->right);
  free(Q->qgram);
  free(Q->n);
  free(Q);
}

/* Lexicographical comparison of two qgrams.
 * output:
 * -1 : q1 < q2
 *  0 : q1 = q2
 *  1 : q1 > q2
 */
static int compare(unsigned int *q1, unsigned int *q2, int q){
  if (q==0) return 0;
  if (q1[0] > q2[0]) return 1;
  if (q1[0] < q2[0]) return -1;
  return compare( q1 + 1, q2 + 1, q - 1 );
}

/* push qgram into binary tree
 * 
 * qtree : see above
 * qgram : see above
 * q     : the 'q' in q-gram
 * iLoc  : To wich count location does this q-gram contribute?
 * nLoc  : how many locations are there?
 */
static qtree *push(qtree *Q, unsigned int *qgram, unsigned int q, int iLoc, int nLoc ){
  int cond;  
  if( Q == NULL ){ // new qgram

    Q = (qtree *) malloc(sizeof(qtree));
    if ( Q == NULL ) return NULL;

    Q->qgram = (unsigned int *) malloc(sizeof(int) * q);
    if (Q->qgram == NULL ) return NULL;

    Q->n = (double *) malloc(sizeof(double) * nLoc);
    if (Q->n == NULL) return NULL;
    for (int i=0; i<nLoc; ++i) Q->n[i] = 0.0;

    Q->n[iLoc]++;
    memcpy(Q->qgram, qgram, sizeof(int) * q);
    Q->left = NULL;
    Q->right= NULL;

  } else if ( ( cond = compare(qgram, Q->qgram, q) ) == 1)  { // qgram larger than the stored qgram
    Q->left = push(Q->left, qgram, q, iLoc, nLoc);
  } else if ( cond == -1 ){ // qgram smaller than the stored qgram
    Q->right = push(Q->right, qgram, q, iLoc, nLoc);
  } else { // qgram equal to stored qgram
    Q->n[iLoc] += 1;
  }
  return Q;
}

/* push qgrams of a string into binary tree */
static qtree *push_string(unsigned int *str, int strlen, unsigned int q, qtree *Q, int iLoc, int nLoc){
  qtree *P;
  for ( int i=0; i < strlen - q + 1; ++i ){
    P = push(Q, str + i, q, iLoc, nLoc);
    if ( P == NULL ){ 
      free_qtree(Q);
      return NULL;
    }
    Q = P;
  }
  return Q;
}



/* The real work starts here */

/* get qgram-distance from tree and set all qgram-freqencies 
 * to 0 (so the tree can be reused).
 */
static void getdist(qtree *Q, double *d){
  if (Q == NULL) return;
  d[0] += (double) abs(Q->n[0] - Q->n[1]);
  Q->n[0] = 0;
  Q->n[1] = 0;
  getdist(Q->left, d);
  getdist(Q->right,d);
}

/* get x.y,||x||and ||y|| for cosine distance from the tree and set all qgram-freqencies 
 * to 0 so the tree van be reused.
 */
static void getcosine(qtree *Q, double *d){
  if ( Q == NULL ) return;
  // inner product
  d[0] += (double) Q->n[0] * Q->n[1];
  // norm of v(s,q)
  d[1] += (double) Q->n[0]*Q->n[0];
  d[2] += (double) Q->n[1]*Q->n[1];
  // clean up and continue
  Q->n[0] = 0;
  Q->n[1] = 0;
  getcosine(Q->left,d);
  getcosine(Q->right,d);
}

/* get jaccard distance from the tree and set all qgram-freqencies 
 * to 0 so the tree van be reused.
 */
static void getjaccard(qtree *Q, double *d){
  if ( Q == NULL ) return;
  // numerator: |x A y|
  if ( Q->n[0] > 0 && Q->n[1] > 0){
    ++d[0];
  } 
  // denominator: |x V y|
  ++d[1];
  // clean up and continue
  Q->n[0] = 0;
  Q->n[1] = 0;
  getjaccard(Q->left,d);
  getjaccard(Q->right,d);
}



/*Get qgram distances 
 * Input
 * s: a string
 * t: a string
 * x: length of s
 * y: length of t
 * q: the 'q' in q-gram
 * Q: a qtree
 * int: distance distance function to compute:
 *  0 : q-gram distance
 *  1 : cosine distance
 *  2 : jaccard distance
 *
 *
 * Return values:
 *  >=0 : qgram distance
 * -1   : infinite distance
 * -2   : Not enough memory
 */
static double qgram_tree(
    unsigned int *s, 
    unsigned int *t, 
    unsigned int x,
    unsigned int y,
    unsigned int q, 
    qtree *Q,
    int distance
  ){
  // return -1 when q is larger than the length of the shortest string.
  if ( q > (x <= y ? x : y) ) return -1.0;
  // rare edge cases.
  if ( q == 0 ){
    if ( x + y > 0 ){ // distance undefined
      return -1.0;
    } else { // x == y == 0.
      return 0.0;
    } 
  }

  double dist[3] = {0,0,0};

  Q = push_string(s, x, q, Q, 0, 2);
  if (Q == NULL) return -2.0;
  Q = push_string(t, y, q, Q, 1, 2);
  if (Q == NULL) return -2.0;

  switch ( distance ){
    case 0:
      getdist(Q,dist);
      break;
    case 1:
      getcosine(Q, dist);
      dist[0] = 1.0 - dist[0]/(sqrt(dist[1]) * sqrt(dist[2]));
      break;
    case 2:
      getjaccard(Q,dist);
      dist[0] = 1.0 - dist[0]/dist[1];
      break;
    default:
      break;
  }
  return dist[0];
}

/* R interface to qgram distance */
SEXP R_qgram_tree(SEXP a, SEXP b, SEXP qq, SEXP distance){
  PROTECT(a);
  PROTECT(b);
  PROTECT(qq);
  PROTECT(distance);
  int q = INTEGER(qq)[0];
  // choose distance function
  int dist = INTEGER(distance)[0];

  int i=0, j=0;
  int na = length(a);
  int nb = length(b);
  int nt = (na > nb) ? na : nb;

  SEXP yy; 
  PROTECT(yy = allocVector(REALSXP, nt));
  double *y = REAL(yy);

  // set up a qtree;
  qtree *Q = NULL;

  for ( int k=0; k < nt; ++k ){
    if (INTEGER(VECTOR_ELT(a,i))[0] == NA_INTEGER || INTEGER(VECTOR_ELT(b,j))[0] == NA_INTEGER){
      y[k] = NA_REAL;
      continue;
    }
    y[k] = qgram_tree(
       (unsigned int *) INTEGER(VECTOR_ELT(a,i)),
       (unsigned int *) INTEGER(VECTOR_ELT(b,j)),
        length(VECTOR_ELT(a,i)),
        length(VECTOR_ELT(b,j)),
        q,
        Q,
        dist
    );
    if (y[k] == -2.0){
      UNPROTECT(5);
      error("Could not allocate enough memory");
    }
    if (y[k] == -1.0){
      y[k] = R_PosInf;
    }
    i = RECYCLE(i+1,na);
    j = RECYCLE(j+1,nb);
  }
  free_qtree(Q);
  UNPROTECT(5);
  return yy;
}


//-- Match function interface with R

SEXP R_match_qgram_tree(SEXP x, SEXP table, SEXP nomatch, SEXP matchNA, SEXP qq, SEXP maxDist, SEXP distance){
  PROTECT(x);
  PROTECT(table);
  PROTECT(nomatch);
  PROTECT(matchNA);
  PROTECT(qq);
  PROTECT(maxDist);
  PROTECT(distance);
  int q = INTEGER(qq)[0];
  double max_dist = REAL(maxDist)[0] == 0.0 ? R_PosInf : REAL(maxDist)[0];
  
  // choose distance function
  int dist = INTEGER(distance)[0];

  int nx = length(x), ntable = length(table);
  int no_match = INTEGER(nomatch)[0];
  int match_na = INTEGER(matchNA)[0];
  

  // set up a qtree;
  qtree *Q = NULL;

  // output vector
  SEXP yy;
  PROTECT(yy = allocVector(INTSXP, nx));
  int *y = INTEGER(yy);
  int *X, *T;


  double d = R_PosInf, d1 = R_PosInf;
  int index, xNA, tNA;

  for ( int i=0; i<nx; i++){
    index = no_match;

    X = INTEGER(VECTOR_ELT(x,i));
    xNA = (X[0] == NA_INTEGER);

    for ( int j=0; j<ntable; j++){

      T = INTEGER(VECTOR_ELT(table,j));
      tNA = (T[0] == NA_INTEGER);

      if ( !xNA && !tNA ){        // both are char (usual case)
        d = qgram_tree(
          (unsigned int *) X,
          (unsigned int *) T,
          length(VECTOR_ELT(x,i)),
          length(VECTOR_ELT(table,j)),
          q,
          Q,
          dist
        );
        if ( d > max_dist ){
          continue;
        } else if ( d > -1 && d < d1){ 
          index = j + 1;
          if ( abs(d) < 1e-14 ) break; 
          d1 = d;
        }
      } else if ( xNA && tNA ) {  // both are NA
        index = match_na ? j + 1 : no_match;
        break;
      }
    }
    
    y[i] = index;
  }  
  UNPROTECT(8);
  return(yy);
}


/* R interface to qgram tabulator */

static void count_qtree(qtree *Q, int *n){
  if (Q == NULL ) return ;
  n[0]++;
  count_qtree(Q->left,n);
  count_qtree(Q->right,n);
}

static void get_counts( qtree *Q, int q, int *qgrams, int nLoc, int *index, double *count ){
  if ( Q == NULL ) return ;
  memcpy(qgrams + q*index[0], Q->qgram, sizeof(int) * q);
  memcpy(count + nLoc*index[0], Q->n, sizeof(double) * nLoc);
  ++index[0];
  get_counts(Q->left, q, qgrams, nLoc, index, count);
  get_counts(Q->right,q, qgrams, nLoc, index, count);
}

/* TODO:
 * 
 * - Detect memory allocation failure.
 */
SEXP R_get_qgrams(SEXP a, SEXP qq){
  PROTECT(a);
  PROTECT(qq);

  int q = INTEGER(qq)[0];

  if ( q < 0 ){
    UNPROTECT(2);
    error("q must be a nonnegative integer");
  }

  // set up a tree; push all the qgrams.
  qtree *Q = NULL;

  SEXP strlist;
  int nstr, nchar, nLoc = length(a);
  unsigned int *str;
  
  for ( int iLoc = 0; iLoc < nLoc; ++iLoc ){
    strlist = VECTOR_ELT(a, iLoc);
    nstr    = length(strlist);

    for ( int i=0; i < nstr; ++i ){
      str   = (unsigned int *) INTEGER(VECTOR_ELT(strlist,i));
      nchar = length(VECTOR_ELT(strlist,i));
      if ( str[0] == NA_INTEGER 
          || q > nchar
          || ( q == 0 && nchar > 0 )
        ){
        continue ;
      }
      Q = push_string(str, nchar, q, Q, iLoc, nLoc);
    }
  }
  // pick and delete the tree

  int nqgram[1] = {0};

  // helper variable for get_counts 
  int index[1] = {0};

  count_qtree(Q,nqgram);  

  SEXP qgrams, qcount;
  PROTECT(qgrams = allocVector(INTSXP, q*nqgram[0]));
  PROTECT(qcount = allocVector(REALSXP, nLoc*nqgram[0]));

  get_counts(Q, q, INTEGER(qgrams), nLoc, index, REAL(qcount));
  
  setAttrib(qcount, install("qgrams"), qgrams);
  
  free_qtree(Q);
  UNPROTECT(4);

  return(qcount);
}



