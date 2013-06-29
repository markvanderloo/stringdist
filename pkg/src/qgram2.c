/* a slightly more advanced implementation of the qgram distance.
 * q-grams are pushed onto a binairy tree, which is kept over one call
 * of stringdist (which loops over string pairs).
 */

//#define USE_RINTERNALS
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rdefines.h>
#include "utils.h"
/* binary tree; dictionary of qgrams */

typedef struct qnode {
  unsigned int *qgram;
  unsigned int n[2]; // nr of occurrences of qgram in s and t
  struct qnode *left;
  struct qnode *right;
} qtree;

static void free_qtree(qtree *Q){
  if ( Q == NULL ) return;
  free_qtree(Q->left);
  free_qtree(Q->right);
  free(Q->qgram);
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

/* push qgram into binary tree */
static qtree *push(qtree *Q, unsigned int *qgram, unsigned int q, int location){
  int cond;  
  if( Q == NULL ){ // new qgram
    Q = (qtree *) malloc(sizeof(qtree));
    if ( Q == NULL ) return NULL;
    Q->qgram = (unsigned int *) malloc(sizeof(int) * q);
    if (Q == NULL ) return(NULL);
    Q->n[0] = 0;
    Q->n[1] = 0;
    Q->n[location]++;
    memcpy(Q->qgram, qgram, sizeof(int) * q);
    Q->left = NULL;
    Q->right= NULL;
  } else if ( ( cond = compare(qgram, Q->qgram, q) ) == 1)  { // qgram larger than the stored qgram
    Q->left = push(Q->left, qgram, q, location);
  } else if ( cond == -1 ){ // qgram smaller than the stored qgram
    Q->right = push(Q->right, qgram, q, location);
  } else { // qgram equal to stored qgram
    Q->n[location] += 1;
  }
  return Q;
}

/* push qgrams of a string into binary tree */
static qtree *push_string(unsigned int *str, int strlen, unsigned int q, qtree *Q, int location){
  qtree *P;
  for ( int i=0; i < strlen - q + 1; ++i ){
    P = push(Q, str + i, q, location);
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
 * return values:
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
  if ( q > (x <= y ? x : y) ) return -1;
  // rare edge cases.
  if ( q == 0 ){
    if ( x + y > 0 ){ // distance undefined
      return -1;
    } else { // x == y == 0.
      return 0;
    } 
  }

  double dist[3] = {0,0,0};

  Q = push_string(s, x, q, Q, 0);
  if (Q == NULL) return -2;
  Q = push_string(t, y, q, Q, 1);
  if (Q == NULL) return -2;

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
  if ( q < 0 ){
    UNPROTECT(4);
    error("q must be a nonnegative integer");
  }
  // choose distance function
  int dist = INTEGER(distance)[0];
  if ( dist < 0 || dist > 2 ){
    UNPROTECT(4);
    error("unkown distance function");
  } 

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
    if (y[k] == -2){
      error("Could not allocate enough memory");
    }
    i = RECYCLE(i+1,na);
    j = RECYCLE(j+1,nb);
  }
  free_qtree(Q);
  UNPROTECT(5);
  return yy;
}


/* R interface to qgram tabulator */

static void count_qtree(qtree *Q, int *n){
  if (Q == NULL ) return ;
  n[0]++;
  count_qtree(Q->left,n);
  count_qtree(Q->right,n);
}

static void get_counts(qtree *Q, int q, int *qgrams, int *count, int *index){
  if ( Q == NULL ) return ;
  memcpy(qgrams + q*index[0], Q->qgram, sizeof(int) * q);
  count[index[0]] = Q->n[0];
  ++index[0];
  get_counts(Q->left, q, qgrams, count, index);
  get_counts(Q->right,q,qgrams, count, index);
}

/* TODO:
 * 
 * - Detect memory allocation failure.
 */
SEXP R_get_qgrams(SEXP a, SEXP qq){
  PROTECT(a);
  PROTECT(qq);

  int q = INTEGER(qq)[0];
  int n = length(a);
  if ( q < 0 ){
    UNPROTECT(2);
    return R_NilValue;
  }

  if ( q < 0 ){
    UNPROTECT(2);
    error("q must be a nonnegative integer");
  }
  // set up a tree
  qtree *Q = NULL;
  for ( int i=0; i < n; ++i ){
    if ( INTEGER(VECTOR_ELT(a,i))[0] == NA_INTEGER 
        || q > length(VECTOR_ELT(a,i)) 
        || ( q == 0 && length(VECTOR_ELT(a,i)) > 0 )
      ){
      continue ;
    }
    Q = push_string(
      (unsigned int *) INTEGER(VECTOR_ELT(a,i)),
      length(VECTOR_ELT(a,i)),
      q, Q, 0 
    );
  }
  // pick and delete the tree

  int nqgram[1] = {0};
  int index[1] = {0};  
  
  count_qtree(Q,nqgram);  

  SEXP qgrams, qcount;
  PROTECT(qgrams = allocVector(INTSXP, q*nqgram[0]));
  PROTECT(qcount = allocVector(INTSXP, nqgram[0]));

  get_counts(Q, q, INTEGER(qgrams), INTEGER(qcount),index);
  
  setAttrib(qcount, install("qgrams"), qgrams);
  
  free_qtree(Q);
  UNPROTECT(4);

  return(qcount);
}



