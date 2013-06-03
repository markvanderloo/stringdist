/* a slightly more advanced implementation of the qgram distance.
 * q-grams are pushed onto a binairy tree, which is kept over one call
 * of stringdist (which loops over string pairs).
 */

#include<stdlib.h>
#include<string.h>
#include<stdio.h>
#include<R.h>
#include<Rdefines.h>

/* binary tree; dictionary of qgrams */

typedef struct qnode {
  unsigned int *qgram;
  unsigned int n[2]; // nr of occurrences of qgram in s and t
  struct qnode *left;
  struct qnode *right;
} qtree;

static void free_qtree(qtree *Q){
  if (Q==NULL) return;
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
static void getdist(qtree *Q, int *d){
  if (Q == NULL) return;
  d[0] = d[0] + abs(Q->n[0] - Q->n[1]);
  Q->n[0] = 0;
  Q->n[1] = 0;
  getdist(Q->left, d);
  getdist(Q->right,d);
}

/*Get qgram distances 
 * return values:
 *  >=0 : qgram distance
 * -1   : infinite distance
 * -2   : Not enough memory
 */
static int qgram_tree(
    unsigned int *s, 
    unsigned int *t, 
    unsigned int x,
    unsigned int y,
    unsigned int q, 
    qtree *Q
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

  int dist[1] = {0};

  Q = push_string(s, x, q, Q, 0);
  if (Q == NULL) return -2;
  Q = push_string(t, y, q, Q, 1);
  if (Q == NULL) return -2;

  getdist(Q,dist);
  return dist[0];
}

/* R interface */
SEXP R_qgram_tree(SEXP a, SEXP b, SEXP qq){
  PROTECT(a);
  PROTECT(b);
  int q = INTEGER(qq)[0];
  if ( q < 0 ){
    UNPROTECT(2);
    error("q must be a nonnegative integer");
  } 
  int i, j, k;
  int na = length(a);
  int nb = length(b);
  int nt = (na > nb) ? na : nb;

  SEXP yy; 
  PROTECT(yy = allocVector(INTSXP, nt));
  int *y = INTEGER(yy);

  // set up a qtree;
  qtree *Q = NULL;

  for ( k=0; k < nt; ++k ){
    i = k % na;
    j = k % nb;
    if (INTEGER(VECTOR_ELT(a,i))[0] == NA_INTEGER || INTEGER(VECTOR_ELT(b,j))[0] == NA_INTEGER){
      y[k] = NA_REAL;
      continue;
    }
    y[k] = qgram_tree(
        INTEGER(VECTOR_ELT(a,i)),
        INTEGER(VECTOR_ELT(b,j)),
        length(VECTOR_ELT(a,i)),
        length(VECTOR_ELT(b,j)),
        q,
        Q
    );
    if (y[k] == -2){
      error("Could not allocate enough memory");
    }
  }
  free_qtree(Q);
  UNPROTECT(3);
  return yy;
}

void count_qtree(qtree *Q, int *n){
  if ( Q == NULL ) return ;
  n[0]++;
  count_qtree(Q->left, n);
  count_qtree(Q->right, n);
}



SEXP R_get_qgrams(SEXP a, SEXP qq){
  PROTECT(a);
  PROTECT(qq);
  
  int q = INTEGER(qq)[0];
  if ( q < 0 ){
    error("q must be a nonnegative integer");
  }
  int n = length(a);
  // set up a tree, push all qgrams of strings in a into tree. 
  qtree *Q = NULL;
  for ( int i=0; i<n; ++i){
    if (INTEGER(VECTOR_ELT(a,i))[0] == NA_INTEGER){
      continue;
    }
    Q = push_string(
      INTEGER(VECTOR_ELT(a,i)),
      length(VECTOR_ELT(a,i)),
      q, Q, 0
   );
   if (Q == NULL){
    error("Could not allocate enough memory");
   }
  }
  // pick q-grams from the tree
  int nqgrams[1] = {0};
  count_qtree(Q,nqgrams);
  // this 1d-vector represents a nqgrams X (q+1) array
  // where the first q colums represent qgrams and the q+1st colum
  // the number of qgrams.
  SEXP yy;
  PROTECT(yy = allocVector(INTSXP, nqgrams*(q+1));
  list_qgrams(Q, INTEGER(yy), 0);
  free_qtree(Q);
  UNPROTECT(3);
  return(yy);
}


