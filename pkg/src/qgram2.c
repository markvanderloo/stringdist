/* a slightly more advanced implementation of the qgram distance.
 * q-grams are pushed onto a binairy tree, which is kept over one call
 * of stringdist (which loops over string pairs).
 */

#include<stdlib.h>
#include<string.h>
#include<stdio.h>
#include<R.h>
#include<Rdefines.h>

/* sorted list; dictionary of qgrams */

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
  compare( q1 + 1, q2 + 1, q - 1 );
}


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

// get qgram-distance from tree and set all qgram-freqencies to 0.
static void getdist(qtree *Q, int *d){
  if (Q == NULL) return;
  d[0] = d[0] + abs(Q->n[0] - Q->n[1]);
  Q->n[0] = 0;
  Q->n[1] = 0;
  getdist(Q->left, d);
  getdist(Q->right,d);
}


static int qgram2(
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
  qtree *P;

  for (int i=0; i < x - q + 1; ++i ){
    P = push(Q, s + i, q, 0);
    
    if ( P == NULL ){
      free_qtree(Q);
      return -2;
    } else {
      Q = P;
    }

  }
  for (int i=0; i < y - q + 1; ++i ){
    Q = push(Q, t + i, q, 1);

    if ( P == NULL ){
      free_qtree(Q);
      return -2;
    } else {
      Q = P;
    }

  }

  getdist(Q,dist);
  return dist[0];
}


SEXP R_qgram2(SEXP a, SEXP b, SEXP qq){
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
  PROTECT(yy = allocVector(REALSXP, nt));
  double *y = REAL(yy);

  // set up a qtree;
  qtree *Q = NULL;

  for ( k=0; k < nt; ++k ){
    i = k % na;
    j = k % nb;
    if (INTEGER(VECTOR_ELT(a,i))[0] == NA_INTEGER || INTEGER(VECTOR_ELT(b,j))[0] == NA_INTEGER){
      y[k] = NA_REAL;
      continue;
    }
    y[k] = (double) qgram2(
        INTEGER(VECTOR_ELT(a,i)),
        INTEGER(VECTOR_ELT(b,j)),
        length(VECTOR_ELT(a,i)),
        length(VECTOR_ELT(b,j)),
        q,
        Q
    );
    if (y[k] == -2){
      error("Could not alocate enough memory");
    }
  }
  free_qtree(Q);
  UNPROTECT(3);
  return yy;
}





/*

static void qprint(qtree *Q, int q){
  for ( int i=0; i < q; ++i ){
    printf("%d ",Q->qgram[i]);
  }
    
    printf(": s(%d) t(%d) \n",Q->n[0], Q->n[1]);
}

static void treeprint(qtree *Q, int q){
  if ( Q == NULL ) return;
  qprint(Q,q);
  treeprint(Q->left, q);
  treeprint(Q->right, q);
}


void main(){
  qtree *Q = NULL;
  int foo [] = {1,2,3};
  int bar [] = {3,2,1};
  Q = push(Q, foo, 3,0);
  Q = push(Q, bar, 3,1);
  Q = push(Q, bar, 3,0);
  Q = push(Q, foo, 3,0);

  treeprint(Q,3);
  int d[1] = {0};
  getdist(Q,d);
  printf("dist : %d\n",d[0]);
}

*/
