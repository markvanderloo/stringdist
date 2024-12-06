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

#ifdef _OPENMP
#include <omp.h>
#endif

#define USE_RINTERNALS
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "qtree.h"



/* -- Simple memory allocator to store nodes of the qtree -- */


/* Nodes are stored in boxes which are stored on a shelf.
 * Every time a new box is added to the shelf, the capacity 
 * for node storage doubles, unless MAXBOXES is surpassed.
 */
#define MAXBOXES 20             // max number of boxes
#define MIN_BOX_SIZE (1<<10)    // nr of nodes in initial box
#define MAX_NUM_THREADS (1<<10) // nr of threads for which alocator is thread-safe.

// A Box of nodes.
typedef struct {
  
  int nnodes; // number of nodes
  int nalloc; // number of nodes allocated

  // Room in the box.
  unsigned int  *intblocks; // store qgrams
  double        *dblblocks; // store qgram-counts
  qtree         *qtrblocks; // store nodes
} Box;


static Box *new_box(int nnodes, int q, int nstr){
  Box *b = malloc(sizeof(Box));
  if ( b == NULL ) 
    return NULL;
  b->nnodes     = nnodes;
  b->nalloc     = 0L;
  b->intblocks  = (unsigned int *) malloc(sizeof(int) * nnodes * q);
  b->dblblocks  = (double *) malloc(sizeof(double) * nnodes * nstr);
  b->qtrblocks  = (qtree *) malloc(sizeof(qtree) * nnodes);
  return b;
}

static void free_box(Box *box){
  free( box->intblocks);
  free( box->dblblocks);
  free( box->qtrblocks);
  free( box);
}

// A shelf can store up to MAXBOXES Boxes.
typedef struct {
  Box *box[MAXBOXES];
  int nboxes;         // number of boxes on the shelf
  int q;              // the q in q-gram
  int nstr;           // the number of strings compared
} Shelf;

// A wall with shelfs: one for each thread.
static Shelf wall[MAX_NUM_THREADS];

// When multithreaded, check what shelf we're storing stuff.
static inline int get_shelf_num(void){
  int thread_num=0;
  #ifdef _OPENMP
  thread_num = omp_get_thread_num();
  #endif
  return thread_num;
}



static void init_shelf(int q, int nstr){

  Shelf *shelf = &wall[get_shelf_num()];

  shelf->q = q;
  shelf->nstr = nstr;
  shelf->nboxes = 0L;
  for ( int i=0; i<MAXBOXES; i++ )
    shelf->box[i] = NULL;
}

/* add box to shelf
 *
 * return values:
 * 0: okidoki
 * 1: no luck: out of memory or MAXBOXES exceeded
 *
 */
static int add_box(int nnodes){

  Shelf *shelf = &wall[get_shelf_num()];
  // is there room for another box?
  if ( shelf->nboxes >= MAXBOXES ) return 1;

  Box *b = new_box(nnodes, shelf->q, shelf->nstr);
  if ( b != NULL ){
    shelf->box[shelf->nboxes] = b;
    shelf->nboxes++;
    return 1L;
  } else {
    return 0L;
  }
  
}

static void clear_shelf(void){
  Shelf *shelf = &wall[get_shelf_num()];
  for ( int i = 0; i < shelf->nboxes; i++ ){
    free_box(shelf->box[i]);
  }
  shelf->nboxes=0L;
}

/* the allocator can store ints (q-grams), doubles (counts, per string) 
 * and qtree (nodes).
 */
typedef enum { uInt, Double, Qtree } type;


// cf. n1256.pdf (C99 std) sect 6.3.2.3 for pointer conversion.
static void *alloc(type t){
  Shelf *shelf = &wall[get_shelf_num()];


  if ( 
    shelf->nboxes == 0L &&
    !add_box(MIN_BOX_SIZE)
  ) return NULL;

  Box *box = shelf->box[shelf->nboxes-1L];
  if ( box->nalloc == box->nnodes ){
    // add box such that storage size is doubled.
    if (
      !add_box( (1 << (shelf->nboxes-1L)) * MIN_BOX_SIZE )
    ) return NULL;
    box = shelf->box[shelf->nboxes-1L];
  }

  void *x;
  switch ( t){
    case uInt:
      x = (void *) (box->intblocks + box->nalloc * shelf->q);
      break;
    case Double:
      x = (void *) (box->dblblocks + box->nalloc * shelf->nstr);
      break;
    case Qtree:

      x = (void *) (box->qtrblocks + box->nalloc);
      box->nalloc += 1;
      break;
    default:
      return NULL;
      break;
  }
  
  return x;
}


/* -- END OF ALLOCATOR -- */


// helper functions
qtree *new_qtree(int q, int nstr){
  init_shelf(q, nstr);
  return NULL;
}

void free_qtree(void){
  clear_shelf();
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
 * node  : int array of length nLoc, will contain contents of a node.
 */
static qtree *push(qtree *Q, unsigned int *qgram, unsigned int q, int iLoc, int nLoc, double *node ){
  int cond;  
  if( Q == NULL ){ // new qgram
    Q = (qtree *) alloc( Qtree);
    if ( Q == NULL ) return NULL;
    Q->qgram = (unsigned int *) alloc( uInt);
    if (Q->qgram == NULL ) return NULL;

    Q->n = (double *) alloc( Double);
    if (Q->n == NULL) return NULL;
    for (int i=0; i<nLoc; ++i) Q->n[i] = 0.0;

    Q->n[iLoc]++;
    memcpy(Q->qgram, qgram, sizeof(int) * q);
    Q->left = NULL;
    Q->right= NULL;
    // copy content for all iLocs to output parameter.
    if (node != NULL ) memcpy(node, Q->n, sizeof(double) * nLoc);

  } else if ( ( cond = compare(qgram, Q->qgram, q) ) == 1)  { // qgram larger than the stored qgram
    Q->left = push(Q->left, qgram, q, iLoc, nLoc, node);
  } else if ( cond == -1 ){ // qgram smaller than the stored qgram
    Q->right = push(Q->right, qgram, q, iLoc, nLoc, node);
  } else { // qgram equal to stored qgram
    Q->n[iLoc] += 1;
    // copy content for all iLocs to output parameter.
    if (node != NULL ) memcpy(node,Q->n, sizeof(double) * nLoc);
  }
  return Q;
}

/* pull qgram from binary tree: decrease valaue for one of the strings.
 * 
 * qtree : see above
 * qgram : see above
 * q     : the 'q' in q-gram
 * iLoc  : To wich count location does this q-gram contribute?
 * nLoc  : how many locations are there?
 * node  : int array of length nLoc, will contain contents of a node.
 */
static qtree *pull(qtree *Q, unsigned int *qgram, unsigned int q, int iLoc, int nLoc, double *node){
  if (Q == NULL) return(NULL);
  int cond = compare(qgram, Q->qgram, q);

  if ( cond == -1 ){ // qgram smaller than stored qgram
    pull(Q->right, qgram, q, iLoc, nLoc, node);
  } else if (cond == 1) { //qram larger than stored qgram
    pull(Q->left, qgram, q, iLoc, nLoc, node);
  } else { // qgram equal to stored qgram
    Q->n[iLoc] -= 1;
    // copy content for all iLocs to output parameter.
    if (node != NULL) memcpy(node, Q->n, sizeof(double) * nLoc);
  }
  return Q;
}


/* push qgrams of a string into binary tree */
static qtree *push_string(unsigned int *str, int strlen, unsigned int q, qtree *Q, int iLoc, int nLoc){
  qtree *P;
  for ( int i=0; i < (int) (strlen - q + 1); i++ ){
    P = push(Q, str + i, q, iLoc, nLoc, NULL);
    if ( P == NULL ){ 
      free_qtree();
      return NULL;
    }
    Q = P;
  }
  return Q;
}




/* get qgram-distance from tree and set all qgram-freqencies 
 * to 0 (so the tree can be reused).
 */
static void getdist(qtree *Q, double *d){
  if (Q == NULL) return;
  d[0] += fabs(Q->n[0] - Q->n[1]);
  Q->n[0] = 0;
  Q->n[1] = 0;
  getdist(Q->left, d);
  getdist(Q->right,d);
}

/* get x.y,||x||and ||y|| for cosine distance from the tree and set all qgram-freqencies 
 * to 0 so the tree van be reused.
 */
static void getcosine(qtree *Q, double *d, int clean){
  if ( Q == NULL ) return;
  // inner product
  d[0] += (double) Q->n[0] * Q->n[1];
  // norm of v(s,q)
  d[1] += (double) Q->n[0]*Q->n[0];
  d[2] += (double) Q->n[1]*Q->n[1];
  // clean up and continue
  if (clean){
    Q->n[0] = 0;
    Q->n[1] = 0;
  }
  getcosine(Q->left, d, clean);
  getcosine(Q->right, d, clean);
}

static double cosdist(double xy, double xx, double yy){
  // x and y are equal: return precisely zero.
  if (xy == xx && xy == yy){
    return 0.0;
  } else {
  // use fabs to avoid numerical -0.
    return ( fabs(1.0- xy/( sqrt(xx) * sqrt(yy))) );
  }
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
  if ( Q->n[0] > 0 || Q->n[1] > 0){
    ++d[1];
  }
  // clean up and continue
  Q->n[0] = 0;
  Q->n[1] = 0;
  getjaccard(Q->left,d);
  getjaccard(Q->right,d);
}

/* for testing purposes only
static void print_qtree(qtree *Q, int q){
  if (Q==NULL) return;
  Rprintf("q=%d ",q);
  Rprintf("qgram = {");
  for(int i = 0; i < q; i++)
    Rprintf("%03d ",Q->qgram[i]);
  Rprintf("}");
  Rprintf("n = [%2.0f %2.0f]\n", Q->n[0], Q->n[1]);
  print_qtree(Q->left,q);
  print_qtree(Q->right,q);
}*/

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
double qgram_dist(
    unsigned int *s, 
    int x,
    unsigned int *t, 
    int y,
    unsigned int q, 
    qtree **Qp,
    int distance
  ){

  // rare edge case: q==0. Note that we return 0 for all cases where
  // q equals zero. In the R journal paper we used Inf for cases where 
  // q=0 and |s| or |t| > 0
  if ( q == 0 ) return 0.0;

  double dist[3] = {0,0,0};
  *Qp = push_string(s, x, q, *Qp, 0, 2);

  *Qp = push_string(t, y, q, *Qp, 1, 2);
  if (*Qp == NULL) return 0;
  

  qtree *Q = *Qp;
  switch ( distance ){
    case 0:
      getdist(Q,dist);
      break;
    case 1:
      getcosine(Q, dist, 1);
      dist[0] = cosdist(dist[0],dist[1],dist[2]);
      break;
    case 2:
      getjaccard(*Qp,dist);
      dist[0] = 1.0 - dist[0]/dist[1];
      break;
    default:
      break;
  }

  return dist[0];
}




/*
* s: text to search
* x: length of s
* t: pattern
* y: length of pattern
* q: size of q-gram
* qtree: a qtree object.
* store: length 3 array to store intermediate values;
*/
double running_cosine_dist(
  unsigned int *s,
  int x,
  unsigned int *t,
  int y,
  unsigned int q,
  qtree **Qp,
  double *store
  ){
  
  double d;

  unsigned int *first_qgram;
  unsigned int *last_qgram;

  // pwi: value of qgram table of pattern and window at location 
  //      where one qgram is removed. 
  // pwj: value of qgram table of pattern and window at location
  //      where one qgram is added.
  double pwi[2] = {0.,0.}, pwj[2] = {0.,0.};

  
  if ( *Qp == NULL ){ // new tree, 
    // push the search pattern, location 0
    *Qp = push_string(t, y, q, *Qp, 0, 2);
    // push the first window
    *Qp = push_string(s, x, q, *Qp, 1, 2);
    store[0] = store[1] = store[2] = 0;
    // store[0]: w.p (inner product)
    // store[1]: p.p (squared norm of pattern)
    // store[2]: w.w (squared norm of window)
    getcosine(*Qp, store, 0);
    d = cosdist(store[0], store[1], store[2]); 
  } else { // we are running
    first_qgram = s - 1;
    last_qgram  = s + y - q;
    // special case: q-gram to remove is equal to qgram to add
    if (compare(first_qgram, last_qgram, q) == 0){
      d = cosdist(store[0], store[1], store[2]); 
    } else { 
      // take first q-gram of the previous window from the table.
      *Qp = pull(*Qp, s-1, q, 1, 2, pwi);
      // add last qgram of the current window to the table.
      *Qp = push(*Qp, s+y-q, q, 1, 2, pwj);
     
      store[0] = store[0] - pwi[0] + pwj[0];
      store[2] = store[2] + 2*(pwj[1] - pwi[1] - 1);
      d  = cosdist(store[0], store[2], store[1]);
    }
  }

  return d;

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

/* 
 * 
 * 
 */
SEXP R_get_qgrams(SEXP a, SEXP qq){

  int q = INTEGER(qq)[0];

  if ( q < 0 ){
    error("q must be a nonnegative integer");
  }


  SEXP strlist;
  int nstr, nchar, nLoc = length(a);
  unsigned int *str;
  
  // set up a tree; push all the qgrams.
  qtree *Q = new_qtree( q, nLoc);
  
  for ( int iLoc = 0; iLoc < nLoc; ++iLoc ){
    strlist = VECTOR_ELT(a, iLoc);
    nstr    = length(strlist);
 
    for ( int i=0; i < nstr; ++i ){
      str   = (unsigned int *) INTEGER(VECTOR_ELT(strlist,i));
      nchar = length(VECTOR_ELT(strlist,i));
      if ( nchar == 0
          || str[0] == NA_INTEGER 
          || q > nchar
          || ( q == 0 && nchar > 0 )
        ){
        continue ;
      }
      Q = push_string(str, nchar, q, Q, iLoc, nLoc);
      if ( Q == NULL ){
        error("could not allocate enough memory");
      }
    }
  }
  // pick and delete the tree

  int nqgram[1] = {0};

  // helper variable for get_counts 
  int index[1] = {0};

  count_qtree(Q,nqgram);  

  SEXP qgrams, qcount;
  qgrams = PROTECT(allocVector(INTSXP, q*nqgram[0]));
  qcount = PROTECT(allocVector(REALSXP, nLoc*nqgram[0]));

  get_counts(Q, q, INTEGER(qgrams), nLoc, index, REAL(qcount));
  
  setAttrib(qcount, install("qgrams"), qgrams);
  
  free_qtree();
  UNPROTECT(2);

  return(qcount);
}



