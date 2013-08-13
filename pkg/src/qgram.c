/* 
 * 
 * 
 * 
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


// only useful in debug
void print_qtree(qtree *Q){
  if (Q==NULL) return;
  Rprintf("qgram %d\n",Q->qgram[0]);
  Rprintf("n     %g\n",Q->n[0]);
  print_qtree(Q->left);
  print_qtree(Q->right);
}


/* -- Simple memory allocator to store nodes of the qtree -- */


/* Nodes are stored in boxes which are stored on a shelve.
 * Every time a new box is added to the shelve, the capacity 
 * for node storage doubles, unless MAXBOXES is surpassed.
 */
#define MAXBOXES 20
// number of nodes to initially allocate space for.
#define MIN_BOX_SIZE (1<<3)  

// A Box can store a number of nodes
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
// TODO: raise hell ico malloc failure.
  b->nnodes     = nnodes;
  b->nalloc     = 0L;
  b->intblocks  = (unsigned int *) malloc(sizeof(int) * nnodes * q);
  b->dblblocks  = (double *) malloc(sizeof(double) * nnodes * nstr);
  b->qtrblocks  = (qtree *) malloc(sizeof(qtree) * nnodes);
  return b;
}

static void free_box(Box *box){
//Rprintf("floopsiee! %d\n", box);
  // empty box
  free( box->intblocks);
// Rprintf("flapsiee!\n");
  free( box->dblblocks);
  free( box->qtrblocks);
  // throw box
  free( box);
}

// A shelve can store up to MAXBOXES Boxes.
typedef struct {
  Box *box[MAXBOXES];
  int nboxes;         // number of boxes on the shelf
  int q;              // the q in q-gram
  int nstr;           // the number of stings compared
} Shelve;

// one shelve for all.
static Shelve shelve;

void print_box(Box *box){
  Rprintf("qgram: ");
  for (int i=0; i<shelve.q; i++ ){
    Rprintf("%d,",  box->intblocks[i]);
  }
  Rprintf("\n");
  Rprintf("count: ");
  for (int i=0; i<shelve.nstr; i++ ){
    Rprintf("%g",box->dblblocks[i]);
  }  
  Rprintf("\n");
}

void print_shelve(){
  Rprintf("Shelve q    :%d\n",shelve.q);
  Rprintf("Shelve nstr :%d\n",shelve.nstr);
  Rprintf("Shelve nbox :%d\n",shelve.nboxes);
  
}

static void init_shelve(int q, int nstr){
  shelve.q = q;
  shelve.nstr = nstr;
  shelve.nboxes = 0L;
//  for ( int i=0; i<MAXBOXES; i++ ){ 
//    shelve.box[i] = NULL; 
//  }
}

/* add box to shelve
 *
 * return values:
 * 0: okidoki
 * 1: no luck: out of memory or MAXBOXES exceeded
 *
 */
static int add_box(int nnodes){
Rprintf("adding box! %d\n",nnodes);
  int nb = shelve.nboxes;
  if ( nb < MAXBOXES ){
    shelve.box[nb] = new_box(nnodes, shelve.q, shelve.nstr);
    ++shelve.nboxes;
  } else {
    return 1;
  }
  return 0;
}

static void clear_shelve(){
  for ( int i = 0; i < shelve.nboxes; i++ ){
    free_box(shelve.box[i]);
  }
  shelve.nboxes=0L;
}

/* the allocator can store ints (q-grams), doubles (counts, per string) 
 * and qtree (nodes).
 */
typedef enum { uInt, Double, Qtree } type;


// cf. n1256.pdf (C99 std) sect 6.3.2.3 for pointer conversion.
static void *alloc(type t){

  if ( shelve.nboxes == 0L ){
    // TODO add check.
    add_box(MIN_BOX_SIZE);
  }

  Box *box = shelve.box[shelve.nboxes-1L];
  if ( box->nalloc == box->nnodes ){
    // add box such that storage size is doubled.
    if ( !add_box(2^(shelve.nboxes-1L) * MIN_BOX_SIZE) ){
      return NULL;
    }
    box = shelve.box[shelve.nboxes-1L];
  }
  
  void *x;
  switch ( t){
    case uInt:
      x = (void *) (box->intblocks + box->nalloc * shelve.q);
      break;
    case Double:
      x = (void *) (box->dblblocks + box->nalloc * shelve.nstr);
      break;
    case Qtree:
  //Rprintf("box->qtrblocks %d\n", box->qtrblocks);
  //Rprintf("box->nalloc %d\n", box->nalloc);
      x = (void *) (box->qtrblocks + box->nalloc);
      break;
    default:
      // TODO: raise hell
      break;
  }
  box->nalloc += 1;
  return x;
}


/* -- END OF ALLOCATOR -- */


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

// helper functions
static qtree *new_qtree(int q, int nstr){
  init_shelve(q, nstr);
  return NULL;
}

// TODO remove Q argument
static void free_qtree(){
  clear_shelve();
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
      free_qtree();
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
  qtree *Q = new_qtree(q, 2L);
 
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


  free_qtree();

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
  qtree *Q = new_qtree(q, 2);

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
  
  free_qtree();
  UNPROTECT(4);

  return(qcount);
}



