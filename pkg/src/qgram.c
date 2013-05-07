/* This is a very simple implementation of the q-gram distance.
 * q-grams for strings s and t are pushed onto an unsorted list where they
 * are counted as well. 
*/

#include <string.h>
#include <R.h>
#include <Rdefines.h>


/* Some workspace stuff */
typedef struct{
    unsigned int *A;
    unsigned int nrec;
    unsigned int recsize;
} workspace;

static workspace *new_workspace(unsigned int nrec, unsigned int recsize){
  workspace *ws   = (workspace *) malloc(sizeof(workspace));
  if (ws == NULL){ 
      return NULL;
  }
  ws->A = (unsigned int *) calloc(sizeof(int),nrec*recsize);
  if (ws->A == NULL){ 
      return NULL;
  }
  ws->nrec            = nrec;
  ws->recsize         = recsize;
  return ws;   
}

static void free_workspace(workspace *ws){
  free(ws->A);
  free(ws);
}

static int double_workspace(workspace *ws){
  unsigned int *A = (unsigned int *) calloc(sizeof(int), 2 * ws->nrec * ws->recsize );
  if (A == NULL){
    return(0);
  }
  unsigned int *a = ws->A;
  ws->A = A;
  memcpy(ws->A, a, sizeof(int) * ws->nrec * ws->recsize);
  free(a);
  ws->nrec    *= 2;
  return 1;
}

static void print_workspace(workspace *ws){
  for ( int i=0; i<ws->nrec; ++i){
    for (int j=0; j<ws->recsize; ++j ){
      Rprintf("%d, ",ws->A[i*ws->recsize + j]);
    }
    Rprintf("\n");
  }
}

static void reset_workspace(workspace *ws){
  memset(ws->A, 0, sizeof(int) * ws->nrec * ws->recsize);
}

/* Actual work starts here */

static int qgram_equals(unsigned int *x, unsigned int *y, int q){
  int i = 0;
  while ( i < q & x[i] == y[i] ) ++i;
  return ( i == q ) ? 1 : 0;
}


static int add_qgram(
        unsigned int *s, 
        int q, 
        int offset, 
        workspace *work,
        unsigned int nfilled
){
  size_t qg_size = sizeof(unsigned int) * q;
  int i = 0; 
  int i_pos = 0;

  // Find out if qgram already exists in workspace
  while ( i < nfilled & !qgram_equals(s, work->A + i_pos, q) ){
      ++i;
      i_pos += q + 2;
  }

  // if not, add qgram
  if ( i == nfilled ){
    memcpy(work->A + i_pos, s, qg_size);
    nfilled += 1;
  }

  // increase counter for qgram.
  work->A[i_pos + q + offset] += 1;

  // double workspace size if workspace is full
  if (nfilled == work->nrec){
    if ( !double_workspace(work) ){
      free_workspace(work);
      error("Could not allocate memory");
    }
  }

  // return number of stored qgrams in workspace.
  return nfilled;
}


// Compute qgram distance using L1-distance between sparse vectors
// This is a naive implementation, but its fast when comparing short strings.
static int qgram(
      unsigned int *s, 
      unsigned int *t, 
      unsigned int x, 
      unsigned int y, 
      unsigned int q, 
      workspace *work
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

  if ( q == 0 && x + y > 0 ) return -1;
  // edge case, comparing two empty strings with q=0 gives 0.
  if ( x == 0 && y == 0 && q == 0 ) return 0;
  int nfill = 0;
  // store and count qgrams.
  for (int i=0; i < x-q+1; ++i){
    nfill = add_qgram(s+i, q, 0, work, nfill);
  }
  for (int i=0; i < y-q+1; ++i){
    nfill = add_qgram(t+i, q, 1, work, nfill);
  } 
  int L1 = 0;
  int ncols = q + 2;
  for ( int j=0; j < nfill; j++ ){
    L1 +=  abs(work->A[j*ncols + q] - work->A[j*ncols + q + 1]);
  }
  reset_workspace(work);
  return L1;
}



SEXP R_qgram(SEXP a, SEXP b, SEXP qq){
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

  // Enough room to store 10k unique q-grams and counters (this is a startvalue; worspace is expanded as needed).
  workspace *work = new_workspace(10000,q+2);
  if (work == NULL){
    error("Could not allocate enough memory");
  }
  for ( k=0; k < nt; ++k ){
    i = k % na;
    j = k % nb;
    if (INTEGER(VECTOR_ELT(a,i))[0] == NA_INTEGER || INTEGER(VECTOR_ELT(b,j))[0] == NA_INTEGER){
      y[k] = NA_REAL;
      continue;
    }
    y[k] = qgram(
        INTEGER(VECTOR_ELT(a,i)),
        INTEGER(VECTOR_ELT(b,j)),
        length(VECTOR_ELT(a,i)),
        length(VECTOR_ELT(b,j)),
        q,
        work
    );
  }
  free(work);
  UNPROTECT(3);
  return yy;
}

