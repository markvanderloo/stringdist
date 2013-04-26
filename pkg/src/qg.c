
//#include<stdlib.h>
#include <string.h>
#include <R.h>
#include <Rdefines.h>

//#include<stdio.h>

static int qgram_equals(unsigned int *x, unsigned int *y, int q){
  int i=0;
  while ( i < q & x[i] == y[i] ) ++i;
  return ( i == q ) ? 1 : 0;
}

static int add_qgram(unsigned int *s, int q, int offset, unsigned int *work, unsigned int nfilled){
  size_t qg_size = sizeof(unsigned int) * q;
  int i = 0; 
  int i_pos = 0;
  // Find out if qgram already exists in workspace
  while ( i < nfilled & !qgram_equals(s, work + i_pos, q) ){
      ++i;
      i_pos += q + 2;
  }
  // if not, add qgram
  if ( i == nfilled ){
    memcpy(work + i_pos, s, qg_size);
    nfilled += 1;
  }
  // increase counter for qgram.
  work[i_pos + q + offset] += 1;
  // return number of stored qgrams in workspace.
  return nfilled;
}


// Compute qgram distance using L1-distance between sparse vectors
// This is the naive implementation, but its fast when comparing short strings.
int qgram(
      unsigned int *s, 
      unsigned int *t, 
      unsigned int x, 
      unsigned int y, 
      unsigned int q, 
      unsigned int *work
    ){

  int nfill = 0;
  // store and count qgrams.
  // TODO: add workspace doubler.
  for (int i=0; i < x-q+1; ++i){
    nfill = add_qgram(s+i, q, 0, work, nfill);
  }
  for (int i=0; i < y-q+1; ++i){
    nfill = add_qgram(t+i, q, 1, work, nfill);
  } 
  int L1 = 0;
  int nrows = q + 2;
  for ( int j=0; j<nfill; j++ ){
    L1 = L1 + abs(work[j*nrows + q] - work[j*nrows + q + 1]);
  }
  return L1;
}



SEXP R_qgram(SEXP a, SEXP b, SEXP qq){
  PROTECT(a);
  PROTECT(b);
   
  int i, j, k;
  int na = length(a);
  int nb = length(b);
  int nt = (na > nb) ? na : nb;
  int q = *INTEGER(qq);

  SEXP yy; 
  PROTECT(yy = allocVector(REALSXP, nt));
  double *y = REAL(yy);

  // ask for 1M of ints
  unsigned int *work = (unsigned int *) malloc(sizeof(unsigned int)*1024*1024);
  for ( k=0; k < nt; ++k ){
    i = k % na;
    j = k % nb;
    if (INTEGER(VECTOR_ELT(a,i))[0] == NA_INTEGER || INTEGER(VECTOR_ELT(b,j))[0] == NA_INTEGER){
      y[k] = NA_INTEGER;
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

/*
void main(){
  unsigned int s[4] = {0,1,2,3};
  unsigned int t[5] = {2,3,4,5,6};
  unsigned int *w = malloc(sizeof(int)*1000);


  int d = qg_dist_def(s,t,4,5,2,w,10);
  printf("%d\n",qgram_equals(w,w,1));

  printf("%d\n",d); 
  free(w);


}

static void print_work(unsigned int *w,int q, int nrow){
  for ( int j=0; j<nrow; ++j){
    for( int i=0; i < q+2; ++i ){
      printf("%d ",w[j*(q+2) + i]);
    }
    printf("\n");
  }

}

*/
