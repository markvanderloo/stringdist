

#include <R.h>
#include <Rdefines.h>
#include "utils.h"
#include <string.h>


static inline int max(int x, int y){
  int m = (x < y) ? y : x;
  return m;
}

static inline int min(int x, int y){
  int m = (x < y) ? x : y;
  return m; 
}


/* First match of a in b[] 
 * Returns -1 if no match is found
 * Parameter 'guard; indicates which elements of b have been matched before to avoid
 * matching two instances of the same character to the same position in b (which we treat read-only).
 */
static int match_int(unsigned int a, unsigned int *b, int *guard, int width){

  int i = 0;
  while ( 
      b[i] && 
      ( i < width ) && 
      ( b[i] != a || (b[i] == a && guard[i])) 
  ){
    ++i;
  }
  if ( b[i] == a ){
    guard[i] = 1;
    return i;
  } 
  return -1;
}


/* jaro distance (see http://en.wikipedia.org/wiki/Jaro%E2%80%93Winkler_distance).
 *
 * a    : string (in uint rep)
 * b    : string (in uint rep)
 * x    : length of a (in uints)
 * y    : length of b (in uints)
 * work : workspace, minimally of length max(x,y)
 *
 */
double jaro( unsigned int *a, 
             unsigned int *b,
             int x,
             int y,
             int *work
        ){

  // swap arguments if necessary, so we always loop over the shortest string
  if ( x > y ){
    unsigned int *c;
    int z;
    c = b;
    b = a;
    a = c;
    z = y;
    y = x;
    x = z;
  }
//printf("(x,y): (%d,%d)\n",x,y);
  
  // edge case
  if ( x==0 && y == 0 ) return 0;

  // max transposition distance
  int M = max(max(x,y)/2 - 1,0);
  // transposition counter
  double t = 0.0;
  // number of matches 
  double m = 0.0;
  
  int left, right, J, jmax=0;
  
  for ( int i=0; i < x; ++i ){
    left  = max(0, i-M);

    if ( left >= y ){
      J = -1;
    } else {
      right = min(y, i+M);
      J =  match_int(a[i], b + left, work + left, right - left);
    }
    if ( J >= 0 ){
      ++m;
      t += (J + left < jmax ) ? 1 : 0; 
      jmax = max(jmax, J + left);
    }
// printf("Left = %2d right = %2d J = %2d, J+LEFT: %2d \n",left,right,J, J+left);
  }
// printf("M = %d, m=%g, t=%g\n",M, m, t);
  double d;
  if ( m < 1 ){
    d = 1.0;
  } else {
    d = 1.0 - (1.0/3.0)*(m/x + m/y + (m-t)/m);
  }
  memset(work,0,sizeof(int) * y);
  return d;
}



/* R interface                                                           */


SEXP R_jaro(SEXP a, SEXP b){
  PROTECT(a);
  PROTECT(b);

  // find the length of longest strings
  int max_a = max_length(a);
  int max_b = max_length(b);
  int max_char = max(max_a,max_b);
  // recycling parameters
  int na = length(a);
  int nb = length(b);
  int nt = max(na,nb);
  
  int *s, *t;
  int length_a, length_b;

  // workspace for worker function
  int *work = (int *) calloc( sizeof(int), max_char );

  // output variable
  SEXP yy;
  PROTECT(yy = allocVector(REALSXP,nt));
  double *y = REAL(yy);

  // compute distances, skipping NA's
  int i,j;
  for ( int k=0; k < nt; ++k ){
    i = k % na;
    j = k % nb;

    length_a = length(VECTOR_ELT(a,i));
    length_b = length(VECTOR_ELT(b,j));
    s = INTEGER(VECTOR_ELT(a,i));
    t = INTEGER(VECTOR_ELT(b,j));
    if ( s[0] == NA_INTEGER || t[0] == NA_INTEGER){
      y[k] = NA_REAL;
      continue;
    } else {
      y[k] = jaro(s, t, length_a, length_b, work);
    }
    //memset(work, 0 , sizeof(int) * max_char);
  }

  
  UNPROTECT(3);
  free(work);
  return yy;

}







/*

#include <stdlib.h>
#include <stdio.h>
void main(){
  unsigned int a1[6] = {4,23,1,25,14,5};       // DWAYNE
  unsigned int b1[5] = {4,21,1,14,5};          // DUANE
  unsigned int a2[8] = {4,9,3,11,19,15,14,24}; // DICKSONX
  unsigned int b2[5] = {4,9,24,15,14};         // DIXON
  unsigned int a3[6] = {13,1,18,20,8,1};       // MARTHA
  unsigned int b3[6] = {13,1,18,8,20,1};       // MARHTA
  unsigned int a4[5] = {1,1,16,10,5};       // MARTHA
  unsigned int b4[1] = {1};       // MARHTA
  int *work = (int *) calloc(sizeof(int),200);
  printf("distance DWAYNE vs DUANE: %g\n", 1-jaro(a1,b1,6,5,work) );
  printf("distance DICKSONX vs DIXON: %g\n", 1-jaro(a2,b2,8,5,work) );
  printf("distance MARTHA vs MARHTA: %g\n", 1-jaro(a3,b3,6,6,work) );
  printf("distance AAPJE vs  A: %g\n",1-jaro(a4,b4,5,1,work));
  printf("distance A vs  AAPJE: %g\n",1-jaro(b4,a4,1,5,work));
  free(work);
}
*/




