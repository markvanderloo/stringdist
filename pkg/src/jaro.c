  
#include <stdlib.h>
#include <stdio.h>


static inline int max(int x, int y){
  int m = (x < y) ? y : x;
  return m;
}

static inline int min(int x, int y){
  int m = (x < y) ? x : y;
  return m; 
}


// First match of a in b[]
// Returns -1 if no match is found
// Parameter 'guard; indicates which elements of b have been matched before to avoid
// matching two instances of the same character to the same position in b (which we treat read-only).
int match(unsigned int a, unsigned int *b, int *guard, int width){

  int i = 0;
  while ( 
      b[i] &&       
      ( i <= width ) && 
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

double jaro( unsigned int *a, 
             unsigned int *b,
             int x,
             int y
        ){

  // max transposition distance
  int M = max(x,y)/2 - 1;
  // transposition counter
  double t = 0.0;
  // number of matches 
  double m = 0.0;
  // workspace to store guard against double matches
  int *guard = (int *) calloc(sizeof(int), 200);
  int left, right;

  int J;

  for ( int i=0; i < x; ++i ){
    left  = max(0, i-M);
    right = min(y, i+M);
    J =  match(a[i], b + left, guard + left, right - left);
    if ( J >= 0 ){
      ++m;
      t += (J + left > i ) ? 1 : 0; 
    } 
     
    
  }

  double d;
  if ( m < 1 ){
    d = 0;
  } else {
    d = (1.0/3.0)*(m/x + m/y + (m-t)/m);
  }
  free(guard);
  return d;
}

/*
void main(){
  unsigned int a1[6] = {4,23,1,25,14,5};       // DWAYNE
  unsigned int b1[5] = {4,21,1,14,5};          // DUANE
  unsigned int a2[8] = {4,9,3,11,19,15,14,24}; // DICKSONX
  unsigned int b2[5] = {4,9,24,15,14};         // DIXON
  unsigned int a3[6] = {13,1,18,20,8,1};       // MARTHA
  unsigned int b3[6] = {13,1,18,8,20,1};       // MARHTA
  printf("distance DWAYNE vs DUANE: %g\n", jaro(a1,b1,6,5) );
  printf("distance DICKSONX vs DIXON: %g\n", jaro(a2,b2,8,5) );
  printf("distance MARTHA vs MARHTA: %g\n", jaro(a3,b3,6,6) );

}
*/





