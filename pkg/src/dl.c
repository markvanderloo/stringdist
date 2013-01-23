/* This code is gratefully copied from Nick Logan's githup repository
 * https://github.com/ugexe/Text--Levenshtein--Damerau--XS/blob/master/damerau-int.c
 * 
 *
 * Changes/additions to original code:
 * - Added R.h, Rdefines.h inclusion
 * - Added R interface function
 * - Added edit weights (function is now of type double)
 *
 * MvdL, 2013.01.23
 * mark.vanderloo@gmail.com
 */

/* ugexe@cpan.org (Nick Logan)    */

#include <R.h>
#include <Rdefines.h>
#include "utils.h"



/* ugexe@cpan.org (Nick Logan)    */


/* Our unsorted dictionary linked list.   */
/* Note we use character ints, not chars. */

struct dictionary{
  unsigned int key;
  unsigned int value;
  struct dictionary* next;
};
typedef struct dictionary item;


static __inline item* push(unsigned int key,item* curr){
  item* head;
  head = malloc(sizeof(item));   
  head->key = key;
  head->value = 0;
  head->next = curr;
  return head;
}


static __inline item* find(item* head,unsigned int key){
  item* iterator = head;
  while(iterator){
    if(iterator->key == key){
      return iterator;
    }
    iterator = iterator->next;
  }
 
  return NULL;
}

static __inline item* uniquePush(item* head,unsigned int key){
  item* iterator = head;

  while(iterator){
    if(iterator->key == key){
      return head;
    }
    iterator = iterator->next;
  }
 
  return push(key,head); 
}

static void dict_free(item* head){
  item* iterator = head;
  while(iterator){
    item* temp = iterator;
    iterator = iterator->next;
    free(temp);
  }

  head = NULL;
}

/* End of Dictionary Stuff */



 
/* All calculations/work are done here */

static double distance(unsigned int src[],unsigned int tgt[],unsigned int x,unsigned int y,double *weight,double maxDistance){
  item *head = NULL;
  unsigned int swapCount, targetCharCount,i,j;
  double delScore, insScore, subScore, swapScore;
  double *scores = malloc( (x + 2) * (y + 2) * sizeof(double) );
  unsigned int score_ceil = x + y;

  /* intialize matrix start values */
  scores[0] = score_ceil;  
  scores[1 * (y + 2) + 0] = score_ceil;
  scores[0 * (y + 2) + 1] = score_ceil;
  scores[1 * (y + 2) + 1] = 0;
  head = uniquePush(uniquePush(head,src[0]),tgt[0]);
  /* work loops    */
  /* i = src index */
  /* j = tgt index */
  for(i=1;i<=x;i++){ 
    head = uniquePush(head,src[i]);
    scores[(i+1) * (y + 2) + 1] = i;
    scores[(i+1) * (y + 2) + 0] = score_ceil;

    swapCount = 0;
    for(j=1;j<=y;j++){
      if(i == 1) {
          head = uniquePush(head,tgt[j]);
          scores[1 * (y + 2) + (j + 1)] = j;
          scores[0 * (y + 2) + (j + 1)] = score_ceil;
      }

      targetCharCount = find(head,tgt[j-1])->value;
      swapScore = scores[targetCharCount * (y + 2) + swapCount] + i - targetCharCount - 2 + j - swapCount + weight[3];

      if(src[i-1] != tgt[j-1]){
        subScore = scores[i * (y + 2) + j] + weight[2];
        insScore = scores[(i+1) * (y + 2) + j] + weight[1];
        delScore = scores[i * (y + 2) + (j + 1)] + weight[0];
        scores[(i+1) * (y + 2) + (j + 1)] = min2(swapScore, min3(delScore, insScore, subScore));
      } else {
        swapCount = j;
        scores[(i+1) * (y + 2) + (j + 1)] = min2(scores[i * (y + 2) + j], swapScore);
      } 
    }

    /* We will return -1 here if the */
    /* current score > maxDistance   */
    if(maxDistance != 0 && maxDistance < scores[(i+1) * (y + 2) + (y+1)]) {
      dict_free(head);
      free(scores);
      return -1;
    }

    
    find(head,src[i-1])->value = i;
  }

  {
  double score = scores[(x+1) * (y + 2) + (y + 1)];
  dict_free(head);
  free(scores);
  return score;
  }
}

// -- interface with R 

// We need to convert R's char type (1 word) to unsigned int (4 words)
static void char2uint(unsigned int *destination, const char *source){
   int i = 0;
   while(source[i]){
      destination[i] = (unsigned int) source[i];
      ++i;
   }
}

SEXP R_dl(SEXP a, SEXP b, SEXP ncharA, SEXP ncharB, SEXP weight, SEXP maxDistance, SEXP maxchar){
   PROTECT(a);
   PROTECT(b);
   PROTECT(ncharA);
   PROTECT(ncharB);
   PROTECT(maxDistance);
   PROTECT(maxchar);
   PROTECT(weight);
   
   int i, j, k;
   int na = length(a);
   int nb = length(b);
   int nt = (na > nb) ? na : nb;
   double maxDist = REAL(maxDistance)[0];
   double *w = REAL(weight);

   SEXP yy; 
   PROTECT(yy = allocVector(REALSXP, nt));
   double *y = REAL(yy);

   unsigned int *aa = (unsigned int *) calloc(INTEGER(maxchar)[0], sizeof(unsigned int));
   unsigned int *bb = (unsigned int *) calloc(INTEGER(maxchar)[0], sizeof(unsigned int));

   for ( k=0; k < nt; ++k ){
      i = k % na;
      j = k % nb;
      char2uint(aa, CHAR(STRING_ELT(a,i)));
      char2uint(bb, CHAR(STRING_ELT(b,j)));
      y[k] = distance(
         aa,
         bb,
         (unsigned int)  INTEGER(ncharA)[i],
         (unsigned int)  INTEGER(ncharB)[j],
         w,
         maxDist
      );
   }

   free(aa);
   free(bb);
   UNPROTECT(8);
   return yy;
} 

