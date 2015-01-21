#ifndef SD_QTREE_H
#define SD_QTREE_H

#ifdef _OPENMP
#include <omp.h>
#endif

/* binary tree; dictionary of qgrams */

typedef struct qnode {
  unsigned int *qgram; // the q-gram.
  double *n;           // (vector of) counts.
  struct qnode *left;
  struct qnode *right;
} qtree;



#endif
