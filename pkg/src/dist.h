
#ifndef SD_DIST_H
#define SD_DIST_H

#include "qtree.h"
#include "dictionary.h"


double osa_dist(unsigned int *, int, unsigned int *, int, double *, double *);
dictionary *new_dictionary(unsigned int);
void free_dictionary(dictionary *);
double dl_dist(unsigned int *, int, unsigned int *, int, double *, dictionary *, double *);
double hamming_dist(unsigned int *, int, unsigned int *, int);
double lcs_dist(unsigned int *, int, unsigned int *, int, double *);
double lv_dist(unsigned int *, int, unsigned int *, int, double *, double *);
double osa_dist(unsigned int *, int, unsigned int *, int, double *, double *);
double jaro_winkler_dist(unsigned int *, int, unsigned int *, int, double, double, double *, double *);
qtree *new_qtree(int, int);
void free_qtree(void);
double qgram_dist(unsigned int *, int, unsigned int *, int, unsigned int, qtree **, int);
double soundex_dist(unsigned int *, int, unsigned int *, int, unsigned int *);
double running_cosine_dist(unsigned int *, int, unsigned int *, int, unsigned int, qtree **, double *);
#endif
