
#ifndef stringdist_stringset
#define stringdist_stringset

/* A structure for storing integer reps of strings */
typedef struct {
  // array of pointers to integer representation of strings, stored in data.
  unsigned int **string;
  //  array of string lengths.
  int *str_len;
  //  storage room for integer representation of strings
  unsigned int *data;
} Stringset;



#endif



