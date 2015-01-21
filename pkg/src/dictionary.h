
#ifndef SD_DICTIONARY_H
#define SD_DICTIONARY_H

/* Unsorted dictionary for dl distance */
typedef struct {
  // character
  unsigned int *key; 
  // number found
  unsigned int *value;
  // size of dictionary
  unsigned int length;
} dictionary;


#endif
