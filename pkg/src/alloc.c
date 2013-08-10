
#include <stdio.h>
#include <stdlib.h>

typedef long long int qtree;

typedef enum { Qgram, Count, Qtree } type;


#define MAXBOXES 20
// number of nodes to initially allocate space for.
#define MIN_BOX_SIZE (1<<10)  

// A Box can store a number of nodes
typedef struct {
  
  int nnodes; // number of nodes
  int nalloc; // number of nodes allocated

  // Room in the box.
  int     *intblocks; // store qgrams
  double  *dblblocks; // store qgram-counts
  qtree   *qtrblocks; // store nodes
} Box;

// A shelve can store up to MAXBOXES Boxes.
typedef struct {
  Box **box[MAXBOXES];
  int nboxes;         // number of boxes on the shelf
  int q;              // the q in q-gram
  int nstr;           // the number of stings compared
} Shelve;


static Box *new_box(int nnodes, int q, int nstr){
  Box *b = malloc(sizeof(Box));
// TODO: raise hell ico malloc failure.
  b->nnodes     = nnodes;
  b->nalloc     = 0L;
  b->intblocks  = (int *) malloc(sizeof(int) * nnodes * q);
  b->dblblocks  = malloc(sizeof(double) * nnodes * nstr);
  return b;
}

static void free_box(Box *box){
  // empty box
  free( box->intblocks);
  free( box->dblblocks);
  free( box->qtrblocks);
  // throw box
  free( box);
}

// one shelve for all.
static Shelve shelve;

static void init_shelve(int q, int nstr){
  
  shelve.nboxes = 0L;
  for ( int i=0; i<MAXBOXES; i++ ){ 
    shelve.box[i] = NULL; 
  }
}

static void add_box(int nnodes){
  int nb = shelve.nboxes;
  if ( nb + 1 < MAXBOXES ){
    shelve.box[nb + 1] = new_box(nnodes, shelve.q, shelve.nstr);
    ++shelve.nboxes;
  } else {
    // TODO: raise hell.
  }
}

static void clear_shelve(){
  for ( int i = 0; i < shelve.nboxes; i++ ){
    free( shelve.box[i] );
  }
  shelve.nboxes=0L;
}


// cf. n1256.pdf (C99 std) sect 6.3.2.3
static void *alloc(type t){

  if ( shelve.nboxes == 0L ){
    add_box(MIN_BOX_SIZE)
  }

  int ibox = shelve.nboxes;
  Box *box = shelve.box[ibox];
  if ( box->nalloc == box->nnodes ){
    // add box such that storage size is doubled.
    add_box(2^(shelve.nboxes-1L) * MIN_BOX_SIZE);
    ibox = shelve.nboxes;
    box = shelve.box[ibox];
  }
  
  void *x;
  switch ( t){
    case Qgram:
      x = (void *) box->intblocks + box.nnodes * shelve.q;
      break;
    case Count:
      x = (void *) box->dblblocks + box.nnodes * shelve.nstr;
      break;
    case Qtree:
      x = (void *) box->qtrblocks + box.nnodes;
      break;
    default:
      // TODO: raise hell
      break;
  }

  return x;
}



int main(){

  int x[1] = {7};
  void *p = (void *) x;
  int *X = (int *) p;
  printf("hebbewe'm? %d\n",X[0]);
  return 0L;
}

