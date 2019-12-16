#ifndef ARRAY_H
#define ARRAY_H

#include <arb.h>

typedef struct array{

  long size;
  void **table;
  void (*free_value)(void *);
  void *(*read_value)(void *,stream_type);
  void (*write_value)(void *,void *,stream_type);
  void *(*copy_value)(void *);
}array;

array initializeArray(long n, void (*freefunc)(void *),
                      void *(*readfunc)(void *,stream_type),
                      void (*writefunc)(void *, void *,stream_type),
                      void *(*copyfunc)(void *));
void freeArray(array *a, int v);
void removePosition(array *a,long pos);
void readArrayFromFile(array *a,char *fname);
void readArrayFromStream(array *a,void *fs,stream_type st);
void writeArrayToFile(array *a,long start,long end,char *fname);
void writeArrayToStream(array *a,long start,long end,void *fs,
                        stream_type st);
void resizeArray(array *a,long n);
long setElement(array *a,long i,void *newval);
void *getElement(array *a,long i);
array concatenateArrays(array *a1,array *a2);
array cropArray(array *a,long start,long end);
void freeElement(array *a,long i);
array copyArray(array *a0);
void appendArray(array *a1,array *a2);

#endif
