#ifndef MATRIX_H
#define MATRIX_H

#include <arb.h>
#include <doubly_linked_list.h>

typedef struct matrix{

  long rank,size;
  void **t;
}matrix;

matrix *initializeMatrix(long *dims, long rank);
matrix *matrixFromVector(cmplx *v);
matrix *matrixFromRealVector(double *v);
matrix *matrixFromMatrix(cmplx **v);
void freeMatrix(void *t1);
matrix *getMatrixElement(matrix *t,long *indxs);
void *getMatrixElement1(matrix *t,long indx);
long setMatrixElement1(matrix *t,long indx,matrix *t1);
long setMatrixElement(matrix *t,long *indxs,matrix *t1);
matrix *initializeMatrix0(cmplx cval);
cmplx *getMatrixValue(matrix *t,long *indxs);
int setMatrixValue(matrix *t,long *indxs,cmplx cval);
void *copyMatrix(void *t1);
void readMatrixFromStream(matrix **t, void *fs, stream_type st);
void *readMatrix(void *fs, stream_type st);
void writeMatrix(void *t,void *fs,stream_type st);
char *printMatrixDimensions(matrix *t);
void printMatrixDimensionsSub(char **outs,matrix *t);
long addMatrix(matrix *tdest,matrix *tnew);
matrix *diagonalMatrix(matrix *elem);
matrix *addMatrices(matrix *t11, matrix *t22);
long testSquareMatrix(matrix *m);
matrix *subtractMatrices(matrix *t11,matrix *t22);
long subtractMatrix(matrix *tdest,matrix *tnew);
long multiplyScalar(matrix *tdest,cmplx tnew);
cmplx innerProduct(matrix *t1, matrix *t2, long conj, long *success);
matrix *innerProductMul(matrix *t1, matrix *t2);
matrix *multiplyMatrices(matrix *t1, matrix *t2);
dl_list PermutationCycleDecomposition(long *indxs);
long TestEqualValuesL(long *indxs);
long TestEqualValuesLSub(long i1,long *indxs,long ni);
long PermutationParity(long *indxs);
long LeviCivita(long *indxs);
char *PrintPermutationCycles(long *indx);
long *getMatrixDimensions(matrix *m);
void getMatrixDimensionsSub(matrix *m,long *dims);
matrix *identityMatrix(long rank,long size);
matrix *matrixFromScalar(cmplx scal,long *dims);
matrix *numberMatrix(long rank,long size,cmplx v);
matrix *directSum(matrix *m11, matrix *m22);
matrix *directProduct(matrix *m11, matrix *m22);
void transposeMatrix(matrix *m, long recurse);
void conjugateMatrix(matrix *m);
matrix *matrixTrace(matrix *m);
cmplx **matrixToMatrix(matrix *v);
double **matrixToRealMatrix(matrix *v);
matrix *matrixDeterminant(matrix *m);
long invertMatrix(matrix *m);
matrix *matrixEigenvalues(matrix *m);
matrix *matrixEigenvectors(matrix *m);
void vectorizeMatrix(matrix *m);
void matricizeMatrix(matrix *m);
cmplx *getMatrixScalar(matrix *m);
matrix *augmentedMatrix(matrix *m, long cpy);
void freeAugmentedMatrix(matrix *am, long cpy);
matrix *augmentedMatrixLayers(matrix *m,long nlayers,long cpy);
void freeAugmentedMatrixLayers(matrix *am,long nlayers,long cpy);
matrix *initializeMatrixWithElements(matrix *elems,long *dims,long orank);
matrix *augmentedMatrixElements(matrix *m, long depth, long cpy);
void freeAugmentedMatrixElements(matrix *am, long depth, long cpy);
void swapElements(matrix *m,long i1,long i2);
void swapColumns(matrix *m,long i1,long i2);
long appendMatricesV(matrix *mdest,matrix *mnew);
long appendMatricesH(matrix *mdest,matrix *mnew);
long cropMatrix(long pos1,long pos2,matrix *m);
matrix *appendMatrix(matrix *mdest, matrix *mnew, long *success);
void enlargeMatrix(matrix *mdest,long newsize);
matrix *appendColumn(matrix *mdest,matrix *mnew,long *success);
matrix *insertColumn(matrix *mdest,matrix *mnew,long pos,long *success);
matrix *insertMatrix(matrix *mdest,matrix *mnew,long pos,long *success);
void removeNulls(matrix *mdest,long start,long end,long nnulls);
long removeElement(long pos,matrix *m);
long removeElements(long *pos,matrix *m);
long removeColumns(long *pos,matrix *m);

#endif
