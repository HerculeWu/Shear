// common includes
#ifndef COMM_H
#define COMM_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
// #include <omp.h>
#include <math.h>
#include <complex.h>
// common used macros and type definition
#define I _Complex_I
#define EXIT(n) MPI_Finalize(); exit(n)
#define memCheck(ptr) if(ptr == NULL) {printf("memory allocation fail\n");EXIT(100);}
#define fileCheck(file) if (file==NULL) {printf("file fail\n"); EXIT(200);}
#define malloc1d(type, size) (type *)malloc(size*sizeof(type))
#define Index2d(i, j, dim2) ((i*dim2) + j)
#define Index3d(i, j, k, dim2, dim3) ((i*dim2*dim3) + (j*dim2) + k)
#define malloc2d(type, dim1, dim2) (type *)malloc(dim1*dim2*sizeof(type))
#define malloc3d(type, dim1, dim2, dim3) (type *)malloc(dim1*dim2*dim3*sizeof(type))
#define readcheck(Nget, Nshould) if (Nget != Nshould) {printf("read fail!\n"); exit(110);}
// typedef complex double cdouble;
// typedef enum{True = 1, False = 0} bool;
#endif