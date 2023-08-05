#include "com.h"

#define Index2d(i, j, dim2) ((i*dim2) + j)
typedef complex double cdouble;
#ifdef DEBUG
void printdata(int Nsource, int Nhalomax, int *Nhalo, int *sid, double *thetalst, double *philst, double *alphalst);
void printdatatofile(
    FILE *file, int Nsource, int Nhalomax, int *Nhalo, 
    int *sid, double *thetalst, double *philst, double *alphalst);
#endif
void get_batch_size(int data_size, int nproc, int **local_size_out, int **displs_out);
void readfile(
    FILE *input, int *Nsourceout, int**Nhaloout, 
    double **thetalst, double **philst, double **alphalst, int *NhaloMaxout, int **sidout);

cdouble Gamma(
    double theta, double phi, double alpha, 
    double delta, double A, double factort, double factorx);

void get_local_shear_profil(
    int *sid,
    double *thetalis, double *philis, double *alphalis,
    cdouble *profil, cdouble *profilconj,
    double A, 
    double delta,double epsilonh, 
    int Nsources, int* Nhalos, int NhaloMax
);

void distribute(
    int Nsource, int nproc, int rank, int root, int NhaloMax,
    int *sid, int *Nhalos, double *thetalst, double *philst, double *alphalst,
    double A, double delta, double epsilonh);

void printdatatofile(
    FILE *file, int Nsource, int Nhalomax, int *Nhalo, 
    int *sid, double *thetalst, double *philst, double *alphalst);