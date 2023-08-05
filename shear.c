#include "shear.h"

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, nproc, root=0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    int *Nhalos, *Nhalo_local, *sid, *sid_local;
    double *thetalst, *philst, *alphalst;
    double delta, A, epsilonh;
    int Nsource, NhaloMax;
    if (rank==0) {
        if (argc!=5) {
            printf("Usage: ./shear [input data] delta A epsilonh\n");
            printf("Data format:\n");
            printf("First line: number of sources, maximal number of involved halos\n");
            printf("Follow lines:\n");
            printf("1&2: source id,  number of involved halos\n");
            printf("follow columns: epsilon_h theta phi alpha, repeat");
            printf("Note, the maximal number of involved halos can be greater than the actual value\n");
            MPI_Finalize();
            return 0;
        }
        delta = atof(argv[2]);
        A = atof(argv[3]);
        epsilonh = atof(argv[4]);
        FILE *file = fopen(argv[1], "r");
        int status;
        readfile(file, &Nsource, &Nhalos, &thetalst, &philst, &alphalst, &NhaloMax, &sid);
        fclose(file);
#ifdef DEBUG
        printf("readed data:\n");
        printdata(Nsource, NhaloMax, Nhalos, sid, thetalst, philst, alphalst);
#endif
    }
    MPI_Bcast(&A, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(&delta, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(&epsilonh, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(&Nsource, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&NhaloMax, 1, MPI_INT, root, MPI_COMM_WORLD);
#ifdef DBUG
    printf("rank %d, A %lg, delta %lg, epsilonh %lg, Nsource %d, NhaloMax %d\n",
        rank, A, delta, epsilonh, Nsource, NhaloMax
    );
#endif
    distribute(Nsource, nproc, rank, root, NhaloMax, sid, Nhalos, thetalst, philst, alphalst, A, delta, epsilonh);

    if (rank == 0) {
        free(thetalst); free(philst); free(alphalst);
    }

    MPI_Finalize();
    return 0;
}

#ifdef DEBUG
void printdata(int Nsource, int Nhalomax, int *Nhalo, int *sid, double *thetalst, double *philst, double *alphalst) {
    for (int i = 0; i < Nsource; i++) {
        printf("%d %d", sid[i], Nhalo[i]);
        for (int j=0; j < Nhalomax; j++) {
            int k = Index2d(i, j, Nhalomax);
            printf(" %lg %lg %lg", thetalst[k], philst[k], alphalst[k]);
        }
        printf("\n");
    }
}

void printdatatofile(
    FILE *file, int Nsource, int Nhalomax, int *Nhalo, 
    int *sid, double *thetalst, double *philst, double *alphalst) {
    for (int i = 0; i < Nsource; i++) {
        fprintf(file, "%d %d", sid[i], Nhalo[i]);
        for (int j=0; j < Nhalomax; j++) {
            int k = Index2d(i, j, Nhalomax);
            fprintf(file, " %lg %lg %lg", thetalst[k], philst[k], alphalst[k]);
        }
        fprintf(file, "\n");
    }
}
#endif

void readfile(
    FILE *input, int *Nsourceout, int**Nhaloout, 
    double **thetalst, double **philst, double **alphalst, int *NhaloMaxout, int **sidout){
    int Nsource, Nhalomax;
    int Nget;
    int i, j, k;
    Nget = fscanf(input, "%d %d\n", &Nsource, &Nhalomax);
    readcheck(Nget, 2);
    *Nsourceout = Nsource;
    double *thetas = malloc1d(double, Nsource*Nhalomax);
    double *phis = malloc1d(double, Nsource*Nhalomax);
    double *alphas = malloc1d(double, Nsource*Nhalomax);
    memCheck(thetas); memCheck(phis); memCheck(alphas);
    int *Nhalos = malloc1d(int, Nsource);
    int *sid = malloc1d(int, Nsource);
    memCheck(Nhalos); memCheck(sid);
    for (i = 0; i < Nsource; i++) {
        Nget = fscanf(input, "%d %d", &sid[i], &Nhalos[i]);
        readcheck(Nget, 2);
        for (j = 0; j < Nhalomax-1; j++) {
            k = Index2d(i, j, Nhalomax);
            Nget = fscanf(input, " %lg %lg %lg", &thetas[k], &phis[k], &alphas[k]);
            // printf("%lg %lg %lg\n", thetas[k], phis[k], alphas[k]);
        }
        k = Index2d(i, Nhalomax-1, Nhalomax);
        Nget = fscanf(input, " %lg %lg %lg\n", &thetas[k], &phis[k], &alphas[k]);
        readcheck(Nget, 3);
        // printf("%d %d\r", sid[i], Nhalos[i]);
        // printf("\n");
    }
    *Nhaloout = Nhalos;
    *thetalst = thetas;
    *philst = phis;
    *alphalst = alphas;
    *NhaloMaxout = Nhalomax;
    *sidout = sid;
}

void distribute(
    /*
    Distribute to different mpi process
    Note every mpi write output file
    use gather script to combine them
    an internal gather is not needed here
    */
    int Nsource, int nproc, int rank, int root, int NhaloMax,
    int *sid, int *Nhalos, double *thetalst, double *philst, double *alphalst,
    double A, double delta, double epsilonh){
    int *local_size = malloc1d(int, nproc);
    int *displs = malloc1d(int, nproc);
    // data size for each source
    int *local_size_halo = malloc1d(int, nproc);
    int *displs_halo = malloc1d(int, nproc);
    if (rank == root) {
        get_batch_size(Nsource, nproc, &local_size, &displs);
    }
    MPI_Bcast(local_size, nproc, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(displs, nproc, MPI_INT, root, MPI_COMM_WORLD);
    for (int i = 0; i < nproc; i++) {
        local_size_halo[i] = local_size[i]*NhaloMax;
        displs_halo[i] = displs[i]*NhaloMax;
    }
    // local batch
    double *thetalst_local = malloc1d(double, local_size[rank]*NhaloMax);
    double *philst_local = malloc1d(double, local_size[rank]*NhaloMax);
    double *alphalst_local = malloc1d(double, local_size[rank]*NhaloMax);
    int *Nhalos_local = malloc1d(int, local_size[rank]);
    int *sid_local = malloc1d(int, local_size[rank]);
    // scatter the data
    if (rank==root) {
        MPI_Scatterv(
            Nhalos, local_size, displs, MPI_INT, 
            Nhalos_local, local_size[rank], MPI_INT, root, MPI_COMM_WORLD);
        MPI_Scatterv(
            sid, local_size, displs, MPI_INT, 
            sid_local, local_size[rank], MPI_INT, root, MPI_COMM_WORLD);
        MPI_Scatterv(
            thetalst, local_size_halo, displs_halo, MPI_DOUBLE,
            thetalst_local, local_size_halo[rank], MPI_DOUBLE, root, MPI_COMM_WORLD);
        MPI_Scatterv(
            philst, local_size_halo, displs_halo,
            MPI_DOUBLE, philst_local, local_size_halo[rank], MPI_DOUBLE, root, MPI_COMM_WORLD);
        MPI_Scatterv(
            alphalst, local_size_halo, displs_halo, MPI_DOUBLE,
            alphalst_local, local_size_halo[rank], MPI_DOUBLE, root, MPI_COMM_WORLD);
    } else {
        MPI_Scatterv(
            NULL, NULL, NULL, MPI_DATATYPE_NULL,
            Nhalos_local, local_size[rank], MPI_INT, root, MPI_COMM_WORLD);
        MPI_Scatterv(
            NULL, NULL, NULL, MPI_DATATYPE_NULL,
            sid_local, local_size[rank], MPI_INT, root, MPI_COMM_WORLD);
        MPI_Scatterv(
            NULL, NULL, NULL, MPI_DATATYPE_NULL,
            thetalst_local, local_size_halo[rank], MPI_DOUBLE, root, MPI_COMM_WORLD);
        MPI_Scatterv(
            NULL, NULL, NULL, MPI_DATATYPE_NULL,
            philst_local, local_size_halo[rank], MPI_DOUBLE, root, MPI_COMM_WORLD);
        MPI_Scatterv(
            NULL, NULL, NULL, MPI_DATATYPE_NULL,
            alphalst_local, local_size_halo[rank], MPI_DOUBLE, root, MPI_COMM_WORLD);
    }
#ifdef DEBUG
    char rankdata[20];
    sprintf(rankdata, "rank_%d.txt", rank);
    FILE *rankdatafile = fopen(rankdata, "w+");
    printdatatofile(
        rankdatafile,
        local_size[rank], NhaloMax, Nhalos_local, sid_local, thetalst_local, philst_local, alphalst_local);
    fclose(rankdatafile);
#endif
    cdouble *profil = malloc1d(cdouble, local_size[rank]);
    cdouble *profilconj = malloc1d(cdouble, local_size[rank]);
    memCheck(profil); memCheck(profilconj);
    get_local_shear_profil(
        sid_local, thetalst_local, philst_local, alphalst_local,
        profil, profilconj, A, delta, epsilonh,
        local_size[rank], Nhalos_local, NhaloMax);
    char outname[30];
    sprintf(outname, "out_%d.txt", rank);
    FILE *outfile = fopen(outname, "w+");
    fileCheck(outfile);
    // fprintf(outfile, "# sid profile[real, imag] profileconj[real, imag] Wtheta Ptheta\n");
    for (int i = 0; i < local_size[rank]; i++) {
        fprintf(outfile, "%d %.15lg %.15lg %.15lg %.15lg\n", 
                sid_local[i],
                creal(profil[i]), cimag(profil[i]),
                creal(profilconj[i]), cimag(profilconj[i])
                );
    }
    fclose(outfile);

    free(displs); free(local_size); free(local_size_halo); free(displs_halo);
    free(thetalst_local); free(alphalst_local); free(philst_local);
    free(Nhalos_local); free(sid_local);
}

// get shear profil in each rank
void get_local_shear_profil(
    int *sid,
    double *thetalis, double *philis, double *alphalis,
    cdouble *profil, cdouble *profilconj,
    double A, 
    double delta,double epsilonh, 
    int Nsources, int* Nhalos, int NhaloMax
) {
    double factort = (delta-2.0)*(4.-(2.*delta)+(delta*delta))/(delta*delta*(delta-4.));
    double factorx = 4.*(2.-delta)*(1.-delta)/(delta*delta*(delta-4.));
    factort = factort*epsilonh*0.5;
    factorx = factorx*epsilonh*0.5;
    int index, s, h;
    double theta, phi;
    
    for (s = 0; s < Nsources; s++)
    {
        profil[s] = 0.0 + 0.0*I;
        profilconj[s] = 0.0+0.0*I;
        for (h = 0; h < NhaloMax; h++)
        {
            if (h < Nhalos[s]){
                index = Index2d(s, h, NhaloMax);
                theta = thetalis[index];
                phi = philis[index];
                profil[s] += Gamma(theta, phi, alphalis[index], delta, A, factort, factorx)*cexp(-2.0*phi*I);
                profilconj[s] += conj(profil[s]);
            }
        }
    }
}

cdouble Gamma(double theta, double phi, double alpha, double delta, double A, double factort, double factorx){
    double gammat = 1. + (factort*cos((2.*phi)+(2.*alpha)));
    double gammax = factorx*sin((2.*phi)+(2.*alpha));
    cdouble factor = -A*delta*cexp(2.*phi*I)/(pow(theta, delta)*(delta-2.));
    return factor*(gammat + (gammax*I));
}

void get_batch_size(int data_size, int nproc, int **local_size_out, int **displs_out) {
    // calc the batch size
    int *local_size = malloc1d(int, nproc);
    int *displs = malloc1d(int, nproc);
    memCheck(local_size); memCheck(displs)
    int lsize = data_size/nproc;
    int mod = (data_size % nproc);
    int lsize1 = lsize + 1;
    for (int i=0; i < nproc; i++) {
        local_size[i] = i < mod ? lsize1 : lsize;
    }
    int sum=local_size[0];
    displs[0] = 0;
    for (int i=1; i < nproc; i++) {
        displs[i] = sum;
        sum += local_size[i];
    }
    *local_size_out = local_size;
    *displs_out = displs;
}
