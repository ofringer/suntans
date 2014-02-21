#ifndef _sunjoin_h
#define _sunjoin_h

typedef struct _ptrT {
    int **mnptr;
    int **eptr;

    //Temporary arrays for netcdf reading
    double *ncscratch;
    int *nctmpint;

} ptrT;
ptrT *ptr;

#endif
