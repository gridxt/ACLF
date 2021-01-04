#ifndef YBUSSPARSE_H_INCLUDED
#define YBUSSPARSE_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "mkl.h"
#include "mkl_service.h"
#include "Structs.h"
#include "NumericalFunctions.h"
#include "mkl_spblas.h"
#include <omp.h>

//int dense2CSR(MKL_Complex16 *A, MKL_INT M, MKL_INT N, MKL_INT nzmax, MKL_Complex16 *A_Vals, MKL_INT *A_Row, MKL_INT *A_Col);
//int CSR2dense(MKL_Complex16 *A, MKL_INT M, MKL_INT N, MKL_Complex16 *A_Vals, MKL_INT *A_Row, MKL_INT *A_Col);
int MakeYbusSparse(struct Ybus_Struct *YbusOut, double BaseMVA, unsigned long long int Num_of_Buses, unsigned long long int Num_of_Branches, int *FromBus, int *ToBus, int *BranchStatus,
	double *r, double *x, double *b, double *TapRatio, double *ShiftAngle, double *Gs, double *Bs, double **timer);


#endif // YBUS_H_INCLUDED

