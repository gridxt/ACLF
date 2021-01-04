#ifndef YBUS_H_INCLUDED
#define YBUS_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "mkl.h"
#include "NumericalFunctions.h"

int MakeYbus(double BaseMVA, int Num_of_Buses, int Num_of_Branches, int *FromBus, int *ToBus, int *BranchStatus, double *r, double *x, double *b, double *TapRatio, double *ShiftAngle, double *Gs, double *Bs, MKL_Complex16 *Ybus_Row, double **CPU_Execution_TIME);


#endif // YBUS_H_INCLUDED

