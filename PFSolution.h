#include <stdio.h>
#include <math.h>
#include <time.h>
#include "NumericalFunctions.h"
#include "mkl_lapacke.h"
#include "mkl.h"
#include "mkl_service.h"
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"
#include "Structs.h"
#include <omp.h>


void PF_Solution(struct PFSol_Struct *PFSolStruct, struct Ybus_Struct *YbusStruct,
	unsigned long long int Num_of_Branches, int *FromBus, int *ToBus, int *BranchStatus, double *b,
	unsigned long long int Num_of_Buses, int *BusNo, int *BusType, double *Pd, double *Qd, double *TapRatio, double *Gs, double *Bs,
	unsigned long long int OnlineGenData_Num_Of_Elements, int *OnlineGen_Bus_Number, unsigned long long int Gen_Buses_Number_Of_Elements, int *GenBusNo, int *GenStatus, double *Pg, double *Qg, double *Qmin, double *Qmax, double *Pmax,
	double BaseMVA, int Ref_Bus, unsigned long long int *Ref_Buses, unsigned long long int Ref_Buses_Size,
	double *Vm_cal, double *Va_cal,
	int OMP_Enbale, int OMP_Cores, int Sparse_Enable, int Enable_Memory_Stat);