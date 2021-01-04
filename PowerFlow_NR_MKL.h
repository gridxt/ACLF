#ifndef POWERFLOW_NR_MKL_H_INCLUDED
#define POWERFLOW_NR_MKL_H_INCLUDED

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
#include <omp.h>

// DENSE FUNCTIONS
void Sbus(unsigned long long int Num_Buses, unsigned long long int Gen_Num_Buses, double BaseMVA, double *Pd, double *Qd, double *Pg, double *Qg, int *GenBusNo, int *GenStatus, double *Sbus_Real, double *Sbus_Imag, int OMP_Enbale, int OMP_Cores);
void Power_Mimatch(MKL_Complex16 *Ybus_Row_, sparse_matrix_t Ybus_CSR, MKL_Complex16 *V_Rec, MKL_Complex16 *YbusV, double *Sbus_re, double *Sbus_img, unsigned long long int Num_Buses, double *Mis_re, double *Mis_img, int OMP_Enbale, int OMP_Cores, int Sparse_Enable);
void F_Build(double *Mis_re, double *Mis_img, unsigned long long int Num_Buses, unsigned long long int PV_Elem, unsigned long long int PQ_Elem, int *PVBuses, int *PQBuses, double *F_, double *F_Norm, int OMP_Enbale, int OMP_Cores);
void dS_dV(MKL_Complex16 *Ybus_Row_, MKL_Complex16 *V_Rec, MKL_Complex16 *Ibus, MKL_Complex16 *diagVnorm, MKL_Complex16 *YbusdiagVnorm, MKL_Complex16 *YbusdiagV, double *V_mag, unsigned long long int Num_Buses, MKL_Complex16 *ds_dVa, MKL_Complex16 *ds_dVm, int OMP_Enbale, int OMP_Cores);
void Jacobian(MKL_Complex16 *ds_dVm, MKL_Complex16 *ds_dVa, int *PVPQ, unsigned long long int Num_Buses, unsigned long long int PV_Elem, unsigned long long int PQ_Elem, int *PVBuses, int *PQBuses, double *J, int OMP_Enbale, int OMP_Cores);
void Voltage_Update(double *V_Mag, double *V_Angl, MKL_Complex16 *V_Rectangular, unsigned long long int Num_of_Buses, unsigned long long int PV_Buses_Num_Elements, int *PVBuses, double *F, unsigned long long int PQ_Buses_Num_Elements, int *PQBuses, int OMP_Enbale, int OMP_Cores);

// SPARSE FUNCTIONS
MKL_Complex16 find(MKL_Complex16 *VAL, MKL_INT *JA, MKL_INT *IA, unsigned long long int row, unsigned long long int col);
int pardiso_solver(double *a, MKL_INT *ia, MKL_INT *ja, double *b, double *x, MKL_INT n, int OMP_Enbale, int OMP_Cores);
int dS_dV_Sparse(MKL_Complex16 *ds_dVa_vals, MKL_Complex16 *ds_dVm_vals, MKL_INT *ds_dV_IA, MKL_INT *ds_dV_JA,
	unsigned long long int Ybus_nnz, sparse_matrix_t Ybus_CSR, sparse_matrix_t  diagIbus_CSR, sparse_matrix_t diagV_CSR, sparse_matrix_t diagVnorm_CSR,
	MKL_Complex16 *V_Rec, MKL_Complex16 *Ibus, double *V_mag, unsigned long long int Num_Buses, int OMP_Enbale, int OMP_Cores);
long long int BinarySearch(int *arr, long long int value, long long int lower, long long int upper);
int Jacobian_Sparse_Initialize(double **J_val, MKL_INT **J_IA, MKL_INT **J_JA, MKL_INT **J_IA_Array, MKL_INT **J_IA_Position, unsigned long long int *Jnnz_, unsigned long long int *Ref_Buses, unsigned long long int Ref_Buses_Size, unsigned long long int Ybusnnz, MKL_Complex16 *ds_dVa_vals, MKL_Complex16 *ds_dVm_vals, MKL_INT *ds_dV_IA, MKL_INT *ds_dV_JA,
	int *PVPQ, unsigned long long int Num_Buses, unsigned long long int PV_Elem, unsigned long long int PQ_Elem, int *PQBuses, int OMP_Enbale, int OMP_Cores);
void Jacobian_Sparse(double *J_val, MKL_INT *J_IA_Position, MKL_INT *J_IA_Array, unsigned long long int Jnnz, MKL_Complex16 *ds_dVa_vals, MKL_Complex16 *ds_dVm_vals, int OMP_Enbale, int OMP_Cores);
int CheckRef(unsigned long long int value, unsigned long long int *Ref_Buses, unsigned long long int Ref_Buses_Size);

// This function solves the power flow using full newton method
int PF_NR_MKL(MKL_Complex16 *Ybus_Row, MKL_Complex16 *Ybus_val, MKL_INT *Ybus_JA, MKL_INT *Ybus_IA, 
		   double *V_Mag_Init, double *V_Angl_Init, unsigned long long int Num_of_Buses, int *PVBuses, int *PQBuses,
           int *GenBusNo, int *GenStatus, double BaseMVA, double *Pd, double *Qd, double *Pg, double *Qg, double *Vg,
		   unsigned long long int OnlineGen_Num_Elements, unsigned long long int PV_Buses_Num_Elements, unsigned long long int PQ_Buses_Num_Elements, unsigned long long int *Ref_Buses, unsigned long long int Ref_Buses_Size,
		   double *V_Mag, double *V_Angl,
           int FlatStart, int Max_iter, double Tol, int OMP_Enbale, int OMP_Cores, int Sparse_Enable, int Enable_Memory_Stat, double *Peak_Mem_Usage, double *timer); // end of function

#endif // POWERFLOW_NR_H_INCLUDED
