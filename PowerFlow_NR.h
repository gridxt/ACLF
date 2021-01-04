#ifndef POWERFLOW_NR_H_INCLUDED
#define POWERFLOW_NR_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "NumericalFunctions.h"
#include "mkl_lapacke.h"
#include "mkl.h"
#include <omp.h>

void Sbus_NR(int NumberOfBuses, int Gen_Num_Buses, double BaseMVA, double *Pd, double *Qd, double *Pg, double *Qg, int *GenBusNo, int *GenStatus, double **Sbus_Real, double **Sbus_Imag);

void Power_Mimatch_NR(MKL_Complex16 *Ybus_Row_, double *V_mag, double *V_ang, double *Sbus_re, double *Sbus_img, int Num_Buses, double **Mis_re, double **Mis_img);

void F_Build_NR(double *Mis_re, double *Mis_img, int Num_Buses, int PV_Elem, int PQ_Elem, int *PVBuses, int *PQBuses, int Ref_Bus_No, double **F_, double *F_Norm);

void dS_dV_NR(MKL_Complex16 *Ybus_Row_, double *V_mag, double *V_ang, int Num_Buses, double **dS_dVm_re, double **dS_dVm_img, double **dS_dVa_re, double **dS_dVa_img);

void Jacobian_NR(double **dS_dVm_re, double **dS_dVm_img, double **dS_dVa_re, double **dS_dVa_img, int Num_Buses, int PV_Elem, int PQ_Elem, int *PVBuses, int *PQBuses, int Ref_Bus_No, double *J);



// This function solves the power flow using full newton method
int PF_NR(MKL_Complex16 *Ybus_Row, double *V_Mag_, double *V_Angl_, int Num_of_Buses,
           int *PVBuses, int *PQBuses,
           int *GenBusNo, int *GenStatus, double BaseMVA, double *Pd, double *Qd, double *Pg, double *Qg, double *Vg,
           int OnlineGen_Num_Elements, int PV_Buses_Num_Elements, int PQ_Buses_Num_Elements, int Ref_Bus,
           int FlatStart, int Max_iter, double Tol, long double **CPU_Execution_TIME); // end of function

#endif // POWERFLOW_NR_H_INCLUDED
