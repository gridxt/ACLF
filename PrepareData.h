#ifndef PREPAREDATA_H_INCLUDED
#define PREPAREDATA_H_INCLUDED

#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include "mkl.h"

// This function can return the new bus number associated with the original bus number after changing bus numbers to and internal indexing
int ReturnNewBusNo(int *Org, int NumberOfBuses, int BusNo);

// This function can return the original bus number given the new bus number
int ReturnOldBusNo(int *New, int NumberOfBuses, int BusNo);


// This function returns the newly assigned bus number given the original bus number
int ReturnNew(int *data, int BusNo, int NumberOfBuses);


int PrepareData(double *RawBusData, double *RawBranchData, double *RawGenData, double *RawGenCostData,
				unsigned long long int RawBusData_Number_Of_Elements, unsigned long long int RawBranchData_Number_Of_Elements, unsigned long long int RawGenData_Number_Of_Elements, unsigned long long int RawGenCostData_Number_Of_Elements,
                int BusDataColumns, int BranchDataColumns, int GenDataColumns, int GenDataCostColumns, int NumberOfBuses, int Convert_Kilo, int Convert_Ohm, double BaseMVA,
                int **BusNo_, int **BusType_, double **Pd_, double **Qd_, double **Gs_, double **Bs_, int **Area_, double **Vm_, double **Va_, double **BaseKV_, int **Zone_, double **Vmax_, double **Vmin_,
                int **FromBus_, int **ToBus_, double **r_, double **x_, double **b_, double **RateA_, double **RateB_, double **RateC_, double **TapRatio_, double **ShiftAngle_, int **BranchStatus_, double **AngMin_, double **AngMax_,
                int **GenBusNo_, double **Pg_, double **Qg_, double **Qmax_, double **Qmin_, double **Vg_, double **mBase_, int **GenStatus_, double **Pmax_, double **Pmin_, double **PC1_, double **PC2_, double **Qc1min_, double **Qc1max_, double **Qc2min_, double **Qc2max_, double **Ramp_AGC_, double **Ramp_10_, double **Ramp_30_, double **Ramp_q_, double **APF_,
                double **GenCostParam1_, double **GenCostParam2_, double **GenCostParam3_, double **GenCostParam4_, double **GenCostParam5_, double **GenCostParam6_, double **GenCostParam7_,
				unsigned long long int **Ref_Buses, unsigned long long int *Ref_Bus_Size, int **PQBuses, int **PVBuses, int **OnlineGenBuses, unsigned long long int *PQ_Buses_Number_Of_Elements, unsigned long long int *PV_Buses_Number_Of_Elements, unsigned long long int *OnlineGen_Buses_Number_Of_Elements, unsigned long long int *Gen_Buses_Number_Of_Elements,
                double **timer);

#endif // PREPAREDATA_H_INCLUDED
