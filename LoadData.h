#ifndef LOADDATA_H_INCLUDED
#define LOADDATA_H_INCLUDED
#define bufSize 10240
#pragma warning (disable : 1786)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "mkl.h"

int LoadCaseData(char *filename, double **Bus_Data, unsigned long long int *Num_of_Bus, double **Branch_Data, unsigned long long int *Num_of_Branch, double **Gen_Data, unsigned long long int *Num_of_Gen,
	double **Gen_Cost_Data, unsigned long long int *Num_of_Gen_Cost, char Bus_Data_Label[], char Branch_Data_Label[], char Gen_Data_Label[], char GenCost_Data_Label[], char BaseMVA_Label[], char Convert_Kilo_Mega_Label[], char Convert_Ohm_PU_Label[], char End_Data_Label[],
	int *BusDataColumns, int *BranchDataColumns, int *GenDataColumns, int *GenDataCostColumns, double *BaseMVA, int *Convert_Kilo, int *Convert_Ohm, double *timer);

#endif // LOADDATA_H_INCLUDED



