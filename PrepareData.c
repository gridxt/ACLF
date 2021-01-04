#include "PrepareData.h"

// This function can return the new bus number associated with the original bus number after changing bus numbers to and internal indexing
int ReturnNewBusNo(int *Org, int NumberOfBuses, int BusNo){
int i = 0;
while (i < NumberOfBuses && BusNo != Org[i]) {
      i++;
}

 if (i < NumberOfBuses) {
      return i; // BusNo_[i] returns the new bus number
   } else {
      return -1;
   }

}

// This function can return the original bus number given the new bus number
int ReturnOldBusNo(int *New, int NumberOfBuses, int BusNo){
int i = 0;
while (i < NumberOfBuses && BusNo != New[i]) {
      i++;
}

 if (i < NumberOfBuses) {
      return i; // Original_Bus_Number[i] return the old bus number
   } else {
      return -1;
   }

}



// Return the newly assigned bus number given the original bus number
int ReturnNew(int *data, int BusNo, int NumberOfBuses) {
	int j = 0;
	
	for (j = NumberOfBuses; j < (2*NumberOfBuses); j++) {
		if (data[j] == BusNo) return data[j - NumberOfBuses]; // return the new bus number if found	
	}

	return -1; // if the bus number does not exist, return -1

}






int PrepareData(double *RawBusData, double *RawBranchData, double *RawGenData, double *RawGenCostData,
				unsigned long long int RawBusData_Number_Of_Elements, unsigned long long int RawBranchData_Number_Of_Elements, unsigned long long int RawGenData_Number_Of_Elements, unsigned long long int RawGenCostData_Number_Of_Elements,
                int BusDataColumns, int BranchDataColumns, int GenDataColumns, int GenDataCostColumns, int NumberOfBuses, int Convert_Kilo, int Convert_Ohm, double BaseMVA,
                int **BusNo_, int **BusType_, double **Pd_, double **Qd_, double **Gs_, double **Bs_, int **Area_, double **Vm_, double **Va_, double **BaseKV_, int **Zone_, double **Vmax_, double **Vmin_,
                int **FromBus_, int **ToBus_, double **r_, double **x_, double **b_, double **RateA_, double **RateB_, double **RateC_, double **TapRatio_, double **ShiftAngle_, int **BranchStatus_, double **AngMin_, double **AngMax_,
                int **GenBusNo_, double **Pg_, double **Qg_, double **Qmax_, double **Qmin_, double **Vg_, double **mBase_, int **GenStatus_, double **Pmax_, double **Pmin_, double **PC1_, double **PC2_, double **Qc1min_, double **Qc1max_, double **Qc2min_, double **Qc2max_, double **Ramp_AGC_, double **Ramp_10_, double **Ramp_30_, double **Ramp_q_, double **APF_,
                double **GenCostParam1_, double **GenCostParam2_, double **GenCostParam3_, double **GenCostParam4_, double **GenCostParam5_, double **GenCostParam6_, double **GenCostParam7_,
				unsigned long long int **Ref_Buses, unsigned long long int *Ref_Bus_Size, int **PQBuses, int **PVBuses, int **OnlineGenBuses, unsigned long long int *PQ_Buses_Number_Of_Elements, unsigned long long int *PV_Buses_Number_Of_Elements, unsigned long long int *OnlineGen_Buses_Number_Of_Elements, unsigned long long int *Gen_Buses_Number_Of_Elements,
                double **timer){

double start, stop;
start = dsecnd();
unsigned long long int i = 0, j = 0, t = 0;
unsigned long long int Max_Bus_Num = 0; // What is the largest bus number used in the bus data
unsigned long long int PQ_Buses_Count = 0;
unsigned long long int PV_Buses_Count = 0;
unsigned long long int OnlineGen_Buses_Count = 0;
unsigned long long int Gen_buses_count = 0;
int IsRefDefined = 0; // Check if the reference bus is defined in the data 0: No -- 1: Yes
int *GenStatusTemp = (int*)calloc(NumberOfBuses, sizeof(int*)); // Temporary keep the number of ON generators at each bus


// Assign RawBusData content to specified variables
// Bus Data
for (i = 1; i <= RawBusData_Number_Of_Elements; i = i + BusDataColumns){
    (*BusNo_)[j] = RawBusData[i-1];
    (*BusType_)[j] = RawBusData[i-1+1];
    if (Convert_Kilo == 1){ // Convert from kW and kVAR to MW and MVAR
        (*Pd_)[j] = RawBusData[i-1+2]/1000;
        (*Qd_)[j] = RawBusData[i-1+3]/1000;
    }
    else{
        (*Pd_)[j] = RawBusData[i-1+2];
        (*Qd_)[j] = RawBusData[i-1+3];
    }
    (*Gs_)[j] = RawBusData[i-1+4];
    (*Bs_)[j] = RawBusData[i-1+5];
    (*Area_)[j] = RawBusData[i-1+6];
    (*Vm_)[j] = RawBusData[i-1+7];
    (*Va_)[j] = RawBusData[i-1+8];
    (*BaseKV_)[j] = RawBusData[i-1+9];
    (*Zone_)[j] = RawBusData[i-1+10];
    (*Vmax_)[j] = RawBusData[i-1+11];
    (*Vmin_)[j] = RawBusData[i-1+12];

    if ((*BusNo_)[j] > Max_Bus_Num) Max_Bus_Num = (*BusNo_)[j]; // Largest bus number in data

    j++;
}


// Check if the bus numbers are sequential. If not, change the order.
if (Max_Bus_Num > NumberOfBuses) {
	printf("\n	Reordering bus numbers\n");	fflush(stdout);
	//*BusNo_ = realloc(*BusNo_, (NumberOfBuses * 2) * sizeof(int));	// Double the size of BusNo_ vector
	for (j = 0; j < NumberOfBuses; j++) {	// create this vector [new_bus_numbers     original_bus_number]		
		(*BusNo_)[j + NumberOfBuses] = (*BusNo_)[j];
		(*BusNo_)[j] = j + 1;
	}

}
else {
	for (j = 0; j < NumberOfBuses; j++) {	// create this vector [original_bus_numbers     original_bus_number]		
		(*BusNo_)[j + NumberOfBuses] = (*BusNo_)[j];		
	}
}



j = 0;

// Branch Data
for (i = 1; i <= RawBranchData_Number_Of_Elements; i = i + BranchDataColumns) {
	if (Max_Bus_Num <= NumberOfBuses) { // Check if internal indexing of bus numbers was performed
		(*FromBus_)[j] = RawBranchData[i - 1];
		(*ToBus_)[j] = RawBranchData[i - 1 + 1];
	}
	else {
		(*FromBus_)[j] = ReturnNew(*BusNo_, RawBranchData[i - 1], NumberOfBuses);		
		(*ToBus_)[j] = ReturnNew(*BusNo_, RawBranchData[i - 1 + 1], NumberOfBuses); 
	}

	(*b_)[j] = RawBranchData[i - 1 + 4];
	if (Convert_Ohm == 1) { // Convert from Ohm to P.U.
		(*r_)[j] = RawBranchData[i - 1 + 2] / (((*BaseKV_)[j])*((*BaseKV_)[j]) / (BaseMVA));
		(*x_)[j] = RawBranchData[i - 1 + 3] / (((*BaseKV_)[j])*((*BaseKV_)[j]) / (BaseMVA));
	}
	else {
		(*r_)[j] = RawBranchData[i - 1 + 2];
		(*x_)[j] = RawBranchData[i - 1 + 3];
	}
	(*RateA_)[j] = RawBranchData[i - 1 + 5];
	(*RateB_)[j] = RawBranchData[i - 1 + 6];
	(*RateC_)[j] = RawBranchData[i - 1 + 7];
	(*TapRatio_)[j] = RawBranchData[i - 1 + 8];
	//if ((*TapRatio_)[j] == 0) { (*TapRatio_)[j] = 1; } // if tap ratio is 0, change it to 1
	(*ShiftAngle_)[j] = RawBranchData[i - 1 + 9];
	(*BranchStatus_)[j] = RawBranchData[i - 1 + 10];
	(*AngMin_)[j] = RawBranchData[i - 1 + 11];
	(*AngMax_)[j] = RawBranchData[i - 1 + 12];
	j++;
}

j = 0;
// Generators Data
for (i = 1; i <= RawGenData_Number_Of_Elements; i = i + GenDataColumns) {	
		if (Max_Bus_Num <= NumberOfBuses) { // Check if internal indexing of bus numbers was performed
			(*GenBusNo_)[j] = RawGenData[i - 1];
		}
		else {
			(*GenBusNo_)[j] = ReturnNew(*BusNo_, RawGenData[i - 1], NumberOfBuses);
		}

		(*Pg_)[j] = RawGenData[i - 1 + 1];
		(*Qg_)[j] = RawGenData[i - 1 + 2];
		(*Qmax_)[j] = RawGenData[i - 1 + 3];
		(*Qmin_)[j] = RawGenData[i - 1 + 4];
		(*Vg_)[j] = RawGenData[i - 1 + 5];
		(*mBase_)[j] = RawGenData[i - 1 + 6];
		(*GenStatus_)[j] = RawGenData[i - 1 + 7];
		if ((*GenStatus_)[j] != 0 && (*GenStatus_)[j] != 1) {
			printf("\nInvalid Status detected in Generators data!\n"); fflush(stdout);
			exit(0);
		}
		(*Pmax_)[j] = RawGenData[i - 1 + 8];
		(*Pmin_)[j] = RawGenData[i - 1 + 9];
		(*PC1_)[j] = RawGenData[i - 1 + 10];
		(*PC2_)[j] = RawGenData[i - 1 + 11];
		(*Qc1min_)[j] = RawGenData[i - 1 + 12];
		(*Qc1max_)[j] = RawGenData[i - 1 + 13];
		(*Qc2min_)[j] = RawGenData[i - 1 + 14];
		(*Qc2max_)[j] = RawGenData[i - 1 + 15];
		(*Ramp_AGC_)[j] = RawGenData[i - 1 + 16];
		(*Ramp_10_)[j] = RawGenData[i - 1 + 17];
		(*Ramp_30_)[j] = RawGenData[i - 1 + 18];
		(*Ramp_q_)[j] = RawGenData[i - 1 + 19];
		(*APF_)[j] = RawGenData[i - 1 + 20];

		if ((*GenStatus_)[j] == 1) { // Save online generators' bus number
			(*OnlineGenBuses)[OnlineGen_Buses_Count] = (*GenBusNo_)[j];
			OnlineGen_Buses_Count = OnlineGen_Buses_Count + 1;
			*OnlineGenBuses = realloc(*OnlineGenBuses, (OnlineGen_Buses_Count + 1) * sizeof(int));
		}
		Gen_buses_count = Gen_buses_count + 1;
		j++;	
}
*OnlineGen_Buses_Number_Of_Elements = OnlineGen_Buses_Count;
*Gen_Buses_Number_Of_Elements = Gen_buses_count;

j = 0;
for (i = 1; i <= RawGenCostData_Number_Of_Elements; i = i + GenDataCostColumns) {
	(*GenCostParam1_)[j] = RawGenCostData[i - 1];
	(*GenCostParam2_)[j] = RawGenCostData[i - 1 + 1];
	(*GenCostParam3_)[j] = RawGenCostData[i - 1 + 2];
	(*GenCostParam4_)[j] = RawGenCostData[i - 1 + 3];
	(*GenCostParam5_)[j] = RawGenCostData[i - 1 + 4];
	(*GenCostParam6_)[j] = RawGenCostData[i - 1 + 5];
	(*GenCostParam7_)[j] = RawGenCostData[i - 1 + 6];
	j++;
}


// Count the number of online generators at each bus
for (j = 0; j < OnlineGen_Buses_Count; j++) {
	i = (*OnlineGenBuses)[j];
	GenStatusTemp[i - 1] = GenStatusTemp[i - 1] + 1;	
}

// First check if the reference bus exist. Also check if there is only one reference bus.
for (j = 0; j < NumberOfBuses; j++) {
	if ((*BusType_)[j] == 3 && GenStatusTemp[(*BusNo_)[j] - 1] > 0) { 
		IsRefDefined = IsRefDefined + 1; 
	}

}

if (IsRefDefined > 1) {
	// more than one valid reference bus
	*Ref_Buses = (unsigned long long int *)calloc((size_t)(IsRefDefined), sizeof(unsigned long long int));
	*Ref_Bus_Size = IsRefDefined;
}
else{
	// no reference bus or exactly one
	*Ref_Buses = (unsigned long long int *)calloc((size_t)(1), sizeof(unsigned long long int));
	*Ref_Bus_Size = 1;
}

if (IsRefDefined > 1) printf("  	Multiple reference Buses detected!\n"); fflush(stdout);

// Save bus types
i = 0;
for (j = 0; j < NumberOfBuses; j++) {

	// Reference Bus (there is exactly one valid ref. bus with an online generator)
	if ((*BusType_)[j] == 3 && IsRefDefined == 1 && GenStatusTemp[(*BusNo_)[j] - 1] > 0) {
		(*Ref_Buses)[0] = (*BusNo_)[j];
	}

	// There is more than one valid reference bus with online generator
	if ((*BusType_)[j] == 3 && IsRefDefined > 1 && GenStatusTemp[(*BusNo_)[j] - 1] > 0) {		
		(*Ref_Buses)[i] = (*BusNo_)[j];		
		i++;
	}


	// PQ: Consider PV or REF buses with OFFLINE generator as PQ bus.
	if ((*BusType_)[j] == 1 || GenStatusTemp[(*BusNo_)[j] - 1] == 0) { // Save the PQ bus number. 
		(*PQBuses)[PQ_Buses_Count] = (*BusNo_)[j];
		PQ_Buses_Count = PQ_Buses_Count + 1;
		*PQBuses = realloc(*PQBuses, (PQ_Buses_Count + 1) * sizeof(int));
	}


	// PV
	if ((*BusType_)[j] == 2 && GenStatusTemp[(*BusNo_)[j] - 1] > 0) { // Save the PV bus number if the generator on this bus is ON
		//Reference bus has no online generator or is not defined
		if (IsRefDefined == 0 && PV_Buses_Count == 0) { // if the generator on ref. bus is OFF, select the FIRST PV bus as ref.
			(*Ref_Buses)[0] = (*BusNo_)[j];
			IsRefDefined = -1; // stop searching for ref.bus
		}
		// PV Bus
		else {
			(*PVBuses)[PV_Buses_Count] = (*BusNo_)[j];
			PV_Buses_Count = PV_Buses_Count + 1;
			*PVBuses = realloc(*PVBuses, (PV_Buses_Count + 1) * sizeof(int));
		}
	}

	
	// If bus type is invalid
	if ((*BusType_)[j] != 1 && (*BusType_)[j] != 2 && (*BusType_)[j] != 3) { 
		printf("\nInvalid bus type detected in Bus data!!\n"); fflush(stdout);
		exit(0);
	}

}

*PQ_Buses_Number_Of_Elements = PQ_Buses_Count;
*PV_Buses_Number_Of_Elements = PV_Buses_Count;


  stop = dsecnd();
  (*timer)[0] = stop - start;

  free(GenStatusTemp);

return 1;
}
