#include "PFSolution.h"
#define ALIGN 64

/* To avoid constantly repeating the part of code that checks different functions' status, using the below macros */
#define CHECK_SPARSE(function)  do {  if(function != SPARSE_STATUS_SUCCESS){ status = 2; goto memory_free; } } while(0) // Check sparse function status
#define CHECK_MEM(variable)  do {  if(variable == NULL){ status = 1; goto memory_free; } } while(0) // check memory allocation


// ****************** Main Function ******************
void PF_Solution(struct PFSol_Struct *PFSolStruct, struct Ybus_Struct *YbusStruct,
	unsigned long long int Num_of_Branches, int *FromBus, int *ToBus, int *BranchStatus, double *b,
	unsigned long long int Num_of_Buses, int *BusNo, int *BusType, double *Pd, double *Qd, double *TapRatio, double *Gs, double *Bs,
	unsigned long long int OnlineGenData_Num_Of_Elements, int *OnlineGen_Bus_Number, unsigned long long int Gen_Buses_Number_Of_Elements, int *GenBusNo, int *GenStatus, double *Pg_, double *Qg_, double *Qmin_, double *Qmax_, double *Pmax,
	double BaseMVA, int Ref_Bus, unsigned long long int *Ref_Buses, unsigned long long int Ref_Buses_Size,
	double *Vm_cal, double *Va_cal,
	int OMP_Enbale, int OMP_Cores, int Sparse_Enable, int Enable_Memory_Stat)
{

	/* **** Initial Variables (used in both sparse and dense functions) **** */
	double timer_stop[10], timer_start[10]; // Holds the execution time of different parts

	timer_start[0] = dsecnd(); // Total execution time
	if (Enable_Memory_Stat == 1) mkl_peak_mem_usage(MKL_PEAK_MEM_ENABLE);

	unsigned long long int i = 0, j = 0, k = 0, m = 0, dummy = 0;	
	unsigned long long int Chunk_Size;
	unsigned long long int Num_Online_Gens_No_PQ = 0;
	long long int Bus_Number = 0;
	unsigned long long int Num_of_Nonzero = 0;
	unsigned long long int Num_Online_Branches = 0;
	unsigned long long int Num_Offline_Branches = 0;
	unsigned long long int nnz = 0;
	unsigned long long int refgen;
	double eps = pow(2,-52);
	double mis = 0, sum = 0, temp = 0;
	int status = 0;

	MKL_Complex16 alpha, beta;
	alpha.real = 1;
	alpha.imag = 0;
	beta.real = 0;
	beta.imag = 0;

	int *GenBus = NULL; // Holds bus number for online generators which are not located at PQ buses				
	int *GenCount = NULL;	
	double *Pd_gbus = NULL;
	double *Qd_gbus = NULL;
	double *Qg_total = NULL; //vector of total Qg at each bus
	double *Qg_min = NULL; // vector of min total Qg at each bus
	double *Qg_max = NULL; //vector of max total Qg at each bus
	double *Qmin = NULL;
	double *Qmax = NULL;
	double P_loss_max = 0; // used in finding the branch with maximum P loss
	double Q_loss_max = 0;
	double Vm_max = 0;
	double Vm_min = 0;
	double Vang_max = 0;
	double Vang_min = 0;

	
	MKL_Complex16 *V_rect = NULL;
	MKL_Complex16 *V_GenBus = NULL;
	MKL_Complex16 *Ybus = NULL;
	MKL_Complex16 *YbusV = NULL;
	MKL_Complex16 *Sbus = NULL;

	MKL_Complex16 *Ybus_val = NULL;
	MKL_INT *Ybus_JA = NULL;
	MKL_INT *Ybus_IA = NULL;
	
	MKL_Complex16 *Sf = NULL; // complex power at "from" bus 
	MKL_Complex16 *St = NULL; // complex power injected at "to" bus
	sparse_matrix_t Ybus_CSR = NULL, Yf_CSR = NULL, Yt_CSR = NULL;

	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	
	MKL_Complex16 *Yf_New = NULL;
	MKL_INT *Yf_New_IA = NULL;
	MKL_INT *Yf_New_JA = NULL;

	MKL_Complex16 *Yt_New = NULL;
	MKL_INT *Yt_New_IA = NULL;
	MKL_INT *Yt_New_JA = NULL;
	
	MKL_Complex16 *Vdrop = NULL;	
	MKL_Complex16 *loss = NULL;

	MKL_Complex16 *Vf = NULL, *Vt = NULL;
	MKL_Complex16 *V_drop_diag = NULL, *Vrect_conj_diag = NULL, *VddA = NULL, *dYsc = NULL, *V_mag = NULL, *Vm_diag_recip = NULL, *Bc_diag = NULL, *Tap_diag = NULL, *B = NULL, *B_conj_B = NULL, *dlosstemp = NULL;
	MKL_Complex16 *BcTap = NULL, *CfVm = NULL, *BcTapCfVm = NULL, *CfVm_diag = NULL;
	MKL_Complex16 *CtVm = NULL, *BcCtVm = NULL, *CtVm_diag = NULL;


	// Sparse values
	MKL_Complex16 *dYsc_val = NULL, *dYsc_mod_val = NULL, *Bc_diag_val = NULL, *Tap_diag_val = NULL;
	MKL_Complex16 *V_mag_val = NULL, *Vm_recip_val = NULL, *Vrect_conj_val = NULL;
	MKL_INT *IA_Bus = NULL, *JA_Bus = NULL;
	MKL_INT *IA_Br = NULL, *JA_Br = NULL;

	// CSR handles
	sparse_matrix_t Tap_CSR_Handle = NULL, dYsc_CSR_Handle = NULL, dYsc_mod_CSR_Handle = NULL, Bc_CSR_Handle = NULL, Vdrop_CSR_Handle = NULL;
	sparse_matrix_t  Vm_reci_CSR_Handle = NULL, Vconj_CSR_Handle = NULL;
	sparse_matrix_t VddA_CSR_Handle = NULL, B_CSR_Handle = NULL, B_m_conjB_CSR_Handle = NULL, B_p_conjB_CSR_Handle = NULL, dYscBp_CSR_Handle = NULL;
	sparse_matrix_t BcTap_CSR_Handle = NULL, CfVm_CSR_Handle = NULL, CfVmCf_CSR_Handle = NULL, CtVm_CSR_Handle = NULL, CtVmCt_CSR_Handle = NULL;

	// CSR export variables
	MKL_INT rows, cols, rows1, cols1, rows2, cols2;
	sparse_index_base_t indexing = 0, indexing1 = 0, indexing2 = 0;
	MKL_INT *Dummy_JA_1 = NULL, *Dummy_pointerE_1 = NULL, *Dummy_IA_1 = NULL;
	MKL_INT *Dummy_JA_2 = NULL, *Dummy_pointerE_2 = NULL, *Dummy_IA_2 = NULL;
	MKL_Complex16 *B_m_conjB_val = NULL, *B_p_conjB_val = NULL;


	// First find the number of online generators which are not located at PQ buses		
	for (i = 0; i < OnlineGenData_Num_Of_Elements; i++) {
		Bus_Number = OnlineGen_Bus_Number[i] - 1;
		if (BusType[Bus_Number] != 1) 	Num_Online_Gens_No_PQ = Num_Online_Gens_No_PQ + 1;	// Exclude generators at PQ buses	
	}
	
	// Allocate memeory		
	GenBus = (int*)calloc((size_t)Num_Online_Gens_No_PQ, sizeof(int*)); // Holds bus number for online generators which are not located at PQ buses	
	CHECK_MEM(GenBus);
	Pd_gbus = (double*)calloc((size_t)Num_Online_Gens_No_PQ, sizeof(double*));
	CHECK_MEM(Pd_gbus);
	Qd_gbus = (double*)calloc((size_t)Num_Online_Gens_No_PQ, sizeof(double*));
	CHECK_MEM(Qd_gbus);
	GenCount = (int*)calloc((size_t)Num_of_Buses, sizeof(int*));
	CHECK_MEM(GenCount);
	Qg_total = (double*)calloc((size_t)Num_of_Buses, sizeof(double*));
	CHECK_MEM(Qg_total);
	Qg_min = (double*)calloc((size_t)Num_of_Buses, sizeof(double*));
	CHECK_MEM(Qg_min);
	Qg_max = (double*)calloc((size_t)Num_of_Buses, sizeof(double*));
	CHECK_MEM(Qg_max);

	YbusV = (MKL_Complex16 *)mkl_calloc((size_t)Num_Online_Gens_No_PQ, sizeof(MKL_Complex16), ALIGN);
	CHECK_MEM(YbusV);
	V_GenBus = (MKL_Complex16 *)mkl_calloc((size_t)Num_Online_Gens_No_PQ, sizeof(MKL_Complex16), ALIGN);
	CHECK_MEM(V_GenBus);
	V_rect = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Buses, sizeof(MKL_Complex16), ALIGN);
	CHECK_MEM(V_rect);
	Sbus = (MKL_Complex16 *)mkl_calloc((size_t)Num_Online_Gens_No_PQ, sizeof(MKL_Complex16), ALIGN);
	CHECK_MEM(Sbus);
	Vdrop = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches, sizeof(MKL_Complex16), ALIGN);
	CHECK_MEM(Vdrop);
	Vf = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches, sizeof(MKL_Complex16), ALIGN);
	CHECK_MEM(Vf);
	Vt = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches, sizeof(MKL_Complex16), ALIGN);
	CHECK_MEM(Vt);
	
	// Copies of Pg, Qg, Qmin, Qmax
	Qmin = (double*)calloc((size_t)Gen_Buses_Number_Of_Elements, sizeof(double*));
	CHECK_MEM(Qmin);
	Qmax = (double*)calloc((size_t)Gen_Buses_Number_Of_Elements, sizeof(double*));	
	CHECK_MEM(Qmax);
	(*PFSolStruct).Pg = (double*)calloc((size_t)Gen_Buses_Number_Of_Elements, sizeof(double*));
	CHECK_MEM((*PFSolStruct).Pg);
	(*PFSolStruct).Qg = (double*)calloc((size_t)Gen_Buses_Number_Of_Elements, sizeof(double*));
	CHECK_MEM((*PFSolStruct).Qg);

	// Allocate memory for the struct variables
	(*PFSolStruct).Timer = (double*)calloc(7, sizeof(double*));
	CHECK_MEM((*PFSolStruct).Timer);
	(*PFSolStruct).Vm = (double*)calloc((size_t)Num_of_Buses, sizeof(double*));
	CHECK_MEM((*PFSolStruct).Vm);
	(*PFSolStruct).Vang_deg = (double*)calloc((size_t)Num_of_Buses, sizeof(double*));
	CHECK_MEM((*PFSolStruct).Vang_deg);
	(*PFSolStruct).PF = (double*)calloc((size_t)Num_of_Branches, sizeof(double*));
	CHECK_MEM((*PFSolStruct).PF);
	(*PFSolStruct).QF = (double*)calloc((size_t)Num_of_Branches, sizeof(double*));
	CHECK_MEM((*PFSolStruct).QF);
	(*PFSolStruct).PT = (double*)calloc((size_t)Num_of_Branches, sizeof(double*));
	CHECK_MEM((*PFSolStruct).PT);
	(*PFSolStruct).QT = (double*)calloc((size_t)Num_of_Branches, sizeof(double*));	
	CHECK_MEM((*PFSolStruct).QT);
	(*PFSolStruct).Loss = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches, sizeof(MKL_Complex16), ALIGN);
	CHECK_MEM((*PFSolStruct).Loss);
	(*PFSolStruct).FCHG = (double*)calloc((size_t)Num_of_Branches, sizeof(double*));
	CHECK_MEM((*PFSolStruct).FCHG);
	(*PFSolStruct).TCHG = (double*)calloc((size_t)Num_of_Branches, sizeof(double*));
	CHECK_MEM((*PFSolStruct).TCHG);	


		
	// Make a copy of original Pg, Qg, Qmin, Qmax because we are going to change them	
	Chunk_Size = round(Gen_Buses_Number_Of_Elements / OMP_Cores);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Gen_Buses_Number_Of_Elements,PFSolStruct,Pg_,Qg_,Qmin,Qmin_,Qmax,Qmax_) if(OMP_Enbale == 1)
	for (i = 0; i < Gen_Buses_Number_Of_Elements; i++) {
		(*PFSolStruct).Pg[i] = Pg_[i];
		(*PFSolStruct).Qg[i] = Qg_[i];
		Qmin[i] = Qmin_[i];
		Qmax[i] = Qmax_[i];		
	}
	
	
	// Save the bus number of online generators that are not located at PQ buses. Create the V vector related to these buses.
	j = 0;
	for (i = 0; i < OnlineGenData_Num_Of_Elements; i++) {
		Bus_Number = OnlineGen_Bus_Number[i] - 1; // -1 because we are using 0 based indexing
		if (BusType[Bus_Number] != 1) {
			GenBus[j] = OnlineGen_Bus_Number[i];
			V_GenBus[j].real = Vm_cal[Bus_Number] * cos(Va_cal[Bus_Number]);
			V_GenBus[j].imag = Vm_cal[Bus_Number] * sin(Va_cal[Bus_Number]);			
			j++;		
		}
	}

	
	// Convert and copy calculated V values to V_rect. Also, convert the voltage anglese to degree and save them in the struct
	Chunk_Size = round(Num_of_Buses / OMP_Cores);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Num_of_Buses,V_rect,Vm_cal,Va_cal,PFSolStruct) if(OMP_Enbale == 1)
	for (i = 0; i < Num_of_Buses; i++) {
		V_rect[i].real = Vm_cal[i] * cos(Va_cal[i]);
		V_rect[i].imag = Vm_cal[i] * sin(Va_cal[i]);		
		(*PFSolStruct).Vm[i] = Vm_cal[i];
		(*PFSolStruct).Vang_deg[i] = rad2deg(Va_cal[i]);
	}
		
	
	/* ****************** Generation at Each Bus ******************** */
	timer_start[1] = dsecnd(); // Generation calculation time
	// Create a new Ybus matrix
	// *** Dense computation ***	
	if (Sparse_Enable == 0) {			
		Ybus = (MKL_Complex16 *)mkl_calloc((size_t)(Num_Online_Gens_No_PQ*Num_of_Buses), sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(Ybus);
		Chunk_Size = round(Num_of_Buses / OMP_Cores);
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i,j,Bus_Number) shared(Chunk_Size,Num_Online_Gens_No_PQ,Num_of_Buses,GenBus,Ybus,YbusStruct) collapse(2) if(OMP_Enbale == 1)
		for (i = 0; i < Num_Online_Gens_No_PQ; i++) {
			for (j = 0; j < Num_of_Buses; j++) {
				Bus_Number = GenBus[i] - 1;
				Ybus[i*Num_of_Buses + j].real = (*YbusStruct).Ybus_Dense[Bus_Number*Num_of_Buses + j].real;
				Ybus[i*Num_of_Buses + j].imag = (*YbusStruct).Ybus_Dense[Bus_Number*Num_of_Buses + j].imag;
			}
		}

	}
	// *** Sparse computation ***
	else {				
		// First count the number of nonzero elements associated with the generator buses in the original sparse Ybus
		Chunk_Size = round(Num_Online_Gens_No_PQ / OMP_Cores);
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i,Bus_Number) shared(Chunk_Size,Num_Online_Gens_No_PQ,GenBus,YbusStruct) reduction(+:Num_of_Nonzero) if(OMP_Enbale == 1)
		for (i = 0; i < Num_Online_Gens_No_PQ; i++) {			
				Bus_Number = GenBus[i] - 1;				
				Num_of_Nonzero = Num_of_Nonzero + ( (*YbusStruct).Ybus_JA[Bus_Number + 1] - (*YbusStruct).Ybus_JA[Bus_Number] );				
		}		
		// Allocate memory for the sparse Ybus matrix
		Ybus_val = (MKL_Complex16 *)mkl_calloc((size_t)(Num_of_Nonzero), sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(Ybus_val);
		Ybus_IA = (MKL_INT *)mkl_calloc((size_t)(Num_of_Nonzero), sizeof(MKL_INT), ALIGN);
		CHECK_MEM(Ybus_IA);
		Ybus_JA = (MKL_INT *)mkl_calloc((size_t)(Num_Online_Gens_No_PQ + 1), sizeof(MKL_INT), ALIGN);
		CHECK_MEM(Ybus_JA);
		// Assign the values to each array
		k = 0;
		for (i = 0; i < Num_Online_Gens_No_PQ; i++) {
			Bus_Number = GenBus[i] - 1;
			Ybus_JA[i + 1] = Ybus_JA[i] + ((*YbusStruct).Ybus_JA[Bus_Number + 1] - (*YbusStruct).Ybus_JA[Bus_Number]);
			for (j = (*YbusStruct).Ybus_JA[Bus_Number]; j < (*YbusStruct).Ybus_JA[Bus_Number + 1]; j++) {
				Ybus_val[k].real = (*YbusStruct).Ybus_Values[j].real;
				Ybus_val[k].imag = (*YbusStruct).Ybus_Values[j].imag;
				Ybus_IA[k] = (*YbusStruct).Ybus_IA[j];
				k = k + 1;
			}
		}				
	}
	
	
	// Sbus = V_GenBus.*conj(Ybus(gbus, :) * V_rect) = [Num_Online_Gens_No_PQ x 1] .* [Num_Online_Gens_No_PQ x Num_of_Buses][Num_of_Buses x 1] = [Num_Online_Gens_No_PQ x 1] 
	// YbusV = Ybus x V_rect
	if (Sparse_Enable == 0) {
		// Dense computation
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Num_Online_Gens_No_PQ, 1, Num_of_Buses, &alpha, Ybus, Num_of_Buses, V_rect, 1, &beta, YbusV, 1);
	}
	else {
		//Sparse computation
		CHECK_SPARSE(mkl_sparse_z_create_csr(&Ybus_CSR, SPARSE_INDEX_BASE_ZERO, Num_Online_Gens_No_PQ, Num_of_Buses, Ybus_JA, Ybus_JA + 1, Ybus_IA, Ybus_val));
		CHECK_SPARSE(mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, Ybus_CSR, descr, V_rect, beta, YbusV));
	}
	

	// Sbus = V_GenBus.*conj(YbusV);
	Chunk_Size = round(Num_Online_Gens_No_PQ / OMP_Cores);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Num_Online_Gens_No_PQ,Sbus,V_GenBus,YbusV) if(OMP_Enbale == 1)
	for (i = 0; i < Num_Online_Gens_No_PQ; i++) {
		Sbus[i].real = Mult_re(V_GenBus[i].real, V_GenBus[i].imag, YbusV[i].real, -1 * YbusV[i].imag);
		Sbus[i].imag = Mult_img(V_GenBus[i].real, V_GenBus[i].imag, YbusV[i].real, -1 * YbusV[i].imag);
	}
	
	
	// Total Load 
	Chunk_Size = round(Num_Online_Gens_No_PQ / OMP_Cores);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i,Bus_Number) shared(Chunk_Size,Num_Online_Gens_No_PQ,GenBus,Pd_gbus,Qd_gbus,Pd,Qd) if(OMP_Enbale == 1)
	for (i = 0; i < Num_Online_Gens_No_PQ; i++) {		
		Bus_Number = GenBus[i] - 1;
		Pd_gbus[i] = Pd[Bus_Number];
		Qd_gbus[i] = Qd[Bus_Number];				
	}

	if (Num_Online_Gens_No_PQ > 1) {
		// Count the number of online generators at each bus.  
		Chunk_Size = round(OnlineGenData_Num_Of_Elements / OMP_Cores);
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i,Bus_Number) shared(Chunk_Size,OnlineGenData_Num_Of_Elements,OnlineGen_Bus_Number) reduction(+:GenCount[:Num_of_Buses]) if(OMP_Enbale == 1)
		for (i = 0; i < OnlineGenData_Num_Of_Elements; i++) {
			Bus_Number = OnlineGen_Bus_Number[i] - 1;
			GenCount[Bus_Number]++;
		}

		// Calcualte Qg and divide by number of generators at the bus to distribute equally.
		j = 0;
		for (i = 0; i < Gen_Buses_Number_Of_Elements; i++) {
			Bus_Number = GenBusNo[i] - 1; // -1 because we are using 0 based indexing
			if (BusType[Bus_Number] != 1 && GenStatus[i] != 0) {
				(*PFSolStruct).Qg[i] = (Sbus[j].imag * BaseMVA + Qd_gbus[j]) / GenCount[Bus_Number];
				j++;
			}
			else if (GenStatus[i] == 0) {
				(*PFSolStruct).Qg[i] = 0;
			}			
		}

		
		// Replace +/- Inf limits 	
		Chunk_Size = round(Gen_Buses_Number_Of_Elements / OMP_Cores);
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i,Bus_Number) shared(Chunk_Size,Gen_Buses_Number_Of_Elements,GenBusNo,BusType,GenStatus,Qmin,Qmax,PFSolStruct) if(OMP_Enbale == 1)
		for (i = 0; i < Gen_Buses_Number_Of_Elements; i++) {
			Bus_Number = GenBusNo[i] - 1; // -1 because we are using 0 based indexing
			if (BusType[Bus_Number] != 1 && GenStatus[i] != 0) {

				// both Qmin and Qmax are inf
				if (isinf(Qmin[i]) && isinf(Qmax[i])) {
					if (Qmin[i] < 0) {
						Qmin[i] = -fabs((*PFSolStruct).Qg[i]); // negative inf
					}
					else {
						Qmin[i] = fabs((*PFSolStruct).Qg[i]); // positive inf
					}

					if (Qmax[i] < 0) {
						Qmax[i] = -fabs((*PFSolStruct).Qg[i]); // negative inf
					}
					else {
						Qmax[i] = fabs((*PFSolStruct).Qg[i]); // positive inf
					}
				}

				// only Qmin is inf
				else if (isinf(Qmin[i])) {
					if (Qmin[i] < 0) { // negaive inf
						Qmin[i] = -(fabs((*PFSolStruct).Qg[i]) + fabs(Qmax[i]));
					}
					else { // positive inf
						Qmin[i] = fabs((*PFSolStruct).Qg[i]) + fabs(Qmax[i]);
					}
				}

				// only Qmax is inf
				else if (isinf(Qmax[i])) {
					if (Qmax[i] < 0) { // negaive inf
						Qmax[i] = -(fabs((*PFSolStruct).Qg[i]) + fabs(Qmin[i]));
					}
					else { // positive inf
						Qmax[i] = fabs((*PFSolStruct).Qg[i]) + fabs(Qmin[i]);
					}
				}


			}
		}

		// Vector of total Qg at each bus 	
		Chunk_Size = round(OnlineGenData_Num_Of_Elements / OMP_Cores);
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i,Bus_Number) shared(Chunk_Size,Gen_Buses_Number_Of_Elements,GenBusNo,BusType,GenStatus,Qmin,Qmax,PFSolStruct) reduction(+:Qg_total[:Num_of_Buses],Qg_min[:Num_of_Buses],Qg_max[:Num_of_Buses]) if(OMP_Enbale == 1)
		for (i = 0; i < Gen_Buses_Number_Of_Elements; i++) {
			Bus_Number = GenBusNo[i] - 1; // -1 because we are using 0 based indexing
			if (BusType[Bus_Number] != 1 && GenStatus[i] != 0) {
				Qg_total[Bus_Number] = Qg_total[Bus_Number] + (*PFSolStruct).Qg[i];
				Qg_min[Bus_Number] = Qg_min[Bus_Number] + Qmin[i];
				Qg_max[Bus_Number] = Qg_max[Bus_Number] + Qmax[i];
			}
		}

		// Divide proportionally
		for (i = 0; i < Gen_Buses_Number_Of_Elements; i++) {
			Bus_Number = GenBusNo[i] - 1; // -1 because we are using 0 based indexing
			if (BusType[Bus_Number] != 1 && GenStatus[i] != 0) {
				(*PFSolStruct).Qg[i] = Qmin[i] + ((Qg_total[Bus_Number] - Qg_min[Bus_Number]) / (Qg_max[Bus_Number] - Qg_min[Bus_Number] + eps)) * (Qmax[i] - Qmin[i]);				
			}
		}

				
		// Fix gens at buses with Qg range = 0 (use equal violation for all)
		for (i = 0; i < Gen_Buses_Number_Of_Elements; i++) {
			Bus_Number = GenBusNo[i] - 1;
			// buses with Qg range = 0
			if (BusType[Bus_Number] != 1 && GenStatus[i] != 0) {
				if (fabs(Qg_min[Bus_Number] - Qg_max[Bus_Number]) < 10 * eps) {
					// total mismatch @ bus div by number of gens
					mis = (Qg_total[Bus_Number] - Qg_min[Bus_Number]) / GenCount[Bus_Number];					
					(*PFSolStruct).Qg[i] = Qmin[i] + mis;
				}
			}
		}
		
		
	}

	
	// update Pg for slack gen(s)
	for (k = 0; k < Ref_Buses_Size; k++) {
		refgen = Ref_Buses[k];
		for (i = 0; i < Num_Online_Gens_No_PQ; i++) {
			Bus_Number = GenBus[i];
			if (Bus_Number == refgen) { // if Bus_Number = RefBuses(k)
				// inj P + local Pd
				temp = Sbus[i].real * BaseMVA + Pd_gbus[i];
				for (j = 0; j < Gen_Buses_Number_Of_Elements; j++) {
					// change the Pg value for the first generator at slack bus
					if (GenBusNo[j] == Bus_Number && GenStatus[j] != 0) {
						(*PFSolStruct).Pg[j] = temp;
						dummy = j;
						break;
					}
				}

				// more than one generator at this ref bus?
				if (GenCount[Bus_Number - 1] > 1) {
					// subtract off what is generated by other gens at this bus
					for (j = dummy + 1; j < Gen_Buses_Number_Of_Elements; j++) {
						if (GenBusNo[j] == Bus_Number && GenStatus[j] != 0) (*PFSolStruct).Pg[dummy] = (*PFSolStruct).Pg[dummy] - (*PFSolStruct).Pg[j];
					}
				}
				break;
			}
		}

	}

	timer_stop[1] = dsecnd(); // Generation calculation time
	(*PFSolStruct).Timer[1] = (timer_stop[1] - timer_start[1]);

	/*for (i = 0; i < Gen_Buses_Number_Of_Elements; i++) {		
			printf("%f\n", (*PFSolStruct).Qg[i]);		
	}*/

	
	/*for (i = 0; i < Gen_Buses_Number_Of_Elements; i++) {
		if (GenStatus[i] != 0) {
			printf("%f\n", (*PFSolStruct).Pg[i]);
		}
	}*/
	


	/* *********************  Calculate branch power flows ****************** */
	timer_start[2] = dsecnd(); // Branch flow calculation time
	// Count the number of out-of-service branches
	for (i = 0; i < Num_of_Branches; i++) {
		if (BranchStatus[i] == 0) 	Num_Offline_Branches++;
	}
	Num_Online_Branches = Num_of_Branches - Num_Offline_Branches;

	
	// ***** if all branches are in-service ****
	if (Num_Offline_Branches == 0) {
		Sf = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(Sf);
		St = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches, sizeof(MKL_Complex16), ALIGN);		
		CHECK_MEM(St);

		// YfV = Yf * V_rect
		// ******* Dense calculations ********	
		if (Sparse_Enable == 0) {
			cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Num_of_Branches, 1, Num_of_Buses, &alpha, (*YbusStruct).Yf_Dense, Num_of_Buses, V_rect, 1, &beta, Sf, 1);
			cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Num_of_Branches, 1, Num_of_Buses, &alpha, (*YbusStruct).Yt_Dense, Num_of_Buses, V_rect, 1, &beta, St, 1);
		}
		// ******* Sparse calculations ********
		else {
			CHECK_SPARSE(mkl_sparse_z_create_csr(&Yf_CSR, SPARSE_INDEX_BASE_ZERO, Num_of_Branches, Num_of_Buses, (*YbusStruct).Yf_JA, (*YbusStruct).Yf_JA + 1, (*YbusStruct).Yf_IA, (*YbusStruct).Yf_Val));
			CHECK_SPARSE(mkl_sparse_z_create_csr(&Yt_CSR, SPARSE_INDEX_BASE_ZERO, Num_of_Branches, Num_of_Buses, (*YbusStruct).Yt_JA, (*YbusStruct).Yt_JA + 1, (*YbusStruct).Yt_IA, (*YbusStruct).Yt_Val));
			CHECK_SPARSE(mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, Yf_CSR, descr, V_rect, beta, Sf));
			CHECK_SPARSE(mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, Yt_CSR, descr, V_rect, beta, St));
		}

		for (i = 0; i < Num_of_Branches; i++) {
			Bus_Number = FromBus[i] - 1;
			(*PFSolStruct).PF[i] = Mult_re(V_rect[Bus_Number].real, V_rect[Bus_Number].imag, Sf[i].real, -1 * Sf[i].imag) * BaseMVA;
			(*PFSolStruct).QF[i] = Mult_img(V_rect[Bus_Number].real, V_rect[Bus_Number].imag, Sf[i].real, -1 * Sf[i].imag) * BaseMVA;
			Bus_Number = ToBus[i] - 1;
			(*PFSolStruct).PT[i] = Mult_re(V_rect[Bus_Number].real, V_rect[Bus_Number].imag, St[i].real, -1 * St[i].imag) * BaseMVA;
			(*PFSolStruct).QT[i] = Mult_img(V_rect[Bus_Number].real, V_rect[Bus_Number].imag, St[i].real, -1 * St[i].imag) * BaseMVA;
		}

	}
	// if there are some out-of-service branches
	else {
		
		Sf = (MKL_Complex16 *)mkl_calloc((size_t)Num_Online_Branches, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(Sf);
		St = (MKL_Complex16 *)mkl_calloc((size_t)Num_Online_Branches, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(St);

		// YfV = Yf * V_rect
		// ******* Dense calculations ********	
		if (Sparse_Enable == 0) {
			Yf_New = (MKL_Complex16 *)mkl_calloc((size_t)(Num_Online_Branches*Num_of_Buses), sizeof(MKL_Complex16), ALIGN);
			CHECK_MEM(Yf_New);
			Yt_New = (MKL_Complex16 *)mkl_calloc((size_t)(Num_Online_Branches*Num_of_Buses), sizeof(MKL_Complex16), ALIGN);
			CHECK_MEM(Yt_New);
			// Copy Yf & Yt of online branches
			k = 0;
			Chunk_Size = round(Num_of_Branches / OMP_Cores);
			for (i = 0; i < Num_of_Buses; i++) {
				if (BranchStatus[i] != 0) {
					#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(j) shared(Chunk_Size,Num_of_Buses,k,i,Yf_New,Yt_New,YbusStruct) if(OMP_Enbale == 1)
					for (j = 0; j < Num_of_Buses; j++) {
						Yf_New[k*Num_of_Buses + j].real = (*YbusStruct).Yf_Dense[i*Num_of_Buses + j].real;
						Yf_New[k*Num_of_Buses + j].imag = (*YbusStruct).Yf_Dense[i*Num_of_Buses + j].imag;
						Yt_New[k*Num_of_Buses + j].real = (*YbusStruct).Yt_Dense[i*Num_of_Buses + j].real;
						Yt_New[k*Num_of_Buses + j].imag = (*YbusStruct).Yt_Dense[i*Num_of_Buses + j].imag;
					}
					k++;
				}				
			}			
			cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Num_Online_Branches, 1, Num_of_Buses, &alpha, Yf_New, Num_of_Buses, V_rect, 1, &beta, Sf, 1);
			cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Num_Online_Branches, 1, Num_of_Buses, &alpha, Yt_New, Num_of_Buses, V_rect, 1, &beta, St, 1);			
		}
		// ******* Sparse calculations ********
		else {
			
			// Count the number of nonzero elements in the Yf, Yt matrices of online branches
			nnz = 0;
			Chunk_Size = round(Num_of_Branches / OMP_Cores);
			#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i,Bus_Number) shared(Chunk_Size,Num_of_Branches,BranchStatus,YbusStruct) reduction(+:nnz) if(OMP_Enbale == 1)
			for (i = 0; i < Num_of_Branches; i++) {
				if (BranchStatus[i] != 0) 	nnz = nnz + ((*YbusStruct).Yt_JA[i + 1] - (*YbusStruct).Yt_JA[i]);				
			}
					
			Yt_New = (MKL_Complex16 *)mkl_calloc((size_t)(nnz), sizeof(MKL_Complex16), ALIGN);
			CHECK_MEM(Yt_New);
			Yt_New_IA = (MKL_INT *)mkl_calloc((size_t)(nnz), sizeof(MKL_INT), ALIGN);
			CHECK_MEM(Yt_New_IA);
			Yt_New_JA = (MKL_INT *)mkl_calloc((size_t)(Num_Online_Branches + 1), sizeof(MKL_INT), ALIGN);			
			CHECK_MEM(Yt_New_JA);
			
			Yf_New = (MKL_Complex16 *)mkl_calloc((size_t)(nnz), sizeof(MKL_Complex16), ALIGN);
			CHECK_MEM(Yf_New);
			Yf_New_IA = (MKL_INT *)mkl_calloc((size_t)(nnz), sizeof(MKL_INT), ALIGN);
			CHECK_MEM(Yf_New_IA);
			Yf_New_JA = (MKL_INT *)mkl_calloc((size_t)(Num_Online_Branches + 1), sizeof(MKL_INT), ALIGN);						
			CHECK_MEM(Yf_New_JA);
			
			// Copy Yf & Yt of online branches
			Yt_New_JA[0] = 0;
			Yf_New_JA[0] = 0;

			k = 0; m = 0;  nnz = 0;
			for (i = 0; i < Num_of_Branches; i++) {
				if (BranchStatus[i] != 0) {					
					for (j = (*YbusStruct).Yt_JA[i]; j < (*YbusStruct).Yt_JA[i + 1]; j++) {
						Yt_New[k].real = (*YbusStruct).Yt_Val[j].real;
						Yt_New[k].imag = (*YbusStruct).Yt_Val[j].imag;
						Yt_New_IA[k] = (*YbusStruct).Yt_IA[j];

						Yf_New[k].real = (*YbusStruct).Yf_Val[j].real;
						Yf_New[k].imag = (*YbusStruct).Yf_Val[j].imag;										  						
						Yf_New_IA[k] = (*YbusStruct).Yf_IA[j];

						k++;
					}
					m++;
					nnz = (*YbusStruct).Yt_JA[i + 1] - (*YbusStruct).Yt_JA[i];
					Yf_New_JA[m] = Yf_New_JA[m - 1] + nnz;
					Yt_New_JA[m] = Yt_New_JA[m - 1] + nnz;
				}
			}

			CHECK_SPARSE(mkl_sparse_z_create_csr(&Yf_CSR, SPARSE_INDEX_BASE_ZERO, Num_Online_Branches, Num_of_Buses, Yf_New_JA, Yf_New_JA + 1, Yf_New_IA, Yf_New));
			CHECK_SPARSE(mkl_sparse_z_create_csr(&Yt_CSR, SPARSE_INDEX_BASE_ZERO, Num_Online_Branches, Num_of_Buses, Yt_New_JA, Yt_New_JA + 1, Yt_New_IA, Yt_New));
			CHECK_SPARSE(mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, Yf_CSR, descr, V_rect, beta, Sf));
			CHECK_SPARSE(mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, Yt_CSR, descr, V_rect, beta, St));
		}

		k = 0;
		for (i = 0; i < Num_of_Branches; i++) {
			if (BranchStatus[i] != 0) {
				Bus_Number = FromBus[i] - 1;
				(*PFSolStruct).PF[i] = Mult_re(V_rect[Bus_Number].real, V_rect[Bus_Number].imag, Sf[k].real, -1 * Sf[k].imag) * BaseMVA;
				(*PFSolStruct).QF[i] = Mult_img(V_rect[Bus_Number].real, V_rect[Bus_Number].imag, Sf[k].real, -1 * Sf[k].imag) * BaseMVA;
				Bus_Number = ToBus[i] - 1;
				(*PFSolStruct).PT[i] = Mult_re(V_rect[Bus_Number].real, V_rect[Bus_Number].imag, St[k].real, -1 * St[k].imag) * BaseMVA;
				(*PFSolStruct).QT[i] = Mult_img(V_rect[Bus_Number].real, V_rect[Bus_Number].imag, St[k].real, -1 * St[k].imag) * BaseMVA;
				k++;
			}
			else {
				(*PFSolStruct).PF[i] = 0;
				(*PFSolStruct).QF[i] = 0;
				(*PFSolStruct).PT[i] = 0;
				(*PFSolStruct).QT[i] = 0;
			}
			//printf("PF[%i] = %f, Qf = %f, PT = %f , QT = %f \n", i, (*PFSolStruct).PF[i], (*PFSolStruct).QF[i], (*PFSolStruct).PT[i], (*PFSolStruct).QT[i]);
		}

	}
	timer_stop[2] = dsecnd(); // Branch flow calculation time
	(*PFSolStruct).Timer[2] = (timer_stop[2] - timer_start[2]);


	/* ************************** Calculate line losses ************************* */
	timer_start[3] = dsecnd(); // Line losses calculation time
	// loss = baseMVA * Ysc.*Vdrop.*conj(Vdrop)
	// Vdrop = A * V : vector of voltage drop across series impedance element		
	
	// Vdrop = A * V
	if (Sparse_Enable == 0) {		
		// ******* Dense calculations ********					
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Num_of_Branches, 1, Num_of_Buses, &alpha, (*YbusStruct).A_Dense, Num_of_Buses, V_rect, 1, &beta, Vdrop, 1); 			
	}
	else {
		// ******* Sparse calculations ********		
		CHECK_SPARSE(mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, (*YbusStruct).A_CSR_Handle, descr, V_rect, beta, Vdrop));
	}
	
	// Calculate losses
	Chunk_Size = round(Num_of_Branches / OMP_Cores);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i,temp) shared(Chunk_Size,Num_of_Branches,Vdrop,BaseMVA,YbusStruct,PFSolStruct) if(OMP_Enbale == 1)
	for (i = 0; i < Num_of_Branches; i++) {		
		temp = Vdrop[i].real*Vdrop[i].real + Vdrop[i].imag*Vdrop[i].imag; // Vdrop .* conj(Vdrop);
		(*PFSolStruct).Loss[i].real = BaseMVA * (*YbusStruct).Ysc[i].real * temp;
		(*PFSolStruct).Loss[i].imag = BaseMVA * (*YbusStruct).Ysc[i].imag * temp;
	}

	timer_stop[3] = dsecnd(); // Line losses calculation time
	(*PFSolStruct).Timer[3] = (timer_stop[3] - timer_start[3]);


	/* ********************* Calculate line charging injections ********************* */
	timer_start[4] = dsecnd(); // Line charging calculation time
	if (Sparse_Enable == 0) {
		// ******* Dense calculations ********			
		cblas_zgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Num_of_Branches, 1, Num_of_Buses, &alpha, (*YbusStruct).CfPrime_Dense, Num_of_Branches, V_rect, 1, &beta, Vf, 1);
		cblas_zgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Num_of_Branches, 1, Num_of_Buses, &alpha, (*YbusStruct).CtPrime_Dense, Num_of_Branches, V_rect, 1, &beta, Vt, 1);
	}
	else {
		// ******* Sparse calculations ********		
		CHECK_SPARSE(mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, (*YbusStruct).Cf_CSR_Handle, descr, V_rect, beta, Vf));
		CHECK_SPARSE(mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, (*YbusStruct).Ct_CSR_Handle, descr, V_rect, beta, Vt));
	}

	Chunk_Size = round(Num_of_Branches / OMP_Cores);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Num_of_Branches,Sparse_Enable,BranchStatus,BaseMVA,b,Vt,Vf,YbusStruct,PFSolStruct) if(OMP_Enbale == 1)
	for (i = 0; i < Num_of_Branches; i++) {
		(*PFSolStruct).TCHG[i] = (BaseMVA / 2) * b[i] * (Vt[i].real*Vt[i].real + Vt[i].imag*Vt[i].imag) * BranchStatus[i];
		// Dense
		if (Sparse_Enable == 0) {
			(*PFSolStruct).FCHG[i] = ((BaseMVA / 2) * b[i] * (Vf[i].real*Vf[i].real + Vf[i].imag*Vf[i].imag) * BranchStatus[i]) / ((*YbusStruct).Tap_Dense[i*Num_of_Branches + i].real*(*YbusStruct).Tap_Dense[i*Num_of_Branches + i].real + (*YbusStruct).Tap_Dense[i*Num_of_Branches + i].imag*(*YbusStruct).Tap_Dense[i*Num_of_Branches + i].imag);
		}
		// Sparse
		else {
			(*PFSolStruct).FCHG[i] = ((BaseMVA / 2) * b[i] * (Vf[i].real*Vf[i].real + Vf[i].imag*Vf[i].imag) * BranchStatus[i]) / ((*YbusStruct).Tap_Val[i].real*(*YbusStruct).Tap_Val[i].real + (*YbusStruct).Tap_Val[i].imag*(*YbusStruct).Tap_Val[i].imag);
		}
				
	}

	timer_stop[4] = dsecnd(); // Line charging calculation time
	(*PFSolStruct).Timer[4] = (timer_stop[4] - timer_start[4]);


	/* ********************* Partial derivatives ********************* */	
	timer_start[5] = dsecnd(); // Partial derivatives calculation time
	if (Sparse_Enable == 0) {
		// ******* Dense calculations ********					
		
		B = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches*Num_of_Buses, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(B);
		B_conj_B = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches*Num_of_Buses, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(B_conj_B);
		V_drop_diag = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches*Num_of_Branches, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(V_drop_diag);
		VddA = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches*Num_of_Branches, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(VddA);
		dYsc = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches*Num_of_Branches, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(dYsc);
		Bc_diag = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches*Num_of_Branches, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(Bc_diag);
		Tap_diag = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches*Num_of_Branches, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(Tap_diag);
		dlosstemp = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches*Num_of_Buses, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(dlosstemp);
		
		V_mag = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Buses, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(V_mag);
		Vm_diag_recip = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Buses*Num_of_Buses, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(Vm_diag_recip);
		Vrect_conj_diag = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Buses*Num_of_Buses, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(Vrect_conj_diag);
				
		BcTap = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches*Num_of_Branches, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(BcTap);
		CfVm = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(CfVm);
		CfVm_diag = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches*Num_of_Branches, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(CfVm_diag);
		BcTapCfVm = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches*Num_of_Branches, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(BcTapCfVm);

					
		CtVm = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(CtVm);
		CtVm_diag = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches*Num_of_Branches, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(CtVm_diag);
		BcCtVm = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches*Num_of_Branches, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(BcCtVm);

		
		(*PFSolStruct).dloss_dV_a_Dense = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches*Num_of_Buses, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM((*PFSolStruct).dloss_dV_a_Dense);
		(*PFSolStruct).dloss_dV_m_Dense = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches*Num_of_Buses, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM((*PFSolStruct).dloss_dV_m_Dense);
		(*PFSolStruct).dchg_dVm_f_Dense = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches*Num_of_Buses, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM((*PFSolStruct).dchg_dVm_f_Dense);
		(*PFSolStruct).dchg_dVm_t_Dense = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches*Num_of_Buses, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM((*PFSolStruct).dchg_dVm_t_Dense);

		
		// Create diag(Vdrop), diag(Ysc), diag(Bc), diag(1./TapTap*)
		Chunk_Size = round(Num_of_Branches / OMP_Cores);
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Num_of_Branches,BaseMVA,V_drop_diag,Vdrop,dYsc,Bc_diag,b,Tap_diag,YbusStruct) if(OMP_Enbale == 1)
		for (i = 0; i < Num_of_Branches; i++) {
			V_drop_diag[i*Num_of_Branches + i].real = Vdrop[i].real;
			V_drop_diag[i*Num_of_Branches + i].imag = Vdrop[i].imag;
			dYsc[i*Num_of_Branches + i].real = (*YbusStruct).Ysc[i].real * BaseMVA;
			dYsc[i*Num_of_Branches + i].imag = (*YbusStruct).Ysc[i].imag * BaseMVA;
			Bc_diag[i*Num_of_Branches + i].real = b[i] * BaseMVA;
			Tap_diag[i*Num_of_Branches + i].real = 1 / ((*YbusStruct).Tap_Dense[i*Num_of_Branches + i].real*(*YbusStruct).Tap_Dense[i*Num_of_Branches + i].real + (*YbusStruct).Tap_Dense[i*Num_of_Branches + i].imag*(*YbusStruct).Tap_Dense[i*Num_of_Branches + i].imag);
		}
		
		// Create diag(conj(V)) and diag(1/Vm)
		Chunk_Size = round(Num_of_Buses / OMP_Cores);
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Num_of_Buses,V_mag,Vm_cal,Vm_diag_recip,Vrect_conj_diag,V_rect) if(OMP_Enbale == 1)
		for (i = 0; i < Num_of_Buses; i++) {
			V_mag[i].real = Vm_cal[i];
			Vm_diag_recip[i*Num_of_Buses + i].real = 1 / Vm_cal[i];
			Vrect_conj_diag[i*Num_of_Buses + i].real = V_rect[i].real;
			Vrect_conj_diag[i*Num_of_Buses + i].imag = -1*V_rect[i].imag;
		}
		
		/* ********* dloss_dV ********* */
		/* B = diag(A * V = Vdrop) * conj(A) * diag(conj(V)) */
		// VddA = diag(Vdrop)*conj(A) : The imaginary part of A is zero, therefore conj(A) = A				
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Num_of_Branches, Num_of_Buses, Num_of_Branches, &alpha, V_drop_diag, Num_of_Branches, (*YbusStruct).A_Dense, Num_of_Buses, &beta, VddA, Num_of_Buses);
		// B = VddA * Vrect_conj_diag
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Num_of_Branches, Num_of_Buses, Num_of_Buses, &alpha, VddA, Num_of_Buses, Vrect_conj_diag, Num_of_Buses, &beta, B, Num_of_Buses);
		// B_conj_B = B - conj(B) = 0 + 2*imag(B)
		Chunk_Size = round(Num_of_Buses*Num_of_Branches / OMP_Cores);
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Num_of_Branches,Num_of_Buses,B,B_conj_B) if(OMP_Enbale == 1)
		for (i = 0; i < Num_of_Branches*Num_of_Buses; i++) {
			B_conj_B[i].real = 0;
			B_conj_B[i].imag = B[i].imag*2;
		}		
		// dloss_dV_a = -1j * baseMVA * dYsc * (B - conj(B)) = -1j * dYsc * B_conj_B
		alpha.real = 0;
		alpha.imag = -1;
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Num_of_Branches, Num_of_Buses, Num_of_Branches, &alpha, dYsc, Num_of_Branches, B_conj_B, Num_of_Buses, &beta, (*PFSolStruct).dloss_dV_a_Dense, Num_of_Buses);		
		// dloss_dV_m = baseMVA * dYsc * (B + conj(B)) * diag(1/Vm) = dYsc * B_conj_B * Vm_diag_recip
		Chunk_Size = round(Num_of_Buses*Num_of_Branches / OMP_Cores);
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Num_of_Branches,Num_of_Buses,B,B_conj_B) if(OMP_Enbale == 1)
		for (i = 0; i < Num_of_Branches*Num_of_Buses; i++) {
			B_conj_B[i].real = B[i].real * 2;
			B_conj_B[i].imag = 0;
		}
		alpha.real = 1;
		alpha.imag = 0;
		// dlosstemp = dYsc * B_conj_B 
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Num_of_Branches, Num_of_Buses, Num_of_Branches, &alpha, dYsc, Num_of_Branches, B_conj_B, Num_of_Buses, &beta, dlosstemp, Num_of_Buses);
		// dloss_dV_m = dlosstemp * Vm_diag_recip
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Num_of_Branches, Num_of_Buses, Num_of_Buses, &alpha, dlosstemp, Num_of_Buses, Vm_diag_recip, Num_of_Buses, &beta, (*PFSolStruct).dloss_dV_m_Dense, Num_of_Buses);


		/* ********* dchg_dVm ********* */
		/* dchg_dVm_f = baseMVA * Bc * diag(1/Tap.conj(Tap)) * diag(Cf * VM) * Cf */
		/* dchg_dVm_t = baseMVA * Bc * diag(Ct * VM) * Ct */
		
		// CfVm = Cf * V_mag
		cblas_zgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Num_of_Branches, 1, Num_of_Buses, &alpha, (*YbusStruct).CfPrime_Dense, Num_of_Branches, V_mag, 1, &beta, CfVm, 1);
		// CtVm = Ct * V_mag
		cblas_zgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Num_of_Branches, 1, Num_of_Buses, &alpha, (*YbusStruct).CtPrime_Dense, Num_of_Branches, V_mag, 1, &beta, CtVm, 1);
		// Create CfVm_diag and CtVm_diag
		Chunk_Size = round(Num_of_Branches / OMP_Cores);
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Num_of_Branches,CfVm_diag,CtVm_diag,CfVm,CtVm) if(OMP_Enbale == 1)
		for (i = 0; i < Num_of_Branches; i++) {
			CfVm_diag[i*Num_of_Branches + i].real = CfVm[i].real;
			CfVm_diag[i*Num_of_Branches + i].imag = CfVm[i].imag;	
			CtVm_diag[i*Num_of_Branches + i].real = CtVm[i].real;
			CtVm_diag[i*Num_of_Branches + i].imag = CtVm[i].imag;
		}
		// BcTap = baseMVA * Bc * diag(1/Tap.conj(Tap)) = Bc * Tap_diag
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Num_of_Branches, Num_of_Branches, Num_of_Branches, &alpha, Bc_diag, Num_of_Branches, Tap_diag, Num_of_Branches, &beta, BcTap, Num_of_Branches);
		// BcTapCfVm = BcTap * CfVm_diag
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Num_of_Branches, Num_of_Branches, Num_of_Branches, &alpha, BcTap, Num_of_Branches, CfVm_diag, Num_of_Branches, &beta, BcTapCfVm, Num_of_Branches);
		// BcCtVm = Bc_diag * CtVm_diag
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Num_of_Branches, Num_of_Branches, Num_of_Branches, &alpha, Bc_diag, Num_of_Branches, CtVm_diag, Num_of_Branches, &beta, BcCtVm, Num_of_Branches);

		// dchg_dVm_f = BcTapCfVm * Cf
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasTrans, Num_of_Branches, Num_of_Buses, Num_of_Branches, &alpha, BcTapCfVm, Num_of_Branches, (*YbusStruct).CfPrime_Dense, Num_of_Branches, &beta, (*PFSolStruct).dchg_dVm_f_Dense, Num_of_Buses);
		// dchg_dVm_t = BcCtVm * Ct
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasTrans, Num_of_Branches, Num_of_Buses, Num_of_Branches, &alpha, BcCtVm, Num_of_Branches, (*YbusStruct).CtPrime_Dense, Num_of_Branches, &beta, (*PFSolStruct).dchg_dVm_t_Dense, Num_of_Buses);
	}
	else {
		// ******* Sparse calculations ********	

		IA_Bus = (MKL_INT *)mkl_calloc((size_t)(Num_of_Buses + 1), sizeof(MKL_INT), ALIGN);
		CHECK_MEM(IA_Bus);
		JA_Bus = (MKL_INT *)mkl_calloc((size_t)(Num_of_Buses + 1), sizeof(MKL_INT), ALIGN);
		CHECK_MEM(JA_Bus);
				
		IA_Br = (MKL_INT *)mkl_calloc((size_t)(Num_of_Branches + 1), sizeof(MKL_INT), ALIGN);
		CHECK_MEM(IA_Br);
		JA_Br = (MKL_INT *)mkl_calloc((size_t)(Num_of_Branches + 1), sizeof(MKL_INT), ALIGN);
		CHECK_MEM(JA_Br);
		

		dYsc_val = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(dYsc_val);
		dYsc_mod_val = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches, sizeof(MKL_Complex16), ALIGN); // -1j * dYsc
		CHECK_MEM(dYsc_mod_val);
		Bc_diag_val = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(Bc_diag_val);
		Tap_diag_val = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(Tap_diag_val);
		CfVm = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(CfVm);
		CtVm = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Branches, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(CtVm);

		V_mag_val = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Buses, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(V_mag_val);
		Vm_recip_val = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Buses, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(Vm_recip_val);
		Vrect_conj_val = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Buses, sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(Vrect_conj_val);

		// Create IA and JA for diagonal matrices with size row=col= Branch or Bus
		Chunk_Size = round(Num_of_Buses / OMP_Cores);
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Num_of_Buses,IA_Bus,JA_Bus) if(OMP_Enbale == 1)
		for (i = 0; i < Num_of_Buses; i++) {
			IA_Bus[i] = i;
			JA_Bus[i] = i;
		}
		JA_Bus[Num_of_Buses] = Num_of_Buses;
		
		Chunk_Size = round(Num_of_Branches / OMP_Cores);
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Num_of_Branches,IA_Br,JA_Br) if(OMP_Enbale == 1)
		for (i = 0; i < Num_of_Branches; i++) {
			IA_Br[i] = i;
			JA_Br[i] = i;
		}
		JA_Br[Num_of_Branches] = Num_of_Branches;

				
		// Create diag(Ysc), diag(Bc), diag(1. / TapTap * )
		Chunk_Size = round(Num_of_Branches / OMP_Cores);
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Num_of_Branches,BaseMVA,dYsc_val,dYsc_mod_val,Bc_diag_val,Tap_diag_val,b,YbusStruct) if(OMP_Enbale == 1)
		for (i = 0; i < Num_of_Branches; i++) {
			dYsc_val[i].real = (*YbusStruct).Ysc[i].real * BaseMVA;
			dYsc_val[i].imag = (*YbusStruct).Ysc[i].imag * BaseMVA;
			dYsc_mod_val[i].real = dYsc_val[i].imag; // -1j * dYsc
			dYsc_mod_val[i].imag = -1*dYsc_val[i].real;
			Bc_diag_val[i].real = b[i] * BaseMVA;
			Tap_diag_val[i].real = 1 / ((*YbusStruct).Tap_Val[i].real*(*YbusStruct).Tap_Val[i].real + (*YbusStruct).Tap_Val[i].imag*(*YbusStruct).Tap_Val[i].imag);
		}

		// Create diag(conj(V)) and diag(1/Vm)
		Chunk_Size = round(Num_of_Buses / OMP_Cores);
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Num_of_Buses,V_mag_val,Vm_recip_val,Vrect_conj_val,V_rect,Vm_cal) if(OMP_Enbale == 1)
		for (i = 0; i < Num_of_Buses; i++) {
			V_mag_val[i].real = Vm_cal[i];
			Vm_recip_val[i].real = 1 / Vm_cal[i];
			Vrect_conj_val[i].real = V_rect[i].real;
			Vrect_conj_val[i].imag = -1 * V_rect[i].imag;			
		}
		
		// Creat CSR handle for all diagonal matrices
		CHECK_SPARSE(mkl_sparse_z_create_csr(&dYsc_CSR_Handle, SPARSE_INDEX_BASE_ZERO, Num_of_Branches, Num_of_Branches, JA_Br, JA_Br + 1, IA_Br, dYsc_val));
		CHECK_SPARSE(mkl_sparse_z_create_csr(&dYsc_mod_CSR_Handle, SPARSE_INDEX_BASE_ZERO, Num_of_Branches, Num_of_Branches, JA_Br, JA_Br + 1, IA_Br, dYsc_mod_val));
		CHECK_SPARSE(mkl_sparse_z_create_csr(&Bc_CSR_Handle, SPARSE_INDEX_BASE_ZERO, Num_of_Branches, Num_of_Branches, JA_Br, JA_Br + 1, IA_Br, Bc_diag_val));
		CHECK_SPARSE(mkl_sparse_z_create_csr(&Tap_CSR_Handle, SPARSE_INDEX_BASE_ZERO, Num_of_Branches, Num_of_Branches, JA_Br, JA_Br + 1, IA_Br, Tap_diag_val));
		CHECK_SPARSE(mkl_sparse_z_create_csr(&Vdrop_CSR_Handle, SPARSE_INDEX_BASE_ZERO, Num_of_Branches, Num_of_Branches, JA_Br, JA_Br + 1, IA_Br, Vdrop));

		CHECK_SPARSE(mkl_sparse_z_create_csr(&Vm_reci_CSR_Handle, SPARSE_INDEX_BASE_ZERO, Num_of_Buses, Num_of_Buses, JA_Bus, JA_Bus + 1, IA_Bus, Vm_recip_val));
		CHECK_SPARSE(mkl_sparse_z_create_csr(&Vconj_CSR_Handle, SPARSE_INDEX_BASE_ZERO, Num_of_Buses, Num_of_Buses, JA_Bus, JA_Bus + 1, IA_Bus, Vrect_conj_val));
		

		/* ********* dloss_dV ********* */
		/* B = diag(A * V = Vdrop) * conj(A) * diag(conj(V)) */
		// VddA = diag(Vdrop)*conj(A) : The imaginary part of A is zero, therefore conj(A) = A						
		CHECK_SPARSE(mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, Vdrop_CSR_Handle, (*YbusStruct).A_CSR_Handle, &VddA_CSR_Handle));				
		// B = VddA * Vrect_conj_diag
		CHECK_SPARSE(mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, VddA_CSR_Handle, Vconj_CSR_Handle, &B_CSR_Handle));	
		// Extract the values: will be used in calculating B-conj(B) and B+conj(B)
		CHECK_SPARSE(mkl_sparse_z_export_csr(B_CSR_Handle, &indexing, &rows, &cols, &Dummy_JA_1, &Dummy_pointerE_1, &Dummy_IA_1, &B_m_conjB_val));		
		B_p_conjB_val = (MKL_Complex16 *)mkl_calloc((size_t)Dummy_JA_1[rows], sizeof(MKL_Complex16), ALIGN);			
		// B_m_conjB_val = B - conj(B) = 0 + 2*imag(B)
		// B_p_conjB_val = B + conj(B) = 2*real(B) + 0
		Chunk_Size = round(Dummy_JA_1[rows] / OMP_Cores);
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Dummy_JA_1,rows,B_p_conjB_val,B_m_conjB_val) if(OMP_Enbale == 1)
		for (i = 0; i < Dummy_JA_1[rows]; i++) {
			B_p_conjB_val[i].real = B_m_conjB_val[i].real * 2;
			B_m_conjB_val[i].real = 0;			
			B_m_conjB_val[i].imag = B_m_conjB_val[i].imag * 2;			
		}
		// Create B-conj(B) and B+conj(B) handles
		CHECK_SPARSE(mkl_sparse_z_create_csr(&B_m_conjB_CSR_Handle, SPARSE_INDEX_BASE_ZERO, rows, cols, Dummy_JA_1, Dummy_JA_1 + 1, Dummy_IA_1, B_m_conjB_val));
		CHECK_SPARSE(mkl_sparse_z_create_csr(&B_p_conjB_CSR_Handle, SPARSE_INDEX_BASE_ZERO, rows, cols, Dummy_JA_1, Dummy_JA_1 + 1, Dummy_IA_1, B_p_conjB_val));
			
		// dloss_dV_a = -1j * baseMVA * dYsc * (B - conj(B)) = dYsc_mod * B_m_conjB_val
		CHECK_SPARSE(mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, dYsc_mod_CSR_Handle, B_m_conjB_CSR_Handle, &(*PFSolStruct).dloss_dV_a_CSR_Handle));
		CHECK_SPARSE(mkl_sparse_order((*PFSolStruct).dloss_dV_a_CSR_Handle));
		CHECK_SPARSE(mkl_sparse_z_export_csr((*PFSolStruct).dloss_dV_a_CSR_Handle, &indexing2, &rows2, &cols2, &(*PFSolStruct).dloss_dV_a_JA, &(*PFSolStruct).dloss_dV_a_pointerE, &(*PFSolStruct).dloss_dV_a_IA, &(*PFSolStruct).dloss_dV_a_val));

		// dloss_dV_m = baseMVA * dYsc * (B + conj(B)) * diag(1/Vm) = dYsc * B_p_conj_B * Vm_diag_recip
		CHECK_SPARSE(mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, dYsc_CSR_Handle, B_p_conjB_CSR_Handle, &dYscBp_CSR_Handle));
		CHECK_SPARSE(mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, dYscBp_CSR_Handle, Vm_reci_CSR_Handle, &(*PFSolStruct).dloss_dV_m_CSR_Handle));
		CHECK_SPARSE(mkl_sparse_order((*PFSolStruct).dloss_dV_m_CSR_Handle));
		CHECK_SPARSE(mkl_sparse_z_export_csr((*PFSolStruct).dloss_dV_m_CSR_Handle, &indexing2, &rows2, &cols2, &(*PFSolStruct).dloss_dV_m_JA, &(*PFSolStruct).dloss_dV_m_pointerE, &(*PFSolStruct).dloss_dV_m_IA, &(*PFSolStruct).dloss_dV_m_val));


		/* ********* dchg_dVm ********* */
		/* dchg_dVm_f = BaseMVA * diag(Bc) * diag(1/Tap.*conj(Tap)) * diag(Cf * Vm) * Cf */
		// BcTap = diag(Bc) * diag(1/Tap)
		CHECK_SPARSE(mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, Bc_CSR_Handle, Tap_CSR_Handle, &BcTap_CSR_Handle));
		// CfVm = Cf * Vm
		CHECK_SPARSE(mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, (*YbusStruct).Cf_CSR_Handle, descr, V_mag_val, beta, CfVm));
		// Create diag(CfVm)
		CHECK_SPARSE(mkl_sparse_z_create_csr(&CfVm_CSR_Handle, SPARSE_INDEX_BASE_ZERO, Num_of_Branches, Num_of_Branches, JA_Br, JA_Br + 1, IA_Br, CfVm));
		// CfVmCf_CSR_Handle = diag(CfVm) * Cf
		CHECK_SPARSE(mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, CfVm_CSR_Handle, (*YbusStruct).Cf_CSR_Handle, &CfVmCf_CSR_Handle));
		// dchg_dVm_f_CSR_Handle = BcTap_CSR_Handle * CfVmCf_CSR_Handle
		CHECK_SPARSE(mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, BcTap_CSR_Handle, CfVmCf_CSR_Handle, &(*PFSolStruct).dchg_dVm_f_CSR_Handle));
		CHECK_SPARSE(mkl_sparse_z_export_csr((*PFSolStruct).dchg_dVm_f_CSR_Handle, &indexing2, &rows2, &cols2, &(*PFSolStruct).dchg_dVm_f_JA, &(*PFSolStruct).dchg_dVm_f_pointerE, &(*PFSolStruct).dchg_dVm_f_IA, &(*PFSolStruct).dchg_dVm_f_val));

		/* dchg_dVm_t = BaseMVA * diag(Bc) * diag(Ct * Vm) * Ct */
		// CfVm = Cf * Vm
		CHECK_SPARSE(mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, (*YbusStruct).Ct_CSR_Handle, descr, V_mag_val, beta, CtVm));
		// Create diag(CtVm)
		CHECK_SPARSE(mkl_sparse_z_create_csr(&CtVm_CSR_Handle, SPARSE_INDEX_BASE_ZERO, Num_of_Branches, Num_of_Branches, JA_Br, JA_Br + 1, IA_Br, CtVm));
		// CtVmCf_CSR_Handle = diag(CtVm) * Ct
		CHECK_SPARSE(mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, CtVm_CSR_Handle, (*YbusStruct).Ct_CSR_Handle, &CtVmCt_CSR_Handle));
		// dchg_dVm_t_CSR_Handle = BcTap_CSR_Handle * CfVmCf_CSR_Handle
		CHECK_SPARSE(mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, Bc_CSR_Handle, CtVmCt_CSR_Handle, &(*PFSolStruct).dchg_dVm_t_CSR_Handle));
		CHECK_SPARSE(mkl_sparse_z_export_csr((*PFSolStruct).dchg_dVm_t_CSR_Handle, &indexing2, &rows2, &cols2, &(*PFSolStruct).dchg_dVm_t_JA, &(*PFSolStruct).dchg_dVm_t_pointerE, &(*PFSolStruct).dchg_dVm_t_IA, &(*PFSolStruct).dchg_dVm_t_val));			

		//for (i = 0; i < (*PFSolStruct).dloss_dV_a_JA[Num_of_Branches]; i++) {
		//	printf("val(%i) = %f , %f \n", i, (*PFSolStruct).dloss_dV_a_val[i].real, (*PFSolStruct).dloss_dV_a_val[i].imag);
		//	printf("val(%i) = %f , %f \n", i, (*PFSolStruct).dloss_dV_m_val[i].real, (*PFSolStruct).dloss_dV_m_val[i].imag);
		//}
			
		//for (i = 0; i < (*PFSolStruct).dchg_dVm_f_JA[Num_of_Branches]; i++) {
		//	//printf("val(%i) = %f , %f \n", i, (*PFSolStruct).dchg_dVm_f_val[i].real, (*PFSolStruct).dchg_dVm_f_val[i].imag);
		//	printf("val(%i) = %f , %f \n", i, (*PFSolStruct).dchg_dVm_t_val[i].real, (*PFSolStruct).dchg_dVm_t_val[i].imag);
		//}
				
	}
	
	timer_stop[5] = dsecnd(); // Partial derivatives calculation time
	(*PFSolStruct).Timer[5] = (timer_stop[5] - timer_start[5]);

	/* ****************** General information ****************** */
	timer_start[6] = dsecnd(); // Misc. calculation time
	for (i = 0; i < Num_of_Branches; i++) {
		// Total losses		
		(*PFSolStruct).P_loss_total = (*PFSolStruct).P_loss_total + (*PFSolStruct).Loss[i].real;
		(*PFSolStruct).Q_loss_total = (*PFSolStruct).Q_loss_total + (*PFSolStruct).Loss[i].imag;
		
		// Maximum loss amount and branch index
		P_loss_max = P_loss_max > (*PFSolStruct).Loss[i].real ? P_loss_max : (*PFSolStruct).Loss[i].real;
		Q_loss_max = Q_loss_max > (*PFSolStruct).Loss[i].imag ? Q_loss_max : (*PFSolStruct).Loss[i].imag;
		// Only save branch index not the value
		if (P_loss_max == (*PFSolStruct).Loss[i].real) (*PFSolStruct).P_loss_max_line_ind = i;
		if (Q_loss_max == (*PFSolStruct).Loss[i].imag) (*PFSolStruct).Q_loss_max_line_ind = i;
		
		// Branch charging injection
		(*PFSolStruct).Branch_Q_inj += (*PFSolStruct).FCHG[i] + (*PFSolStruct).TCHG[i];

		// Number of transformers
		if (TapRatio[i] > 0) (*PFSolStruct).Transformers_number++;
	}
	

	for (i = 0; i < Gen_Buses_Number_Of_Elements; i++) {
		// Total gen capacity
		(*PFSolStruct).Pg_total = (*PFSolStruct).Pg_total + Pmax[i];
		(*PFSolStruct).Qg_total_min = (*PFSolStruct).Qg_total_min + Qmin_[i];
		(*PFSolStruct).Qg_total_max = (*PFSolStruct).Qg_total_max + Qmax_[i];

		// Total on-line capacity
		(*PFSolStruct).Pg_total_online = (*PFSolStruct).Pg_total_online + (Pmax[i]*GenStatus[i]);
		(*PFSolStruct).Qg_total_min_online = (*PFSolStruct).Qg_total_min_online + (Qmin_[i] * GenStatus[i]);
		(*PFSolStruct).Qg_total_max_online = (*PFSolStruct).Qg_total_max_online + (Qmax_[i] * GenStatus[i]);

		//Actual generation
		(*PFSolStruct).Pg_actual = (*PFSolStruct).Pg_actual + (*PFSolStruct).Pg[i];
		(*PFSolStruct).Qg_actual = (*PFSolStruct).Qg_actual + (*PFSolStruct).Qg[i];
	}
	
	
	// if it's zero by default then we will not be able to find minimum Vm
	Vm_min = Vm_cal[0]; 
	Vang_min = (*PFSolStruct).Vang_deg[0];
	Vm_max = Vm_cal[0];
	Vang_max = (*PFSolStruct).Vang_deg[0];
	for (i = 0; i < Num_of_Buses; i++) {		
		// Total Load
		(*PFSolStruct).Pd_total = (*PFSolStruct).Pd_total + Pd[i];
		(*PFSolStruct).Qd_total = (*PFSolStruct).Qd_total + Qd[i];	

		// Min and max of voltage mag. and angle
		Vm_max = Vm_max > Vm_cal[i] ? Vm_max : Vm_cal[i];
		Vm_min = Vm_min < Vm_cal[i] ? Vm_min : Vm_cal[i];
		Vang_max = Vang_max > (*PFSolStruct).Vang_deg[i] ? Vang_max : (*PFSolStruct).Vang_deg[i];
		Vang_min = Vang_min < (*PFSolStruct).Vang_deg[i] ? Vang_min : (*PFSolStruct).Vang_deg[i];		
		// Save bus index
		if (Vm_max == Vm_cal[i]) (*PFSolStruct).Vmax_mag_bus = i;		
		if (Vm_min == Vm_cal[i]) (*PFSolStruct).Vmin_mag_bus = i;
		if (Vang_min == (*PFSolStruct).Vang_deg[i]) (*PFSolStruct).Vmin_ang_bus = i;
		if (Vang_max == (*PFSolStruct).Vang_deg[i]) (*PFSolStruct).Vmax_ang_bus = i;

		// Shunt injection
		(*PFSolStruct).Shunt_P_inj = (*PFSolStruct).Shunt_P_inj + (Vm_cal[i] * Vm_cal[i] * Gs[i]);
		(*PFSolStruct).Shunt_Q_inj = (*PFSolStruct).Shunt_Q_inj + (Vm_cal[i] * Vm_cal[i] * Bs[i]);	

		// Number of shunts
		if (Gs[i] != 0 || Bs[i] != 0) (*PFSolStruct).Shunt_number++;

		// Number of loads
		if (Pd[i]!=0 || Qd[i]!=0) (*PFSolStruct).Loads_number++;
	}
	
	(*PFSolStruct).Shunt_P_inj = (*PFSolStruct).Shunt_P_inj * -1;
	(*PFSolStruct).Gen_total_number = Gen_Buses_Number_Of_Elements;
	(*PFSolStruct).Gen_online_number = OnlineGenData_Num_Of_Elements;
	(*PFSolStruct).Gen_offline_number = Gen_Buses_Number_Of_Elements - OnlineGenData_Num_Of_Elements;
	(*PFSolStruct).Bus_number = Num_of_Buses;
	(*PFSolStruct).Branch_in_number = Num_Online_Branches;
	(*PFSolStruct).Branch_out_number = Num_Offline_Branches;

	timer_stop[6] = dsecnd(); // Misc. calculation time
	(*PFSolStruct).Timer[6] = (timer_stop[6] - timer_start[6]);

	/*printf("Online Gens = %llu  Online Gens - No PQ = %llu\n", OnlineGenData_Num_Of_Elements, Num_Online_Gens_No_PQ);
	printf("No. of in-service branches = %llu , No. of out-of-service branches = %llu\n\n", Num_Online_Branches, Num_Offline_Branches);
	*/



memory_free:
	
	if (Enable_Memory_Stat == 1) {
		(*PFSolStruct).Peak_Mem_Usage = mkl_peak_mem_usage(MKL_PEAK_MEM) / 1024;
		mkl_peak_mem_usage(MKL_PEAK_MEM_DISABLE);
	}

	mkl_free(V_rect);
	mkl_free(V_GenBus);
	mkl_free(Ybus);
	mkl_free(YbusV);
	mkl_free(Sbus);
	mkl_free(Sf);
	mkl_free(St);
	mkl_free(Ybus_val);
	mkl_free(Ybus_IA);
	mkl_free(Ybus_JA);
	mkl_free(Yf_New);	
	mkl_free(Yf_New_IA);
	mkl_free(Yf_New_JA);
	mkl_free(Yt_New);
	mkl_free(Yt_New_IA);
	mkl_free(Yt_New_JA);
	mkl_free(Vdrop);
	mkl_free(Vf);
	mkl_free(Vt);
	mkl_free(V_drop_diag);
	mkl_free(VddA);
	mkl_free(dYsc);
	mkl_free(V_mag);
	mkl_free(Vm_diag_recip);
	mkl_free(Bc_diag);
	mkl_free(Vrect_conj_diag);
	mkl_free(Tap_diag);
	mkl_free(B);
	mkl_free(B_conj_B);
	mkl_free(dlosstemp);
	mkl_free(BcTap);
	mkl_free(CfVm);	
	mkl_free(CfVm_diag);
	mkl_free(CtVm);
	mkl_free(BcTapCfVm);
	mkl_free(BcCtVm);
	mkl_free(CtVm_diag);
	
			
	mkl_sparse_destroy(Yf_CSR);
	mkl_sparse_destroy(Yt_CSR);
	mkl_sparse_destroy(Ybus_CSR);
	mkl_sparse_destroy(dYsc_CSR_Handle);
	mkl_sparse_destroy(Bc_CSR_Handle);
	mkl_sparse_destroy(Tap_CSR_Handle);
	mkl_sparse_destroy(Vdrop_CSR_Handle);
	mkl_sparse_destroy(Vm_reci_CSR_Handle);
	mkl_sparse_destroy(Vconj_CSR_Handle);
	mkl_sparse_destroy(VddA_CSR_Handle);
	mkl_sparse_destroy(B_CSR_Handle);
	mkl_sparse_destroy(B_m_conjB_CSR_Handle);
	mkl_sparse_destroy(B_p_conjB_CSR_Handle);
	mkl_sparse_destroy(dYscBp_CSR_Handle);
	mkl_sparse_destroy(BcTap_CSR_Handle);
	mkl_sparse_destroy(CfVm_CSR_Handle);
	mkl_sparse_destroy(CfVmCf_CSR_Handle);
	mkl_sparse_destroy(CtVm_CSR_Handle);
	mkl_sparse_destroy(CtVmCt_CSR_Handle);
	
	free(GenBus);
	free(Pd_gbus);
	free(Qd_gbus);	
	free(GenCount);
	free(Qg_total);
	free(Qg_min);
	free(Qg_max);
	free(Qmin);
	free(Qmax);	

	if (status != 0) {
		
		free((*PFSolStruct).Vm);
		free((*PFSolStruct).Vang_deg);
		free((*PFSolStruct).Pg);
		free((*PFSolStruct).Qg);
		free((*PFSolStruct).PF);
		free((*PFSolStruct).QF);
		free((*PFSolStruct).PT);
		free((*PFSolStruct).QT);
		free((*PFSolStruct).FCHG);
		free((*PFSolStruct).TCHG);

		mkl_free((*PFSolStruct).Loss);
		mkl_free((*PFSolStruct).dloss_dV_a_Dense);
		mkl_free((*PFSolStruct).dloss_dV_a_IA);
		mkl_free((*PFSolStruct).dloss_dV_a_JA);
		mkl_free((*PFSolStruct).dloss_dV_a_pointerE);
		mkl_free((*PFSolStruct).dloss_dV_a_val);
		mkl_free((*PFSolStruct).dloss_dV_m_Dense);
		mkl_free((*PFSolStruct).dloss_dV_m_IA);
		mkl_free((*PFSolStruct).dloss_dV_m_JA);
		mkl_free((*PFSolStruct).dloss_dV_m_pointerE);
		mkl_free((*PFSolStruct).dloss_dV_m_val);
		mkl_free((*PFSolStruct).dchg_dVm_f_Dense);
		mkl_free((*PFSolStruct).dchg_dVm_f_IA);
		mkl_free((*PFSolStruct).dchg_dVm_f_JA);
		mkl_free((*PFSolStruct).dchg_dVm_f_pointerE);
		mkl_free((*PFSolStruct).dchg_dVm_f_val);
		mkl_free((*PFSolStruct).dchg_dVm_t_Dense);
		mkl_free((*PFSolStruct).dchg_dVm_t_IA);
		mkl_free((*PFSolStruct).dchg_dVm_t_JA);
		mkl_free((*PFSolStruct).dchg_dVm_t_pointerE);
		mkl_free((*PFSolStruct).dchg_dVm_t_val);

		mkl_sparse_destroy((*PFSolStruct).dloss_dV_a_CSR_Handle);
		mkl_sparse_destroy((*PFSolStruct).dloss_dV_m_CSR_Handle);
		mkl_sparse_destroy((*PFSolStruct).dchg_dVm_f_CSR_Handle);
		mkl_sparse_destroy((*PFSolStruct).dchg_dVm_t_CSR_Handle);
					   		 		
		mkl_free_buffers();
	}

	(*PFSolStruct).status = status; 

	timer_stop[0] = dsecnd(); //Total execution time
	(*PFSolStruct).Timer[0] = (timer_stop[0] - timer_start[0]); //Total execution time

	
} //end of PF_NR_MKL function