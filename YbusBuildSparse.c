#include "YbusBuildSparse.h"
#define ALIGN 64

/* To avoid constantly repeating the part of code that checks different functions' status, using the below macros */
#define CHECK_SPARSE(function)  do {  if(function != SPARSE_STATUS_SUCCESS){ status = 3; goto memory_free; } } while(0) // Check sparse function status
#define CHECK_MEM(variable)  do {  if(variable == NULL){ status = 5; goto memory_free; } } while(0) // check memory allocation

//// This finction converts a dense matrix to CSR format
//int dense2CSR(MKL_Complex16 *A_, MKL_INT M, MKL_INT N, MKL_INT nzmax, MKL_Complex16 *A_Vals, MKL_INT *A_Col, MKL_INT *A_Row)
//{ 
//// A : Original dense matrix
//// M : Number of rows of the matrix A.
//// N : Number of columns of the matrix A.
//// lda: Leading dimension of A
//// A_Vals : Array containing non-zero values of A
//// A_Row: Contains the row indices plus one for each non - zero element of the matrix A.
//// A_Col: Contains the column indices plus one for each non - zero element of the matrix A.
//	MKL_INT info = 0; // If info=0, the execution is successful.
//	MKL_INT job[6];
//	job[0] = 0; // the rectangular matrix A is converted to the CSR format;
//	job[1] = 0; // zero-based indexing for the rectangular matrix A is used;
//	job[2] = 0; // zero-based indexing for the matrix in CSR format is used;
//	job[3] = 2; // whole matrix
//	job[4] = nzmax; // maximum number of the non-zero elements allowed if job[0] = 0
//	job[5] = 5; // If job[5]>0, arrays acsr, ia, ja are generated for the output storage. If job[5]=0, only array ia is generated for the output storage.
//	mkl_zdnscsr(job, &M, &N, A_, &N, A_Vals, A_Col, A_Row, &info);
//	if (info == 0) { return 0; } // if info = 0, operation of mkl_zdnscsr was successful
//	else { return 1; };
//}
//
//
//// This function converts a matrix in CSR format to dense matrix
//int CSR2dense(MKL_Complex16 *A, MKL_INT M, MKL_INT N, MKL_Complex16 *A_Vals, MKL_INT *A_Row, MKL_INT *A_Col)
//{
//// A : Original dense matrix
//// M : Number of rows of the matrix A.
//// N : Number of columns of the matrix A.
//// lda: Leading dimension of A
//// A_Vals : Array containing non-zero values of A
//// A_Row: Contains the row indices plus one for each non - zero element of the matrix A.
//// A_Col: Contains the column indices plus one for each non - zero element of the matrix A.
//
//	MKL_INT info = 0; // If info=0, the execution is successful.
//	MKL_INT job[6];
//	job[0] = 1; // the rectangular matrix A is restored from the CSR format;
//	job[1] = 0; // zero-based indexing for the rectangular matrix A is used;
//	job[2] = 0; // zero-based indexing for the matrix in CSR format is used;
//	job[3] = 2; // whole matrix
//	//job[4] = nzmax; // maximum number of the non-zero elements allowed if job[0] = 0
//	job[5] = 1; // If job[5]>0, arrays acsr, ia, ja are generated for the output storage. If job[5]=0, only array ia is generated for the output storage.
//	mkl_zdnscsr(job, &M, &N, A, &N, A_Vals, A_Col, A_Row, &info);
//	if (info == 0) { return 0; } // if info = 0, operation of mkl_zdnscsr was successful
//	else { return 1; };
//}
//


// ******************************************************************************** Main Routine ***************************************************************
int MakeYbusSparse(struct Ybus_Struct *YbusOut, double BaseMVA, unsigned long long int Num_of_Buses, unsigned long long int Num_of_Branches, int *FromBus, int *ToBus, int *BranchStatus, double *r, double *x, double *b, double *TapRatio, double *ShiftAngle, double *Gs, double *Bs, double **timer)
{
	
double start, stop;
start = dsecnd();

// ******************************************************** Variables ********************************************************

// Enable memory usage statistics
int N_AllocatedBuffers = 0;
if ( YbusOut->Enable_Memory_Stat == 1 ) mkl_peak_mem_usage(MKL_PEAK_MEM_ENABLE);

//int Max_Num_Threads_Avl = omp_get_max_threads(); // Find the maximum number of available threads
int OMP_Num_Threads_To_Use = YbusOut->Num_Threads; // Number of threads to be used in this function
int OMP_Enable = YbusOut->OMP_Enable; // 0: Disable OpenMP 1: Enable OpenMP

// **** Initial Variables ****
unsigned long long int i = 0, j = 0, k = 0;
unsigned long long int Fbus = 0, Tbus = 0;
double temp1 = 0, temp2 = 0;
int status = 0;
int Sym_Mat = 1; // 0: non-symmetric 1: symmetric : Symmetric Ybus or not? If ShiftAngle is non-zero then Ybus is not symmetric
unsigned long int Chunk_Size;
double *Ys_Re = NULL;
double *Ys_Img = NULL;
double *Tap_Re = NULL;
double *Tap_Img = NULL;
double *Ytt_Re = NULL;
double *Ytt_Img = NULL;
double *Yff_Re = NULL;
double *Yff_Img = NULL;
double *Yft_Re = NULL;
double *Yft_Img = NULL;
double *Ytf_Re = NULL;
double *Ytf_Img = NULL;

// scalars to be used in matrix-matrix product
MKL_Complex16 alpha, beta;
alpha.real = 1;
alpha.imag = 0;
beta.real = 0;
beta.imag = 0;


// **** Variables to be used inside MKL Dense ***
MKL_Complex16 *Ysh_Dense = NULL;
MKL_Complex16 *Ysc = NULL;


// **** Variables to be used inside MKL Sparse ***
// number of non-zero elements in each array
unsigned long long int Yf_nonzero = 0, Yt_nonzero = 0;
unsigned long long int Cf_nonzero = 0, Ct_nonzero = 0;


// Sparse Handles
sparse_matrix_t Yf_CSR = NULL, CfPrimeYf_CSR = NULL;
sparse_matrix_t Yt_CSR = NULL, CtPrimeYt_CSR = NULL;
sparse_matrix_t TapCf_CSR = NULL;

// CSR export variables
MKL_INT rows, cols, rows1, cols1;
sparse_index_base_t indexing = 0, indexing1 = 0;



// ****************************************************** Initialize Data ********************************************************
/*
%% for each branch, compute the elements of the branch admittance matrix where
%%
%%      | If |   | Yff  Yft |   | Vf |
%%      |    | = |          | * |    |
%%      | It |   | Ytf  Ytt |   | Vt |
%%
*/

// Allocate memeory for inital variables
(*YbusOut).Ysc = (MKL_Complex16 *)mkl_calloc((size_t)(Num_of_Branches), sizeof(MKL_Complex16), ALIGN);
CHECK_MEM((*YbusOut).Ysc);
Ys_Re = calloc((size_t)Num_of_Branches, sizeof(double));
CHECK_MEM(Ys_Re);
Ys_Img = calloc((size_t)Num_of_Branches, sizeof(double));
CHECK_MEM(Ys_Img);
Tap_Re = calloc((size_t)Num_of_Branches, sizeof(double));
CHECK_MEM(Tap_Re);
Tap_Img = calloc((size_t)Num_of_Branches, sizeof(double));
CHECK_MEM(Tap_Img);
Ytt_Re = calloc((size_t)Num_of_Branches, sizeof(double));
CHECK_MEM(Ytt_Re);
Ytt_Img = calloc((size_t)Num_of_Branches, sizeof(double));
CHECK_MEM(Ytt_Img);
Yff_Re = calloc((size_t)Num_of_Branches, sizeof(double));
CHECK_MEM(Yff_Re);
Yff_Img = calloc((size_t)Num_of_Branches, sizeof(double));
CHECK_MEM(Yff_Img);
Yft_Re = calloc((size_t)Num_of_Branches, sizeof(double));
CHECK_MEM(Yft_Re);
Yft_Img = calloc((size_t)Num_of_Branches, sizeof(double));
CHECK_MEM(Yft_Img);
Ytf_Re = calloc((size_t)Num_of_Branches, sizeof(double));
CHECK_MEM(Ytf_Re);
Ytf_Img = calloc((size_t)Num_of_Branches, sizeof(double));
CHECK_MEM(Ytf_Img);


Chunk_Size = round(Num_of_Branches / OMP_Num_Threads_To_Use);
#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Num_Threads_To_Use) private(i,temp1,temp2) shared(Chunk_Size,Num_of_Branches,BranchStatus,Sym_Mat,r,x,b,Ys_Re,Ys_Img,Tap_Re,Tap_Img,TapRatio,ShiftAngle,Ytt_Re,Ytt_Img,YbusOut) if(OMP_Enable == 1)
for (i = 0; i < Num_of_Branches; ++i) {

	if (BranchStatus[i] == 1) { // Series admittance: BranchStatus[i]/(r[i]+x[i]*I)
		temp1 = 1 / rec2pol_mag(r[i], x[i]);
		temp2 = (-1)*rec2pol_ang(r[i], x[i]);
		Ys_Re[i] = temp1 * cos(temp2);
		(*YbusOut).Ysc[i].real = Ys_Re[i];
		Ys_Img[i] = temp1 * sin(temp2);
		(*YbusOut).Ysc[i].imag = -Ys_Img[i];
	}


	if (TapRatio[i] != 0) {
		Tap_Re[i] = TapRatio[i] * cos(deg2rad(ShiftAngle[i])); // Add phase shifter to tap ratio: TapRatio[i] * cexpl(deg2rad(ShiftAngle[i])*I)
		Tap_Img[i] = TapRatio[i] * sin(deg2rad(ShiftAngle[i]));
	}
	else {
		Tap_Re[i] = 1 * cos(deg2rad(ShiftAngle[i])); // Add phase shifter to tap ratio: TapRatio[i] * cexpl(deg2rad(ShiftAngle[i])*I)
		Tap_Img[i] = 1 * sin(deg2rad(ShiftAngle[i]));
	}

	Ytt_Re[i] = Ys_Re[i]; //Ys[i] + (Bc[i]/2)*j
	Ytt_Img[i] = Ys_Img[i] + (b[i] / 2) * BranchStatus[i];

	if (ShiftAngle[i] != 0) Sym_Mat = 0;
}


Chunk_Size = round(Num_of_Branches / OMP_Num_Threads_To_Use);
#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Num_Threads_To_Use) private(i,temp1,temp2,Fbus,Tbus) shared(Chunk_Size,Num_of_Buses,Num_of_Branches,Ys_Re,Ys_Img,Tap_Re,Tap_Img,Ytt_Re,Ytt_Img, Yff_Re, Yff_Img, Yft_Re, Yft_Img, Ytf_Re, Ytf_Img,FromBus,ToBus) if(OMP_Enable == 1)
for (i = 0; i < Num_of_Branches; ++i) {

	Yff_Re[i] = Ytt_Re[i] / ((Tap_Re[i] * Tap_Re[i]) + (Tap_Img[i] * Tap_Img[i])); // Ytt[i] / (Tap[i] * conjl(Tap[i]))
	Yff_Img[i] = Ytt_Img[i] / ((Tap_Re[i] * Tap_Re[i]) + (Tap_Img[i] * Tap_Img[i]));

	temp1 = rec2pol_mag((-1)*Ys_Re[i], (-1)*Ys_Img[i]) / rec2pol_mag(Tap_Re[i], (-1)*Tap_Img[i]);
	temp2 = rec2pol_ang((-1)*Ys_Re[i], (-1)*Ys_Img[i]) - rec2pol_ang(Tap_Re[i], (-1)*Tap_Img[i]);
	Yft_Re[i] = temp1 * cos(temp2); //(Ys[i]*(-1)) / conjl(Tap[i])
	Yft_Img[i] = temp1 * sin(temp2);

	temp1 = rec2pol_mag((-1)*Ys_Re[i], (-1)*Ys_Img[i]) / rec2pol_mag(Tap_Re[i], Tap_Img[i]);
	temp2 = rec2pol_ang((-1)*Ys_Re[i], (-1)*Ys_Img[i]) - rec2pol_ang(Tap_Re[i], Tap_Img[i]);
	Ytf_Re[i] = temp1 * cos(temp2); //(Ys[i]*(-1)) / conjl(Tap[i])
	Ytf_Img[i] = temp1 * sin(temp2);

}





// ****************************************************** Build Ybus ********************************************************

// Calculate Ybus = Cf' * Yf + Ct' * Yt  +.... ::Branch admittances
//                  + Ysh :: Shunt admittance

switch ((*YbusOut).CoreEngine)
{

case 0: {

	/* ******************************** MKL Dense ******************************** */

	// **** Variables to be used in Initialize Data routine ****	
	(*YbusOut).Yf_Dense = (MKL_Complex16 *)mkl_calloc((size_t)(Num_of_Buses*Num_of_Branches), sizeof(MKL_Complex16), ALIGN); // Yf = [ Yff Yft]
	CHECK_MEM((*YbusOut).Yf_Dense);
	(*YbusOut).Yt_Dense = (MKL_Complex16 *)mkl_calloc((size_t)(Num_of_Buses*Num_of_Branches), sizeof(MKL_Complex16), ALIGN); // Yt = [Ytf Ytt]
	CHECK_MEM((*YbusOut).Yt_Dense);
	
	// Size of Cf is Num_of_Branches*Num_of_Buses, so the size of CfPrime (transposed) is Num_of_Buses * Num_of_Branches
	(*YbusOut).CfPrime_Dense = (MKL_Complex16 *)mkl_calloc((size_t)(Num_of_Buses * Num_of_Branches), sizeof(MKL_Complex16), ALIGN); // connection matrix for line & from buses (CfPrime actually contains transposed connection matrix to improve computation performance)
	CHECK_MEM((*YbusOut).CfPrime_Dense);
	(*YbusOut).CtPrime_Dense = (MKL_Complex16 *)mkl_calloc((size_t)(Num_of_Buses * Num_of_Branches), sizeof(MKL_Complex16), ALIGN); // connection matrix for line & to buses (CtPrime actually contains transposed connection matrix to improve computation performance)
	CHECK_MEM((*YbusOut).CtPrime_Dense);

	// Size of Tap is Num_of_Branches * Num_of_Branches
	(*YbusOut).Tap_Dense = (MKL_Complex16 *)mkl_calloc((size_t)(Num_of_Branches*Num_of_Branches), sizeof(MKL_Complex16), ALIGN);
	CHECK_MEM((*YbusOut).Tap_Dense);

	// Form the A matrix for use in PFsolution (to calculate line losses) : A = Tap * Cf - Ct : Num_of_Branches * Num_of_Buses	
	(*YbusOut).A_Dense = (MKL_Complex16 *)mkl_calloc((size_t)(Num_of_Branches*Num_of_Buses), sizeof(MKL_Complex16), ALIGN);
	CHECK_MEM((*YbusOut).A_Dense);

	Chunk_Size = round(Num_of_Branches / OMP_Num_Threads_To_Use);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Num_Threads_To_Use) private(i,Fbus,Tbus) shared(Chunk_Size,Num_of_Buses,Num_of_Branches,Ytt_Re,Ytt_Img,Yff_Re,Yff_Img,Yft_Re,Yft_Img,Ytf_Re,Ytf_Img,FromBus,ToBus,YbusOut,Tap_Re,Tap_Img,BranchStatus) if(OMP_Enable == 1)
	for (i = 0; i < Num_of_Branches; ++i) {
		Fbus = FromBus[i] - 1;
		Tbus = ToBus[i] - 1;

		// Form CfPrime, CtPrime, Yf, and Yt
		(*YbusOut).CfPrime_Dense[Fbus*Num_of_Branches + i].real = 1;
		(*YbusOut).CfPrime_Dense[Fbus*Num_of_Branches + i].imag = 0;

		(*YbusOut).CtPrime_Dense[Tbus*Num_of_Branches + i].real = 1;
		(*YbusOut).CtPrime_Dense[Tbus*Num_of_Branches + i].imag = 0;

		// A = Ct
		(*YbusOut).A_Dense[i*Num_of_Buses + Tbus].real = 1;
		(*YbusOut).A_Dense[i*Num_of_Buses + Tbus].imag = 0;
		
		(*YbusOut).Yf_Dense[i*Num_of_Buses + Fbus].real = Yff_Re[i];
		(*YbusOut).Yf_Dense[i*Num_of_Buses + Fbus].imag = Yff_Img[i];
		(*YbusOut).Yf_Dense[i*Num_of_Buses + Tbus].real = Yft_Re[i];
		(*YbusOut).Yf_Dense[i*Num_of_Buses + Tbus].imag = Yft_Img[i];

		(*YbusOut).Yt_Dense[i*Num_of_Buses + Fbus].real = Ytf_Re[i];
		(*YbusOut).Yt_Dense[i*Num_of_Buses + Fbus].imag = Ytf_Img[i];
		(*YbusOut).Yt_Dense[i*Num_of_Buses + Tbus].real = Ytt_Re[i];
		(*YbusOut).Yt_Dense[i*Num_of_Buses + Tbus].imag = Ytt_Img[i];

		// Form the diagonal matrix Tap = 1/(Tap_Re + jTap_Img) = (a-jb)/(a^2+b^2)
		if (BranchStatus[i] == 1) {
			(*YbusOut).Tap_Dense[i*Num_of_Branches + i].real = Tap_Re[i] / (Tap_Re[i] * Tap_Re[i] + Tap_Img[i] * Tap_Img[i]);
			(*YbusOut).Tap_Dense[i*Num_of_Branches + i].imag = -Tap_Img[i] / (Tap_Re[i] * Tap_Re[i] + Tap_Img[i] * Tap_Img[i]);
		}
		
	}
	
	// Use MKL routines for dense matrix multiplication and addition		 
	(*YbusOut).Ybus_Dense = (MKL_Complex16 *)mkl_calloc((size_t)(Num_of_Buses*Num_of_Buses), sizeof(MKL_Complex16), ALIGN);
	CHECK_MEM((*YbusOut).Ybus_Dense);
	
	// Ybus = Cf' * Yf + Ct' * Yt	
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Num_of_Buses, Num_of_Buses, Num_of_Branches, &alpha, (*YbusOut).CfPrime_Dense, Num_of_Branches, (*YbusOut).Yf_Dense, Num_of_Buses, &beta, (*YbusOut).Ybus_Dense, Num_of_Buses); // Ybus_Row_Format_Dense = Cf' * Yf
	beta.real = 1;
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Num_of_Buses, Num_of_Buses, Num_of_Branches, &alpha, (*YbusOut).CtPrime_Dense, Num_of_Branches, (*YbusOut).Yt_Dense, Num_of_Buses, &beta, (*YbusOut).Ybus_Dense, Num_of_Buses);	// Ybus_Row_Format_Dense = Ybus_Row_Format_Dense + (Ct' * Yt)
	
	// Add Ysh to diagonal elements (Shunt admittance)
	Chunk_Size = round(Num_of_Buses / OMP_Num_Threads_To_Use);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Num_Threads_To_Use) private(i) shared(Chunk_Size,Num_of_Buses,Gs,Bs,BaseMVA,YbusOut) if(OMP_Enable == 1)
	for (i = 0; i < Num_of_Buses; ++i)
	{
		(*YbusOut).Ybus_Dense[i*Num_of_Buses + i].real = (*YbusOut).Ybus_Dense[i*Num_of_Buses + i].real + Gs[i] / BaseMVA;
		(*YbusOut).Ybus_Dense[i*Num_of_Buses + i].imag = (*YbusOut).Ybus_Dense[i*Num_of_Buses + i].imag + Bs[i] / BaseMVA;
	}
	
	// A = Tap * Cf - A
	beta.real = -1;
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasTrans, Num_of_Branches, Num_of_Buses, Num_of_Branches, &alpha, (*YbusOut).Tap_Dense, Num_of_Branches, (*YbusOut).CfPrime_Dense, Num_of_Branches, &beta, (*YbusOut).A_Dense, Num_of_Buses);

	break; } // End of Case 0



case 1: {
	/* ********************************* MKL Sparse CSR *****************************/
	   	 
	Yf_nonzero = Num_of_Branches * 2;
	Yt_nonzero = Num_of_Branches * 2;
	Cf_nonzero = Num_of_Branches;
	Ct_nonzero = Num_of_Branches;
	
	(*YbusOut).Yf_Val = (MKL_Complex16 *)mkl_calloc((size_t)Yf_nonzero, sizeof(MKL_Complex16), ALIGN);
	CHECK_MEM((*YbusOut).Yf_Val);
	(*YbusOut).Yf_IA = (MKL_INT *)mkl_calloc((size_t)Yf_nonzero, sizeof(MKL_INT), ALIGN);
	CHECK_MEM((*YbusOut).Yf_IA);
	(*YbusOut).Yf_JA = (MKL_INT *)mkl_calloc((size_t)(Num_of_Branches + 1), sizeof(MKL_INT), ALIGN);
	CHECK_MEM((*YbusOut).Yf_JA);

	(*YbusOut).Yt_Val = (MKL_Complex16 *)mkl_calloc((size_t)Yt_nonzero, sizeof(MKL_Complex16), ALIGN);
	CHECK_MEM((*YbusOut).Yt_Val);
	(*YbusOut).Yt_IA = (MKL_INT *)mkl_calloc((size_t)Yt_nonzero, sizeof(MKL_INT), ALIGN);
	CHECK_MEM((*YbusOut).Yt_IA);
	(*YbusOut).Yt_JA = (MKL_INT *)mkl_calloc((size_t)(Num_of_Branches + 1), sizeof(MKL_INT), ALIGN);
	CHECK_MEM((*YbusOut).Yt_JA);

	(*YbusOut).Cf_Val = (MKL_Complex16 *)mkl_calloc((size_t)Cf_nonzero, sizeof(MKL_Complex16), ALIGN);
	CHECK_MEM((*YbusOut).Cf_Val);
	(*YbusOut).Cf_IA = (MKL_INT *)mkl_calloc((size_t)Cf_nonzero, sizeof(MKL_INT), ALIGN);
	CHECK_MEM((*YbusOut).Cf_IA);
	(*YbusOut).Cf_JA = (MKL_INT *)mkl_calloc((size_t)(Num_of_Branches + 1), sizeof(MKL_INT), ALIGN);
	CHECK_MEM((*YbusOut).Cf_JA);

	(*YbusOut).Ct_Val = (MKL_Complex16 *)mkl_calloc((size_t)Cf_nonzero, sizeof(MKL_Complex16), ALIGN);
	CHECK_MEM((*YbusOut).Ct_Val);
	(*YbusOut).Ct_IA = (MKL_INT *)mkl_calloc((size_t)Cf_nonzero, sizeof(MKL_INT), ALIGN);
	CHECK_MEM((*YbusOut).Ct_IA);
	(*YbusOut).Ct_JA = (MKL_INT *)mkl_calloc((size_t)(Num_of_Branches + 1), sizeof(MKL_INT), ALIGN);
	CHECK_MEM((*YbusOut).Ct_JA);

	// Size of Tap is Num_of_Branches * Num_of_Branches
	(*YbusOut).Tap_Val = (MKL_Complex16 *)mkl_calloc((size_t)(Num_of_Branches), sizeof(MKL_Complex16), ALIGN);
	CHECK_MEM((*YbusOut).Tap_Val);
	(*YbusOut).Tap_IA = (MKL_INT *)mkl_calloc((size_t)(Num_of_Branches), sizeof(MKL_INT), ALIGN);
	CHECK_MEM((*YbusOut).Tap_IA);
	(*YbusOut).Tap_JA = (MKL_INT *)mkl_calloc((size_t)(Num_of_Branches + 1), sizeof(MKL_INT), ALIGN);
	CHECK_MEM((*YbusOut).Tap_JA);
	
	(*YbusOut).A_Val = (MKL_Complex16 *)mkl_calloc((size_t)Cf_nonzero, sizeof(MKL_Complex16), ALIGN);
	(*YbusOut).A_IA = (MKL_INT *)mkl_calloc((size_t)Cf_nonzero, sizeof(MKL_INT), ALIGN);
	(*YbusOut).A_JA = (MKL_INT *)mkl_calloc((size_t)(Num_of_Branches + 1), sizeof(MKL_INT), ALIGN);

	// General IA and JA for diagonal matrix Tap
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Num_Threads_To_Use) private(i) shared(Chunk_Size,Num_of_Branches,YbusOut) if(OMP_Enable == 1)
	for (i = 0; i < (Num_of_Branches); i++) {
		(*YbusOut).Tap_IA[i] = i;
		(*YbusOut).Tap_JA[i] = i;
	}
	(*YbusOut).Tap_JA[Num_of_Branches] = Num_of_Branches;

	
	Chunk_Size = round(Num_of_Branches / OMP_Num_Threads_To_Use);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Num_Threads_To_Use) private(i, Fbus, Tbus) shared(Chunk_Size,Num_of_Branches,FromBus,ToBus,Yff_Re,Yff_Img,Yft_Re,Yft_Img,Ytf_Re,Ytf_Img,Ytt_Re,Ytt_Img,YbusOut,BranchStatus,Tap_Re, Tap_Img) if(OMP_Enable == 1)
	for (i = 0; i < Num_of_Branches; i++) {

		Fbus = FromBus[i] - 1; // From Bus
		Tbus = ToBus[i] - 1; // To Bus

		// Cf		
		(*YbusOut).Cf_JA[i] = i;
		(*YbusOut).Cf_IA[i] = Fbus;
		(*YbusOut).Cf_Val[i].real = 1;
		(*YbusOut).Cf_Val[i].imag = 0;
		// Ct
		(*YbusOut).Ct_JA[i] = i;
		(*YbusOut).Ct_IA[i] = Tbus;
		(*YbusOut).Ct_Val[i].real = 1;
		(*YbusOut).Ct_Val[i].imag = 0;
		
		// Yf 
		(*YbusOut).Yf_JA[i] = i * 2;
		(*YbusOut).Yf_Val[i * 2].real = Yff_Re[i];
		(*YbusOut).Yf_Val[i * 2].imag = Yff_Img[i];
		(*YbusOut).Yf_Val[(i * 2) + 1].real = Yft_Re[i];
		(*YbusOut).Yf_Val[(i * 2) + 1].imag = Yft_Img[i];
		(*YbusOut).Yf_IA[i * 2] = Fbus;
		(*YbusOut).Yf_IA[(i * 2) + 1] = Tbus;

		// Yt
		(*YbusOut).Yt_JA[i] = i * 2;
		(*YbusOut).Yt_Val[i * 2].real = Ytf_Re[i];
		(*YbusOut).Yt_Val[i * 2].imag = Ytf_Img[i];
		(*YbusOut).Yt_Val[(i * 2) + 1].real = Ytt_Re[i];
		(*YbusOut).Yt_Val[(i * 2) + 1].imag = Ytt_Img[i];
		(*YbusOut).Yt_IA[i * 2] = Fbus;
		(*YbusOut).Yt_IA[(i * 2) + 1] = Tbus;

		// Form the diagonal matrix Tap = 1/(Tap_Re + jTap_Img) = (a-jb)/(a^2+b^2)
		if (BranchStatus[i] == 1) {
			(*YbusOut).Tap_Val[i].real = Tap_Re[i] / (Tap_Re[i] * Tap_Re[i] + Tap_Img[i] * Tap_Img[i]);
			(*YbusOut).Tap_Val[i].imag = -Tap_Img[i] / (Tap_Re[i] * Tap_Re[i] + Tap_Img[i] * Tap_Img[i]);
		}

	}
	(*YbusOut).Cf_JA[Num_of_Branches] = Cf_nonzero;
	(*YbusOut).Ct_JA[Num_of_Branches] = Ct_nonzero;
	(*YbusOut).Yf_JA[Num_of_Branches] = Yf_nonzero;
	(*YbusOut).Yt_JA[Num_of_Branches] = Yt_nonzero;	

	CHECK_SPARSE(mkl_sparse_z_create_csr(&(*YbusOut).Cf_CSR_Handle, SPARSE_INDEX_BASE_ZERO, Num_of_Branches, Num_of_Buses, (*YbusOut).Cf_JA, (*YbusOut).Cf_JA + 1, (*YbusOut).Cf_IA, (*YbusOut).Cf_Val));
	CHECK_SPARSE(mkl_sparse_z_create_csr(&Yf_CSR, SPARSE_INDEX_BASE_ZERO, Num_of_Branches, Num_of_Buses, (*YbusOut).Yf_JA, (*YbusOut).Yf_JA + 1, (*YbusOut).Yf_IA, (*YbusOut).Yf_Val));
	CHECK_SPARSE(mkl_sparse_spmm(SPARSE_OPERATION_TRANSPOSE, (*YbusOut).Cf_CSR_Handle, Yf_CSR, &CfPrimeYf_CSR));
	CHECK_SPARSE(mkl_sparse_order(CfPrimeYf_CSR));

	CHECK_SPARSE(mkl_sparse_z_create_csr(&(*YbusOut).Ct_CSR_Handle, SPARSE_INDEX_BASE_ZERO, Num_of_Branches, Num_of_Buses, (*YbusOut).Ct_JA, (*YbusOut).Ct_JA + 1, (*YbusOut).Ct_IA, (*YbusOut).Ct_Val));
	CHECK_SPARSE(mkl_sparse_z_create_csr(&Yt_CSR, SPARSE_INDEX_BASE_ZERO, Num_of_Branches, Num_of_Buses, (*YbusOut).Yt_JA, (*YbusOut).Yt_JA + 1, (*YbusOut).Yt_IA, (*YbusOut).Yt_Val));
	CHECK_SPARSE(mkl_sparse_spmm(SPARSE_OPERATION_TRANSPOSE, (*YbusOut).Ct_CSR_Handle, Yt_CSR, &CtPrimeYt_CSR));
	CHECK_SPARSE(mkl_sparse_order(CtPrimeYt_CSR));
	
	alpha.real = 1;
	alpha.imag = 0;
	CHECK_SPARSE(mkl_sparse_z_add(SPARSE_OPERATION_NON_TRANSPOSE, CfPrimeYf_CSR, alpha, CtPrimeYt_CSR, &(*YbusOut).Ybus_CSR_Handle));
	CHECK_SPARSE(mkl_sparse_z_export_csr((*YbusOut).Ybus_CSR_Handle, &indexing, &rows, &cols, &(*YbusOut).Ybus_JA, &(*YbusOut).pointerE_Ybus, &(*YbusOut).Ybus_IA, &(*YbusOut).Ybus_Values));

	// A = Tap * Cf - Ct
	CHECK_SPARSE(mkl_sparse_z_create_csr(&(*YbusOut).Tap_CSR_Handle, SPARSE_INDEX_BASE_ZERO, Num_of_Branches, Num_of_Branches, (*YbusOut).Tap_JA, (*YbusOut).Tap_JA + 1, (*YbusOut).Tap_IA, (*YbusOut).Tap_Val));
	CHECK_SPARSE(mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, (*YbusOut).Tap_CSR_Handle, (*YbusOut).Cf_CSR_Handle, &TapCf_CSR));
	alpha.real = -1;
	CHECK_SPARSE(mkl_sparse_z_add(SPARSE_OPERATION_NON_TRANSPOSE, (*YbusOut).Ct_CSR_Handle, alpha, TapCf_CSR, &(*YbusOut).A_CSR_Handle));
	CHECK_SPARSE(mkl_sparse_z_export_csr((*YbusOut).A_CSR_Handle, &indexing1, &rows1, &cols1, &(*YbusOut).A_JA, &(*YbusOut).pointerE_A, &(*YbusOut).A_IA, &(*YbusOut).A_Val));
	   

	// Add Ysh to diagonal elements of Ybus
	Chunk_Size = round(Num_of_Buses / OMP_Num_Threads_To_Use);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Num_Threads_To_Use) private(i,j) shared(Chunk_Size,Num_of_Buses,Gs, Bs, YbusOut, BaseMVA) if(OMP_Enable == 1)
	for (i = 0; i < Num_of_Buses; i++) {
		if (Gs[i] != 0 || Bs[i] != 0) {
			for (j = (*YbusOut).Ybus_JA[i]; j < (*YbusOut).Ybus_JA[i + 1]; j++) {
				if ((*YbusOut).Ybus_IA[j] == i) {
					YbusOut->Ybus_Values[j].real = YbusOut->Ybus_Values[j].real + (Gs[i] / BaseMVA);
					YbusOut->Ybus_Values[j].imag = YbusOut->Ybus_Values[j].imag + (Bs[i] / BaseMVA);
				}
			}
		}
	}

	
	break; } // End of Case 1


}



memory_free:

	if (YbusOut->Enable_Memory_Stat == 1) {
		YbusOut->Peak_Mem_Used = mkl_peak_mem_usage(MKL_PEAK_MEM) / 1024; // Return peak memory used in KB					
	}
	
	
	free(Ys_Re);
	free(Ys_Img);
	free(Tap_Re);
	free(Tap_Img);
	free(Ytt_Re);
	free(Ytt_Img);
	free(Ytf_Re);
	free(Ytf_Img);
	free(Yff_Re);
	free(Yff_Img);
	free(Yft_Re);
	free(Yft_Img);

	mkl_sparse_destroy(Yf_CSR);
	mkl_sparse_destroy(Yt_CSR);

	mkl_sparse_destroy(CfPrimeYf_CSR);
	mkl_sparse_destroy(CtPrimeYt_CSR);	
	mkl_sparse_destroy(TapCf_CSR);


if (status != 0) { // unsuccessful operation	

	mkl_free((*YbusOut).Yt_Dense);
	mkl_free((*YbusOut).Yf_Dense);
	mkl_free((*YbusOut).Ysc);
	mkl_free((*YbusOut).Yf_Val);
	mkl_free((*YbusOut).Yt_Val);
	mkl_free((*YbusOut).Yf_IA);
	mkl_free((*YbusOut).Yf_JA);
	mkl_free((*YbusOut).Yt_IA);
	mkl_free((*YbusOut).Yt_JA);	
	mkl_free((*YbusOut).Ybus_Dense);
	mkl_free((*YbusOut).A_Dense);
	mkl_free((*YbusOut).CfPrime_Dense);
	mkl_free((*YbusOut).CtPrime_Dense);
	mkl_free((*YbusOut).Cf_Val);
	mkl_free((*YbusOut).Cf_IA);
	mkl_free((*YbusOut).Cf_JA);
	mkl_free((*YbusOut).Ct_Val);
	mkl_free((*YbusOut).Ct_IA);
	mkl_free((*YbusOut).Ct_JA);
	mkl_sparse_destroy((*YbusOut).Ybus_CSR_Handle);	
	mkl_sparse_destroy((*YbusOut).Cf_CSR_Handle);
	mkl_sparse_destroy((*YbusOut).Ct_CSR_Handle);
	mkl_sparse_destroy((*YbusOut).A_CSR_Handle);
	mkl_free_buffers();
	
}


if (YbusOut->Enable_Memory_Stat == 1) {	
	mkl_peak_mem_usage(MKL_PEAK_MEM_DISABLE);
}

stop = dsecnd();
(*timer)[0] = stop - start;

return status;
 

}
