#include "PowerFlow_NR_MKL.h"
#define ALIGN 64

/* To avoid constantly repeating the part of code that checks different functions' status, using the below macros */
#define CHECK_SPARSE(function)  do {  if(function != SPARSE_STATUS_SUCCESS){ status = 3; goto memory_free; } } while(0) // Check sparse function status
#define CHECK_OPERATION(function)  do {  if(function != 0){ converged = 3; goto memory_free; } } while(0) // check if the operation of an internal function was successful
#define CHECK_SOLVER(function)  do {  if(function != 0){ converged = 4; goto memory_free; } } while(0) // check if the operation of an internal function was successful
#define CHECK_MEM(variable)  do {  if(variable == NULL){ converged = 5; goto memory_free; } } while(0) // check memory allocation
#define CHECK_MEM_JACOBIAN(variable)  do {  if(variable == NULL){ status = 1; goto mem_free; } } while(0) // check memory allocation



// ***********************************************************************************************************************************************
// ************************************************************** DENSE FUNCTIONS ***************************************************************

// ********************************* Sbus ***************************
// This function calculates complex bus power injection
void Sbus(unsigned long long int Num_Buses, unsigned long long int Gen_Num_Buses, double BaseMVA, double *Pd, double *Qd, double *Pg, double *Qg, int *GenBusNo, int *GenStatus, double *Sbus_Real, double *Sbus_Imag, int OMP_Enbale, int OMP_Cores){
	
	unsigned long long int i = 0;
	unsigned long long int Bus_No = 0;

	// Negative of Demand
	unsigned long long int Chunk_Size = round(Num_Buses / OMP_Cores);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Num_Buses,Sbus_Real,Sbus_Imag,Pd,Qd,BaseMVA) if(OMP_Enbale == 1)
	for (i=0;i< Num_Buses;i++){
		Sbus_Real[i] = (-1)*Pd[i]/BaseMVA;
		Sbus_Imag[i] = (-1)*Qd[i]/BaseMVA;
	}

	//Online Generation - Demand
	Chunk_Size = round(Gen_Num_Buses / OMP_Cores);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i,Bus_No) shared(Chunk_Size,Gen_Num_Buses,Sbus_Real,Sbus_Imag,Pg,Qg,GenStatus,GenBusNo,BaseMVA) if(OMP_Enbale == 1)
	for (i=0;i<Gen_Num_Buses;i++){
		Bus_No = GenBusNo[i] - 1;
		Sbus_Real[Bus_No] = (Pg[i]*GenStatus[i]/BaseMVA)+Sbus_Real[Bus_No];
		Sbus_Imag[Bus_No] = (Qg[i]*GenStatus[i]/BaseMVA)+Sbus_Imag[Bus_No];
	}


}



// **************************** Injected Bus Power ********************
// This function calculates the injected bus powers
void Power_Mimatch(MKL_Complex16 *Ybus_Row_, sparse_matrix_t Ybus_CSR, MKL_Complex16 *V_Rec, MKL_Complex16 *YbusV, double *Sbus_re, double *Sbus_img, unsigned long long int Num_Buses, double *Mis_re, double *Mis_img, int OMP_Enbale, int OMP_Cores, int Sparse_Enable) {
	
	unsigned long int i = 0;
	int Chunk_Size = round(Num_Buses / OMP_Cores);
	MKL_Complex16 alpha, beta;
	alpha.real = 1;
	alpha.imag = 0;
	beta.real = 0;
	beta.imag = 0;
	struct matrix_descr descrA;
	descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

	/* Power Mismatch = V .* conj(Ybus * V) - Sbus(Vm) */

	// YbusV = Ybus * V
	if (Sparse_Enable == 0) {
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Num_Buses, 1, Num_Buses, &alpha, Ybus_Row_, Num_Buses, V_Rec, 1, &beta, YbusV, 1);
	}
	else{
		mkl_sparse_z_mm(SPARSE_OPERATION_NON_TRANSPOSE, alpha, Ybus_CSR, descrA, SPARSE_LAYOUT_ROW_MAJOR, V_Rec, 1, Num_Buses, beta, YbusV, Num_Buses);
	}


	// V .* conj(YbusV) - Sbus(Vm)
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size, Num_Buses, YbusV, Mis_re, Mis_img, V_Rec, Sbus_re, Sbus_img) if(OMP_Enbale == 1)
	for (i = 0; i < Num_Buses; i++) {
		Mis_re[i] = Mult_re(V_Rec[i].real, V_Rec[i].imag, YbusV[i].real, -1*YbusV[i].imag) - Sbus_re[i];
		Mis_img[i] = Mult_img(V_Rec[i].real, V_Rec[i].imag, YbusV[i].real, -1*YbusV[i].imag) - Sbus_img[i];
	}

}


// ***************************** Objective Function ****************
// This function builds norm vector [ Mis_re[PV;PQ] ; Mis_img[PQ] ] and calculate the infinity norm (largest number in the vector)
void F_Build(double *Mis_re, double *Mis_img, unsigned long long int Num_Buses, unsigned long long int PV_Elem, unsigned long long int PQ_Elem, int *PVBuses, int *PQBuses, double *F_, double *F_Norm, int OMP_Enbale, int OMP_Cores){
	
	unsigned long long int i = 0, F_index = 0;
	double max_val = 0, max_val1 = 0;
	unsigned long long int Chunk_Size = round(PV_Elem / OMP_Cores);

	// setting scheduling to guided caused problem in linux
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,F_,PV_Elem,Mis_re,PVBuses) reduction(max:max_val) if(OMP_Enbale == 1)
	for (i = 0; i < PV_Elem; i++) {
		F_[i] = Mis_re[PVBuses[i] - 1];
		max_val = max_val > fabs(F_[i]) ? max_val : fabs(F_[i]);
	}


	// setting scheduling to guided caused problem in linux
	Chunk_Size = round(PQ_Elem / OMP_Cores);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i,F_index) shared(Chunk_Size,F_,PQ_Elem,PV_Elem,Mis_re,Mis_img,PQBuses,PVBuses) reduction(max:max_val1) if(OMP_Enbale == 1)
	for (i = 0; i < PQ_Elem; i++) {
		F_index = i + PV_Elem;
		F_[F_index] = Mis_re[PQBuses[i] - 1];
		max_val1 = max_val1 > fabs(F_[F_index]) ? max_val1 : fabs(F_[F_index]);	

		F_index = i + PV_Elem + PQ_Elem;
		F_[F_index] = Mis_img[PQBuses[i] - 1];
		max_val1 = max_val1 > fabs(F_[F_index]) ? max_val1 : fabs(F_[F_index]);	
	}


	if (max_val > max_val1) {
		*F_Norm = max_val;
	}
	else {
		*F_Norm = max_val1;
	}


}

// ***************************** dSbus_dV ***************************
//  Computes partial derivatives of power injection with respect to voltage.
// Calculate Ibus (m x 1) = Ybus (m x m) * V (m x 1)
// Calculate dS_dVm_re = diagV * conj(Ybus * diagVnorm) + conj(diagIbus) * diagVnorm
// Calculate dSbus_dVa = 1j * diagV * conj(diagIbus - Ybus * diagV);

void dS_dV(MKL_Complex16 *Ybus_Row_, MKL_Complex16 *V_Rec, MKL_Complex16 *Ibus, MKL_Complex16 *diagVnorm, 
	MKL_Complex16 *YbusdiagVnorm, MKL_Complex16 *YbusdiagV , double *V_mag, unsigned long long int Num_Buses, MKL_Complex16 *ds_dVa, MKL_Complex16 *ds_dVm, int OMP_Enbale, int OMP_Cores){

	unsigned long long int i = 0, j = 0;
	MKL_Complex16 alpha, beta;
	alpha.real = 1;
	alpha.imag = 0;
	beta.real = 0;
	beta.imag = 0;

	// Ibus = Ybus * V
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Num_Buses, 1, Num_Buses, &alpha, Ybus_Row_, Num_Buses, V_Rec, 1, &beta, Ibus, 1); 

	// Now form diagonal matrices
	unsigned long long int Chunk_Size = round(Num_Buses / OMP_Cores);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Num_Buses,V_Rec,V_mag,Ibus,diagVnorm) if(OMP_Enbale == 1)
	for (i = 0; i < Num_Buses; i++) {
		// Form conjugate complex of diagonal Ibus	
		Ibus[i].imag = Ibus[i].imag * (-1);
		// Form diagnal Vnorm
		diagVnorm[i].real = V_Rec[i].real / V_mag[i];
		diagVnorm[i].imag = V_Rec[i].imag / V_mag[i];
	}

	// *** ds_dVm & ds_dVa***
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i,j) shared(Chunk_Size,Num_Buses,YbusdiagVnorm,Ybus_Row_,V_Rec,diagVnorm,YbusdiagV,Ibus,ds_dVm,ds_dVa) collapse(2) if(OMP_Enbale == 1)
	for (i = 0; i < Num_Buses; i++) {
		for (j = 0; j < Num_Buses; j++) {
			/* YbusdiagVnorm = Ybus * diagVnorm */
			YbusdiagVnorm[i*Num_Buses + j].real = Mult_re(Ybus_Row_[i*Num_Buses + j].real, Ybus_Row_[i*Num_Buses + j].imag, diagVnorm[j].real, diagVnorm[j].imag);
			YbusdiagVnorm[i*Num_Buses + j].imag = (Mult_img(Ybus_Row_[i*Num_Buses + j].real, Ybus_Row_[i*Num_Buses + j].imag, diagVnorm[j].real, diagVnorm[j].imag))*(-1);
			/* YbusdiagV = Ybus * diagV */
			YbusdiagV[i*Num_Buses + j].real = Mult_re(Ybus_Row_[i*Num_Buses + j].real, Ybus_Row_[i*Num_Buses + j].imag, V_Rec[j].real, V_Rec[j].imag) * (-1);
			YbusdiagV[i*Num_Buses + j].imag = Mult_img(Ybus_Row_[i*Num_Buses + j].real, Ybus_Row_[i*Num_Buses + j].imag, V_Rec[j].real, V_Rec[j].imag);
			/* *** ds_dVm = diagV * conj(YbusdiagVnorm) *** */
			ds_dVm[i*Num_Buses + j].real = Mult_re(V_Rec[i].real, V_Rec[i].imag, YbusdiagVnorm[i*Num_Buses + j].real, YbusdiagVnorm[i*Num_Buses + j].imag);
			ds_dVm[i*Num_Buses + j].imag = Mult_img(V_Rec[i].real, V_Rec[i].imag, YbusdiagVnorm[i*Num_Buses + j].real, YbusdiagVnorm[i*Num_Buses + j].imag);
			/* Diagonal Elements */
			if (i == j) {
				ds_dVm[i*Num_Buses + j].real = ds_dVm[i*Num_Buses + j].real + Mult_re(Ibus[i].real, Ibus[i].imag, diagVnorm[j].real, diagVnorm[j].imag);
				ds_dVm[i*Num_Buses + j].imag = ds_dVm[i*Num_Buses + j].imag + Mult_img(Ibus[i].real, Ibus[i].imag, diagVnorm[j].real, diagVnorm[j].imag);
				YbusdiagV[i*Num_Buses + j].real = Ibus[i].real + YbusdiagV[i*Num_Buses + j].real;
				YbusdiagV[i*Num_Buses + j].imag = Ibus[i].imag + YbusdiagV[i*Num_Buses + j].imag;		
			}	
			/* *** ds_dVa = 1j * diagV * ds_dVa *** */
			ds_dVa[i*Num_Buses + j].imag = Mult_re(V_Rec[i].real, V_Rec[i].imag, YbusdiagV[i*Num_Buses + j].real, YbusdiagV[i*Num_Buses + j].imag);
			ds_dVa[i*Num_Buses + j].real = (-1)*Mult_img(V_Rec[i].real, V_Rec[i].imag, YbusdiagV[i*Num_Buses + j].real, YbusdiagV[i*Num_Buses + j].imag);
		}
	}


}



// ********************************* Jacobian Matrix *********************
// This function forms the Jacobian matrix
void Jacobian(MKL_Complex16 *ds_dVm, MKL_Complex16 *ds_dVa, int *PVPQ, unsigned long long int Num_Buses, unsigned long long int PV_Elem, unsigned long long int PQ_Elem, int *PVBuses, int *PQBuses, double *J_, int OMP_Enbale, int OMP_Cores){
	
	unsigned long long int i = 0, j = 0;
	unsigned long long int dVa_real_ind, dVa_imag_ind, dVm_real_ind, dVm_imag_ind;
	unsigned long long int PVPQsize = PV_Elem + PQ_Elem;
	unsigned long long int Jsize = PV_Elem + (PQ_Elem * 2);
	unsigned long long int Chunk_Size = round(PVPQsize / OMP_Cores);


	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i,j,dVa_real_ind,dVa_imag_ind,dVm_real_ind,dVm_imag_ind) shared(Chunk_Size,PVPQsize,J_,Num_Buses,Jsize,PVPQ,PV_Elem,PQ_Elem,PVBuses,PQBuses,ds_dVa,ds_dVm) collapse(2) if(OMP_Enbale == 1)
	for (i=0; i<PVPQsize; i++){ // J = [J1 J2; J3 J4]
		for (j=0;j<PVPQsize;j++){
			dVa_real_ind = (PVPQ[i] - 1)*Num_Buses + (PVPQ[j] - 1);
			dVa_imag_ind = (PQBuses[i] - 1)*Num_Buses + (PVPQ[j] - 1);
			dVm_real_ind = (PVPQ[i] - 1)*Num_Buses + (PQBuses[j] - 1);
			dVm_imag_ind = (PQBuses[i] - 1)*Num_Buses + (PQBuses[j] - 1);

				J_[(i*Jsize) + j] = ds_dVa[dVa_real_ind].real;			

				if (j < PQ_Elem)
					J_[(i*Jsize) + (j + PVPQsize)] = ds_dVm[dVm_real_ind].real;

				if (i < PQ_Elem)
					J_[(i + PVPQsize)*Jsize + j] = ds_dVa[dVa_imag_ind].imag;

				if ((i < PQ_Elem) && (j < PQ_Elem))
					J_[(i + PVPQsize)*Jsize + (j + PVPQsize)] = ds_dVm[dVm_imag_ind].imag;		
		}
	}


}



// ********************************* Update Voltages *********************
// This function updates the voltage vectors
void Voltage_Update(double *V_Mag, double *V_Angl, MKL_Complex16 *V_Rectangular, unsigned long long int Num_of_Buses, unsigned long long int PV_Buses_Num_Elements, int *PVBuses, double *F, unsigned long long int PQ_Buses_Num_Elements, int *PQBuses, int OMP_Enbale, int OMP_Cores) {

	unsigned long long int i, j, Chunk_Size;
	double temp1, temp2;
	
	Chunk_Size = round(PV_Buses_Num_Elements / OMP_Cores);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i,j) shared(Chunk_Size,PV_Buses_Num_Elements,PVBuses,V_Mag,V_Angl,V_Rectangular,F) if(OMP_Enbale == 1)
	for (i = 0; i < PV_Buses_Num_Elements; i++) {
		j = PVBuses[i] - 1; // PV bus number
		V_Angl[j] = V_Angl[j] + F[i] * (-1); // updated angle
		V_Rectangular[j].real = V_Mag[j] * cos(V_Angl[j]);
		V_Rectangular[j].imag = V_Mag[j] * sin(V_Angl[j]);
	}
	
	Chunk_Size = round(PQ_Buses_Num_Elements / OMP_Cores);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i,j) shared(Chunk_Size,PQ_Buses_Num_Elements,PQBuses,PV_Buses_Num_Elements,PVBuses,V_Mag,V_Angl,V_Rectangular,F) if(OMP_Enbale == 1)
	for (i = 0; i < PQ_Buses_Num_Elements; i++) {
		j = PQBuses[i] - 1; // PQ bus number
		V_Angl[j] = V_Angl[j] + F[i + PV_Buses_Num_Elements] * (-1); // update angle
		V_Mag[j] = V_Mag[j] + F[i + PV_Buses_Num_Elements + PQ_Buses_Num_Elements] * (-1); // update magnitude
		V_Rectangular[j].real = V_Mag[j] * cos(V_Angl[j]);
		V_Rectangular[j].imag = V_Mag[j] * sin(V_Angl[j]);
	}
	
	// Update Vm and Va again in case we wrapped around with a negative Vm : V = Vm .* exp(1j * Va) = Vm * [cos(Va) + j*sin(Va)] = Vm*cos(va) + j*Vm*sin(Va) = Vre + j*Vimg
	Chunk_Size = round(Num_of_Buses / OMP_Cores);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i,temp1,temp2) shared(Chunk_Size,Num_of_Buses,V_Mag,V_Angl,V_Rectangular) if(OMP_Enbale == 1)
	for (i = 0; i < Num_of_Buses; i++) {
		temp1 = V_Mag[i] * cos(V_Angl[i]); // real part
		temp2 = V_Mag[i] * sin(V_Angl[i]); // imaginary part
		V_Mag[i] = rec2pol_mag(temp1, temp2);
		V_Angl[i] = rec2pol_ang(temp1, temp2);
		V_Rectangular[i].real = temp1;
		V_Rectangular[i].imag = temp2;
	}

}


// ***********************************************************************************************************************************************
// ************************************************************** SPARSE FUNCTIONS ***************************************************************

MKL_Complex16 find(MKL_Complex16 *VAL, MKL_INT *JA, MKL_INT *IA, unsigned long long int row, unsigned long long int col) {
	unsigned long long int i;
	MKL_Complex16 result;
	result.real = 0;
	result.imag = 0;    
		for (i = JA[row]; i < JA[row + 1]; i++) {
			if (IA[i] == col) {
				result.real = VAL[i].real;
				result.imag = VAL[i].imag;
				break;
			}			
		}
	return result;
}


int pardiso_solver(double *a, MKL_INT *ia, MKL_INT *ja, double *b, double *x, MKL_INT n, int OMP_Enbale, int OMP_Cores) {

	/* Auxiliary variables. */
	MKL_INT i, j;
	double ddum;          /* Double dummy */
	MKL_INT idum;         /* Integer dummy. */
	double res = 0, res0 = 0;
	int status = 0;

	// Real unsymmetric matrix 
	MKL_INT mtype = 11;
	// Real and structurally symmetric
	//MKL_INT mtype = 1;
	// Descriptor of main sparse matrix properties
	struct matrix_descr descrA;
	// Structure with sparse matrix stored in CSR format
	sparse_matrix_t       csrA;
	sparse_operation_t    transA;
	/* RHS and solution vectors. */
	MKL_INT nrhs = 1;     /* Number of right hand sides. */
	/* Internal solver memory pointer pt, */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	/* or void *pt[64] should be OK on both architectures */
	void *pt[64];
	/* Pardiso control parameters. */
	MKL_INT iparm[64];
	MKL_INT maxfct, mnum, phase, error, msglvl;


	/* -------------------------------------------------------------------- */
	/* .. Setup Pardiso control parameters. */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++) iparm[i] = 0;

	iparm[0] = 1;         /* No solver default */
	iparm[1] = 3;         /* Fill-in reordering using nested dissection algorithm */
	iparm[3] = 0;         /* No iterative-direct algorithm */
	iparm[4] = 0;         /* No user fill-in reducing permutation */
	iparm[5] = 0;         /* Write solution into x */
	iparm[6] = 0;         /* Not in use */
	iparm[7] = 2;         /* Max numbers of iterative refinement steps */
	iparm[8] = 0;         /* Not in use */
	iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
	//iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0;        /* Conjugate transposed/transpose solve */
	//iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric = 1) */
	iparm[13] = 0;        /* Output: Number of perturbed pivots */
	iparm[14] = 0;        /* Not in use */
	iparm[15] = 0;        /* Not in use */
	iparm[16] = 0;        /* Not in use */
	iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1;       /* Output: Mflops for LU factorization */
	iparm[19] = 0;        /* Output: Numbers of CG Iterations */
	iparm[23] = 1;        /* Use a two-level factorization algorithm */
	iparm[24] = 2;        /* Use a parallel algorithm for the solve step */
	iparm[34] = 1;        /* Zero indexing */
	maxfct = 1;           /* Maximum number of numerical factorizations. */
	mnum = 1;             /* Which factorization to use. */
	msglvl = 0;           /* Print statistical information  */
	error = 0;            /* Initialize error flag */

	/* -------------------------------------------------------------------- */
	/* .. Initialize the internal solver memory pointer. This is only */
	/* necessary for the FIRST call of the PARDISO solver. */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++)  pt[i] = 0;


	/* -------------------------------------------------------------------- */
	/* Analysis, numerical factorization, solve, iterative refinement.      */
	/* -------------------------------------------------------------------- */
	phase = 13;
	descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
	//descrA.mode = SPARSE_FILL_MODE_UPPER;
	//descrA.diag = SPARSE_DIAG_NON_UNIT;
	iparm[11] = 0;        /* Conjugate transposed/transpose solve */
	PARDISO_64(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
	if (error != 0)
	{
		status = 3;
		goto mem;
	}

	/* -------------------------------------------------------------------- */
	/* Copy the result													    */
	/* -------------------------------------------------------------------- */
	unsigned long long int Chunk_Size = round(n / OMP_Cores);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,n,b,x) if(OMP_Enbale == 1)
	for (i = 0; i < n; i++) b[i] = x[i];

	/* -------------------------------------------------------------------- */
	/* .. Termination and release of memory. */
	/* -------------------------------------------------------------------- */
mem:
	phase = -1;           /* Release internal memory. */
	PARDISO_64(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0) status = 3;

	return status;

}




int dS_dV_Sparse(MKL_Complex16 *ds_dVa_vals, MKL_Complex16 *ds_dVm_vals, MKL_INT *ds_dV_IA, MKL_INT *ds_dV_JA,
	unsigned long long int Ybus_nnz, sparse_matrix_t Ybus_CSR, sparse_matrix_t  diagIbus_CSR, sparse_matrix_t diagV_CSR, sparse_matrix_t diagVnorm_CSR,
	MKL_Complex16 *V_Rec, MKL_Complex16 *Ibus, double *V_mag, unsigned long long int Num_Buses, int OMP_Enbale, int OMP_Cores) {

	int status = 0;
	unsigned long long int i = 0;
	unsigned long long int j = 0;
	MKL_Complex16 alpha, beta;
	alpha.real = 1;
	alpha.imag = 0;
	beta.real = 0;
	beta.imag = 0;
	MKL_Complex16 temp;
	
	sparse_index_base_t indexing = 0;
	
	sparse_matrix_t YbusdiagVnorm_CSR = NULL, YbusdiagV_CSR = NULL;
	sparse_matrix_t  ds_dVm_CSR = NULL, ds_dVa_CSR = NULL;
	sparse_matrix_t ds_dVm_CSR1 = NULL, ds_dVm_CSR2 = NULL;	
	sparse_matrix_t diagIbusMinusYbusdiagV_CSR = NULL;
	
	
	MKL_INT YbusdiagVnorm_row, YbusdiagVnorm_col, YbusdiagV_row, YbusdiagV_col, ds_dVm_row, ds_dVm_col, ds_dVa_row, ds_dVa_col;
	MKL_INT *YbusdiagVnorm_pointerB = NULL, *YbusdiagVnorm_pointerE = NULL, *YbusdiagVnorm_columns = NULL;
	MKL_INT *YbusdiagV_pointerB = NULL, *YbusdiagV_pointerE = NULL, *YbusdiagV_columns = NULL;
	MKL_Complex16  *YbusdiagVnorm_values = NULL, *YbusdiagV_values = NULL;
	
	MKL_INT *ds_dVm_JA = NULL, *ds_dVm_pointerE = NULL, *ds_dVm_IA = NULL;
	MKL_Complex16 *ds_dVm_values = NULL;
	MKL_INT *ds_dVa_JA = NULL, *ds_dVa_pointerE = NULL, *ds_dVa_IA = NULL;
	MKL_Complex16 *ds_dVa_values = NULL;


	struct matrix_descr descr_type_gen, descrA, descrB;
	descr_type_gen.type = SPARSE_MATRIX_TYPE_GENERAL;  //or SPARSE_MATRIX_TYPE_SYMMETRIC or ...
	descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
	descrB.type = SPARSE_MATRIX_TYPE_DIAGONAL;

	// Ibus = Ybus * V (result is a dense vector)
	mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, Ybus_CSR, descr_type_gen, V_Rec, beta, Ibus); // multiply a sparse matrix by a dense vector and return the result in dense vector format
	
	// Now form diagonal matrices	
	unsigned long long int Chunk_Size = round(Num_Buses / OMP_Cores);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i,temp) shared(Chunk_Size,Num_Buses,V_mag,V_Rec,diagV_CSR,diagIbus_CSR,diagVnorm_CSR,Ibus) if(OMP_Enbale == 1)
	for (i = 0; i < Num_Buses; i++) {
		// Form diagonal V (conervt polar form to rectangular first)
		mkl_sparse_z_set_value(diagV_CSR, i, i, V_Rec[i]);
		
		// Form conjugate complex of diagonal Ibus
		temp.real = Ibus[i].real;
		temp.imag = Ibus[i].imag * (-1);
		mkl_sparse_z_set_value(diagIbus_CSR, i, i, temp);
		
		// Form diagnal Vnorm
		temp.real = V_Rec[i].real / V_mag[i]; // real V/mag V
		temp.imag = V_Rec[i].imag / V_mag[i]; // imag V/mag v
		mkl_sparse_z_set_value(diagVnorm_CSR, i, i, temp);		
	}

	
	// ****** ds_dVm ******
	// YbusdiagVnorm = Ybus * diagVnorm		
	CHECK_SPARSE( mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, Ybus_CSR, diagVnorm_CSR, &YbusdiagVnorm_CSR) );
	// YbusdiagV = Ybus * diagV (This is used in ds_dVa)		
	CHECK_SPARSE( mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, Ybus_CSR, diagV_CSR, &YbusdiagV_CSR) );
	// Convert YbusdiagVnorm_CSR and YbusdiagV_CSR handles to 4-array CSR
	CHECK_SPARSE( mkl_sparse_z_export_csr(YbusdiagVnorm_CSR, &indexing, &YbusdiagVnorm_row, &YbusdiagVnorm_col, &YbusdiagVnorm_pointerB, &YbusdiagVnorm_pointerE, &YbusdiagVnorm_columns, &YbusdiagVnorm_values) );
	CHECK_SPARSE( mkl_sparse_z_export_csr(YbusdiagV_CSR, &indexing, &YbusdiagV_row, &YbusdiagV_col, &YbusdiagV_pointerB, &YbusdiagV_pointerE, &YbusdiagV_columns, &YbusdiagV_values) );
	// Calculate complex conjugate of YbusdiagVnorm and YbusdiagV
	Chunk_Size = round(Ybus_nnz / OMP_Cores);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Ybus_nnz,YbusdiagVnorm_values,YbusdiagV_values) if(OMP_Enbale == 1)
	for (i = 0; i < Ybus_nnz; i++) {
		YbusdiagVnorm_values[i].imag = YbusdiagVnorm_values[i].imag * (-1);
		YbusdiagV_values[i].imag = YbusdiagV_values[i].imag * (-1);		
	}	
	CHECK_SPARSE( mkl_sparse_z_create_csr(&YbusdiagVnorm_CSR, SPARSE_INDEX_BASE_ZERO, Num_Buses, Num_Buses, YbusdiagVnorm_pointerB, YbusdiagVnorm_pointerB + 1, YbusdiagVnorm_columns, YbusdiagVnorm_values) );
	CHECK_SPARSE( mkl_sparse_z_create_csr(&YbusdiagV_CSR, SPARSE_INDEX_BASE_ZERO, Num_Buses, Num_Buses, YbusdiagV_pointerB, YbusdiagV_pointerB + 1, YbusdiagV_columns, YbusdiagV_values) );
	// calculate ds_dVm_1 = diagV * conj(YbusdiagVnorm)	
	CHECK_SPARSE( mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, diagV_CSR, YbusdiagVnorm_CSR, &ds_dVm_CSR1) );
	// calculate ds_dVm_2 = conj(diagIbus)*diagVnorm
	CHECK_SPARSE( mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, diagIbus_CSR, diagVnorm_CSR, &ds_dVm_CSR2) );
	// calculate ds_dVm = ds_dVm_1 + ds_dVm_2
	alpha.real = 1;
	alpha.imag = 0;
	CHECK_SPARSE( mkl_sparse_z_add(SPARSE_OPERATION_NON_TRANSPOSE, ds_dVm_CSR1, alpha, ds_dVm_CSR2, &ds_dVm_CSR) );
	CHECK_SPARSE( mkl_sparse_z_export_csr(ds_dVm_CSR, &indexing, &ds_dVm_row, &ds_dVm_col, &ds_dVm_JA, &ds_dVm_pointerE, &ds_dVm_IA, &ds_dVm_values) );
		


	// ******* ds_dVa ******
	// diagIbusMinusYbusdiagV = conj(diagIbus) - conj(YbusdiagV)	
	alpha.real = -1;
	alpha.imag = 0;
	CHECK_SPARSE( mkl_sparse_z_add(SPARSE_OPERATION_NON_TRANSPOSE, YbusdiagV_CSR, alpha, diagIbus_CSR, &diagIbusMinusYbusdiagV_CSR) );
	// ds_dVa = 1j * diagV * ds_dVa : mkl_sparse_spmm does not support multiplication by scalar.. so we have to extract the values and multiply them by 1j.	
	CHECK_SPARSE( mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, diagV_CSR, diagIbusMinusYbusdiagV_CSR, &ds_dVa_CSR) );
	// Convert ds_dVa_CSR to 4-array version
	CHECK_SPARSE( mkl_sparse_z_export_csr(ds_dVa_CSR, &indexing, &ds_dVa_row, &ds_dVa_col, &ds_dVa_JA, &ds_dVa_pointerE, &ds_dVa_IA, &ds_dVa_values) );
	// x 1j operation
	Chunk_Size = round(Ybus_nnz / OMP_Cores);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Ybus_nnz,ds_dVa_values,ds_dVm_values,ds_dVa_vals,ds_dVm_vals,ds_dV_IA,ds_dVa_IA) if(OMP_Enbale == 1)
	for (i = 0; i < Ybus_nnz; i++) {
		ds_dVa_vals[i].real = ds_dVa_values[i].imag * (-1);
		ds_dVa_vals[i].imag = ds_dVa_values[i].real;
		ds_dVm_vals[i].real = ds_dVm_values[i].real;
		ds_dVm_vals[i].imag = ds_dVm_values[i].imag;
		ds_dV_IA[i] = ds_dVa_IA[i];
	}
	
	Chunk_Size = round((Num_Buses+1) / OMP_Cores);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Num_Buses,ds_dV_JA, ds_dVa_JA) if(OMP_Enbale == 1)
	for (i = 0; i < (Num_Buses+1); i++) {
		ds_dV_JA[i] = ds_dVa_JA[i];
	}

memory_free:
	mkl_sparse_destroy(YbusdiagVnorm_CSR);
	mkl_sparse_destroy(YbusdiagV_CSR);
	mkl_sparse_destroy(ds_dVm_CSR1);
	mkl_sparse_destroy(ds_dVm_CSR2);
	mkl_sparse_destroy(diagIbusMinusYbusdiagV_CSR);
	mkl_sparse_destroy(ds_dVm_CSR);
	mkl_sparse_destroy(ds_dVa_CSR);

	return status;
}


long long int BinarySearch(int *arr, long long int value, long long int lower, long long int upper) {
	long long int middle = 0;
	long long int index = -1;
	
	while (lower <= upper) {
		middle = (lower + upper) / 2;
		
		// if the elemnt is present at the middle 
		if (arr[middle] == value) {
			index = middle;
			break;
		} //return middle;

		// if element is smaller than middle, search the left sub array
		if (arr[middle] > value) {
			upper = middle - 1;
		}
		// if element is larger than middle, search the right sub array
		else {
			lower = middle + 1;
		}
		 
	}
	
	return index;
}


// ***** Check if Value is inside Ref_Buses-1: used in Jacobian_Sparse_Initialize()  ********
int CheckRef(unsigned long long int value, unsigned long long int *Ref_Buses, unsigned long long int Ref_Buses_Size) {
	unsigned long long int i;
	int status = 0;
	for (i = 0; i < Ref_Buses_Size; i++) {
		if ((Ref_Buses[i] - 1) == value) {
			status = 1; // value meets the condition
			break;
		}
	}
	return status;
}

// ********************************* Initialize Jacobian Matrix - Sparse *********************
// This function initializes the sparse Jacobian matrices
int Jacobian_Sparse_Initialize(double **J_val, MKL_INT **J_IA, MKL_INT **J_JA, MKL_INT **J_IA_Array, MKL_INT **J_IA_Position, unsigned long long int *Jnnz_, unsigned long long int *Ref_Buses, unsigned long long int Ref_Buses_Size, unsigned long long int Ybusnnz, MKL_Complex16 *ds_dVa_vals, MKL_Complex16 *ds_dVm_vals, MKL_INT *ds_dV_IA, MKL_INT *ds_dV_JA, int *PVPQ, unsigned long long int Num_Buses, unsigned long long int PV_Elem, unsigned long long int PQ_Elem, int *PQBuses, int OMP_Enbale, int OMP_Cores) {
	unsigned long long int i = 0, j = 0, k = 0, l = 0;
	unsigned long long int PVPQsize = PV_Elem + PQ_Elem;
	unsigned long long int Jsize = PV_Elem + (PQ_Elem * 2);
	unsigned long long int row = 0, col = 0;
	int PQ_true = 0, status = 0;
	unsigned long long int Chunk_Size = round(Num_Buses / OMP_Cores);
	unsigned long long int Jnnz = 0; // totoal number of non-zero elements in J
	MKL_INT *Jnnz_row = NULL; // this array contains the number of non-zero elemenst in each row of J	
	long long int temp = 0;
		
	/* we need to find the number of nonzero elements in each row of J in order to correctly allocate memory for sparse matrices*/	
		Jnnz_row = (MKL_INT *)mkl_calloc((size_t)(Jsize), sizeof(MKL_INT), ALIGN);
		CHECK_MEM_JACOBIAN(Jnnz_row);
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i,j,row,col,k,PQ_true,temp) shared(Chunk_Size,Ref_Buses,Ref_Buses_Size,PVPQsize,Num_Buses,Jsize,PVPQ,PV_Elem,PQ_Elem,PQBuses,ds_dVa_vals, ds_dVm_vals, ds_dV_JA, ds_dV_IA) reduction(+:Jnnz_row[:Jsize],Jnnz) if(OMP_Enbale == 1)		
		for (i = 0; i < Num_Buses; i++) {
			for (j = ds_dV_JA[i]; j < ds_dV_JA[i + 1]; j++) {

				if (CheckRef(i,Ref_Buses,Ref_Buses_Size) == 0 && CheckRef(ds_dV_IA[j], Ref_Buses, Ref_Buses_Size) == 0) { 
					PQ_true = 0;
					
					//j12 [PVPQ:PQ]
					//J_[(row*Jsize) + (col + PVPQsize)] = ds_dVm_vals[j].real;										
					temp = BinarySearch(PVPQ, i + 1, 0, PV_Elem - 1);
					if (temp == -1) temp = BinarySearch(PVPQ, i + 1, PV_Elem, PVPQsize - 1);
					row = temp;

					temp = BinarySearch(PVPQ, ds_dV_IA[j] + 1, 0, PV_Elem - 1);
					if (temp == -1) temp = BinarySearch(PVPQ, ds_dV_IA[j] + 1, PV_Elem, PVPQsize - 1);
					col = temp;

					Jnnz_row[row] = Jnnz_row[row] + 1;
					Jnnz = Jnnz + 1;

					//j12 [PVPQ:PQ]
					//J_[(row*Jsize) + (col + PVPQsize)] = ds_dVm_vals[j].real;
					temp = BinarySearch(PQBuses, ds_dV_IA[j] + 1, 0, PQ_Elem - 1);
					if ( temp != (-1) ) {
						col = temp;						
						Jnnz_row[row] = Jnnz_row[row] + 1;
						Jnnz = Jnnz + 1;
					}
					

					//j21 [PQ:PVPQ]
					//J_[(row + PVPQsize)*Jsize + col] = ds_dVa_vals[j].imag;
					PQ_true = 0;
					temp = BinarySearch(PQBuses, i + 1, 0, PQ_Elem - 1);
					if (temp != (-1)) {
						row = temp;
						PQ_true = 1;
					}

					if (PQ_true == 1) {
						// First search the [0 PV_Elements] of PVPQ array
						temp = BinarySearch(PVPQ, ds_dV_IA[j] + 1, 0, PV_Elem - 1);
						// If not found, search the remaining PVPQ
						if (temp == -1) temp = BinarySearch(PVPQ, ds_dV_IA[j] + 1, PV_Elem, PVPQsize - 1);
						
						if (temp != (-1)) {
							col = temp;
							Jnnz_row[row + PVPQsize] = Jnnz_row[row + PVPQsize] + 1;
							Jnnz = Jnnz + 1;
						}												
					}
							   										

					//j22 [PQ:PQ]				
					//J_[(row + PVPQsize)*Jsize + (col + PVPQsize)] = ds_dVm_vals[j].imag;
					if (PQ_true == 1) {
						temp = BinarySearch(PQBuses, ds_dV_IA[j] + 1, 0, PQ_Elem - 1);
						if (temp != (-1)) {
							col = temp;							
							Jnnz_row[row + PVPQsize] = Jnnz_row[row + PVPQsize] + 1;
							Jnnz = Jnnz + 1;
						}						
					}
															

				} // end of ref. check


			} // end of j loop
		} //end of i loop

		
		// Jnnz_ is used in the main function
		(*Jnnz_) = Jnnz;
		
		// allocate memory for sparse matrices
		*J_val = (double *)calloc((size_t)(Jnnz), sizeof(double)); // Values
		CHECK_MEM_JACOBIAN(*J_val);
		*J_IA = (MKL_INT *)mkl_calloc((size_t)(Jnnz), sizeof(MKL_INT), ALIGN); // Column index	
		CHECK_MEM_JACOBIAN(*J_IA);
		*J_JA = (MKL_INT *)mkl_calloc((size_t)(Jsize+1), sizeof(MKL_INT), ALIGN); // Row pointer	
		CHECK_MEM_JACOBIAN(*J_JA);
		*J_IA_Array = (MKL_INT *)mkl_calloc((size_t)(Jnnz), sizeof(MKL_INT), ALIGN);
		CHECK_MEM_JACOBIAN(*J_IA_Array);
		*J_IA_Position = (MKL_INT *)mkl_calloc((size_t)(Jnnz), sizeof(MKL_INT), ALIGN);						
		CHECK_MEM_JACOBIAN(*J_IA_Position);
		
		// Set J_JA values
		(*J_JA)[0] = 0;				
		#pragma nounroll
		for (i = 1; i < Jsize+1; i++) {			
			(*J_JA)[i] = (*J_JA)[i-1] + Jnnz_row[i-1];			
		}

		// Set J_IA_Temp values to -1 for locking purpose
		Chunk_Size = round(Jnnz / OMP_Cores);
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Jnnz,J_IA) if(OMP_Enbale == 1)
		for (i = 0; i < Jnnz + 1; i++) {
			(*J_IA)[i] = -1;
		}
		
		// OMP lock declaration and initializing
		//omp_lock_t lock[Jnnz];
		omp_lock_t *lock = (omp_lock_t *)malloc(Jnnz * sizeof(omp_lock_t));
		for (i = 0; i < Jnnz; i++) {
			omp_init_lock(&(lock[i]));
		}

				
		/* now find the column number of non-zero elements in each row */
		Chunk_Size = round(Num_Buses / OMP_Cores);
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i,j,k,l,row,col,PQ_true,temp) shared(Chunk_Size,lock,Ref_Buses,Ref_Buses_Size,PVPQsize,Num_Buses,J_IA_Array,J_IA_Position,J_IA,J_JA,Jsize,PVPQ,PV_Elem,PQ_Elem,PQBuses,ds_dVa_vals, ds_dVm_vals, ds_dV_JA, ds_dV_IA) if(OMP_Enbale == 1)
		for (i = 0; i < Num_Buses; i++) {
			for (j = ds_dV_JA[i]; j < ds_dV_JA[i + 1]; j++) {

				if (CheckRef(i, Ref_Buses, Ref_Buses_Size) == 0 && CheckRef(ds_dV_IA[j], Ref_Buses, Ref_Buses_Size) == 0) {
					PQ_true = 0;					
	
					//j11 [PVPQ:PVPQ]
					//J_[(row*Jsize) + col] = ds_dVa_vals[j].real;				
					temp = BinarySearch(PVPQ, i + 1, 0, PV_Elem - 1);
					if (temp == -1) temp = BinarySearch(PVPQ, i + 1, PV_Elem, PVPQsize - 1);
					row = temp;

					temp = BinarySearch(PVPQ, ds_dV_IA[j] + 1, 0, PV_Elem - 1);
					if (temp == -1) temp = BinarySearch(PVPQ, ds_dV_IA[j] + 1, PV_Elem, PVPQsize - 1);
					col = temp;

					for (l = (*J_JA)[row]; l < (*J_JA)[row+1]; l++) {						
						omp_set_lock(&(lock[l]));
						if ((*J_IA)[l] == -1) { // if J_IA is not set						
							(*J_IA)[l] = col;
							(*J_IA_Array)[l] = 1;
							(*J_IA_Position)[l] = j;
							omp_unset_lock(&(lock[l]));
							break;
						}													
						omp_unset_lock(&(lock[l]));						
					}

				
					//j12 [PVPQ:PQ]	
					//J_[(row*Jsize) + (col + PVPQsize)] = ds_dVm_vals[j].real;
					temp = BinarySearch(PQBuses, ds_dV_IA[j] + 1, 0, PQ_Elem - 1);
					if (temp != (-1)) {
						col = temp;
						for (l = (*J_JA)[row]; l < (*J_JA)[row + 1]; l++) {
							omp_set_lock(&(lock[l]));
							if ((*J_IA)[l] == -1) { // if J_IA is not set						
								(*J_IA)[l] = (col + PVPQsize);
								(*J_IA_Array)[l] = 3;
								(*J_IA_Position)[l] = j;
								omp_unset_lock(&(lock[l]));
								break;
							}
							omp_unset_lock(&(lock[l]));
						}
					}
												

					//j21 [PQ:PVPQ]
					//J_[(row + PVPQsize)*Jsize + col] = ds_dVa_vals[j].imag;
					PQ_true = 0;
					temp = BinarySearch(PQBuses, i + 1, 0, PQ_Elem - 1);
					if (temp != (-1)) {
						row = temp;
						PQ_true = 1;
					}

					if (PQ_true == 1) {
						// First search the [0 PV_Elements] of PVPQ array
						temp = BinarySearch(PVPQ, ds_dV_IA[j] + 1, 0, PV_Elem - 1);
						// If not found, search the remaining PVPQ
						if (temp == -1) temp = BinarySearch(PVPQ, ds_dV_IA[j] + 1, PV_Elem, PVPQsize - 1);

						if (temp != (-1)) {
							col = temp;

							for (l = (*J_JA)[row + PVPQsize]; l < (*J_JA)[row + PVPQsize + 1]; l++) {
								omp_set_lock(&(lock[l]));
								if ((*J_IA)[l] == -1) { // if J_IA is not set						
									(*J_IA)[l] = col;
									(*J_IA_Array)[l] = 2;
									(*J_IA_Position)[l] = j;
									omp_unset_lock(&(lock[l]));
									break;
								}
								omp_unset_lock(&(lock[l]));
							}

						}
					}
					
					
					//j22 [PQ:PQ]
					//J_[(row + PVPQsize)*Jsize + (col + PVPQsize)] = ds_dVm_vals[j].imag;
					if (PQ_true == 1) {
						temp = BinarySearch(PQBuses, ds_dV_IA[j] + 1, 0, PQ_Elem - 1);
						if (temp != (-1)) {
							col = temp;
							for (l = (*J_JA)[row + PVPQsize]; l < (*J_JA)[row + PVPQsize + 1]; l++) {
								omp_set_lock(&(lock[l]));
								if ((*J_IA)[l] == -1) { // if J_IA is not set						
									(*J_IA)[l] = (col + PVPQsize);
									(*J_IA_Array)[l] = 4;
									(*J_IA_Position)[l] = j;
									omp_unset_lock(&(lock[l]));
									break;
								}
								omp_unset_lock(&(lock[l]));
							}

						}
					}

				

				} // end of ref. check


			} // end of j loop
		} //end of i loop
		
				
		/* We now need to sort the values in IA arrays*/
		Chunk_Size = round(Jsize / OMP_Cores);
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i,j,l,temp) shared(Chunk_Size,Jsize,J_IA_Array,J_IA_Position,J_IA,J_JA) if(OMP_Enbale == 1)
		for (l = 0; l < Jsize; l++) {

			for (i = (*J_JA)[l]; i < (*J_JA)[l + 1]; i++) {
				for (j = (*J_JA)[l]; j < (*J_JA)[l + 1]-1; j++) {
					if ( (*J_IA)[j] > (*J_IA)[j + 1]) {
						//J_IA
						temp = (*J_IA)[j];
						(*J_IA)[j] = (*J_IA)[j + 1];
						(*J_IA)[j + 1] = temp;
						//J_IA_Array
						temp = (*J_IA_Array)[j];
						(*J_IA_Array)[j] = (*J_IA_Array)[j + 1];
						(*J_IA_Array)[j + 1] = temp;
						//J_IA_Position
						temp = (*J_IA_Position)[j];
						(*J_IA_Position)[j] = (*J_IA_Position)[j + 1];
						(*J_IA_Position)[j + 1] = temp;
					}
					
				}						
			}
		
		
		}

mem_free:
	free(lock);
	mkl_free(Jnnz_row);	
	return status;
}


// ********************************* Jacobian Matrix - Sparse *********************
// This function forms the Jacobian matrix
void Jacobian_Sparse(double *J_val, MKL_INT *J_IA_Position, MKL_INT *J_IA_Array, unsigned long long int Jnnz, MKL_Complex16 *ds_dVa_vals, MKL_Complex16 *ds_dVm_vals, int OMP_Enbale, int OMP_Cores) {

	unsigned long long int i = 0, temp = 0;
	unsigned long long int Chunk_Size = round(Jnnz / OMP_Cores);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i,temp) shared(Chunk_Size,Jnnz,J_IA_Array,J_IA_Position,J_val,ds_dVa_vals,ds_dVm_vals) if(OMP_Enbale == 1)
	for (i = 0; i < Jnnz; i++) {
		temp = J_IA_Position[i]; //position in dvx array
		if (J_IA_Array[i] == 1) J_val[i] = ds_dVa_vals[temp].real; //dVa.real					
		if (J_IA_Array[i] == 2) J_val[i] = ds_dVa_vals[temp].imag; //dVa.imag
		if (J_IA_Array[i] == 3) J_val[i] = ds_dVm_vals[temp].real; //dVm.real					
		if (J_IA_Array[i] == 4) J_val[i] = ds_dVm_vals[temp].imag; //dVm.imag				
	}



}

// ***********************************************************************************************************************************************
// ************************************************************** Main Function ******************************************************************

// This function solves the power flow using full Newton-Raphson method
int PF_NR_MKL(MKL_Complex16 *Ybus_Row, MKL_Complex16 *Ybus_val, MKL_INT *Ybus_JA, MKL_INT *Ybus_IA, double *V_Mag_Init, double *V_Angl_Init, unsigned long long int Num_of_Buses,
	int *PVBuses, int *PQBuses,
	int *GenBusNo, int *GenStatus, double BaseMVA, double *Pd, double *Qd, double *Pg, double *Qg, double *Vg,
	unsigned long long int Gen_Buses_Num_Elements, unsigned long long int PV_Buses_Num_Elements, unsigned long long int PQ_Buses_Num_Elements, unsigned long long int *Ref_Buses, unsigned long long int Ref_Buses_Size,
	double *V_Mag, double *V_Angl,
	int FlatStart, int Max_iter, double Tol, int OMP_Enbale, int OMP_Cores, int Sparse_Enable, int Enable_Memory_Stat, double *Peak_Mem_Usage, double *timer) {
	
	/* **** Initial Variables (used in both sparse and dense functions) **** */
	double timer_stop[10], timer_start[10]; // Holds the execution time of different parts
		
	timer_start[0] = dsecnd(); // Total execution time
	if (Enable_Memory_Stat == 1) mkl_peak_mem_usage(MKL_PEAK_MEM_ENABLE);

	timer_start[1] = dsecnd(); // Memory allocation time
	unsigned long long int Chunk_Size = round(Num_of_Buses / OMP_Cores);
	unsigned long long int i = 0, j = 0, k = 0;
	int converged = 6; // Convergence flag (0:success)
	int iter = 0; // iteration counter
	unsigned long long int temp = 0; // dummy
	double temp1 = 0; // dummy
	double temp2 = 0; // dummy

	double *Sbus_Real = NULL;
	double *Sbus_Imag = NULL;
	double *Mis_Real = NULL;
	double *Mis_Imag = NULL;
	double *F = NULL;
	double *J = NULL;
	int *PVPQ = NULL;
	MKL_INT *ipiv = NULL;

	double F_Norm = 0; // contains the infinity norm of vector F
	unsigned long long int Jsize = PV_Buses_Num_Elements + (PQ_Buses_Num_Elements * 2);
	unsigned long long int PVPQsize = PV_Buses_Num_Elements + PQ_Buses_Num_Elements;		
	Sbus_Real = malloc((size_t)Num_of_Buses * sizeof(double)); // real part of complex bus power injection (P)
	CHECK_MEM(Sbus_Real);
	Sbus_Imag = malloc((size_t)Num_of_Buses * sizeof(double)); // imaginary part of complex bus power injection (Q)
	CHECK_MEM(Sbus_Imag);
	Mis_Real = malloc((size_t)Num_of_Buses * sizeof(double)); // real part of power mismatch vector
	CHECK_MEM(Mis_Real);
	Mis_Imag = malloc((size_t)Num_of_Buses * sizeof(double)); // imaginary part of power mismatch vector
	CHECK_MEM(Mis_Imag);
	F = malloc((size_t)Jsize * sizeof(double)); // update state: imaginary part of power mismatch vector
	CHECK_MEM(F);	
	PVPQ = malloc((size_t)PVPQsize * sizeof(int)); // contains all PV and PQ bus numbers [PV;PQ]
	CHECK_MEM(PVPQ);
	ipiv = (MKL_INT *)mkl_malloc((size_t)Jsize * sizeof(MKL_INT), ALIGN);
	CHECK_MEM(ipiv);
	
	// These variable are used in dSbus_dV function
	// Partial derivatives of power injection with respect to voltage magnitude and angle (for all buses)
	MKL_Complex16 *V_Rectangular = NULL; // V in rectangular format
	MKL_Complex16 *Ibus = NULL;
	MKL_Complex16 *diagVnorm = NULL;
	MKL_Complex16 *YbusdiagVnorm = NULL;
	MKL_Complex16 *YbusdiagV = NULL;
	MKL_Complex16 *ds_dVm = NULL;
	MKL_Complex16 *ds_dVa = NULL;
	MKL_Complex16 *YbusV = NULL; // This is used in Power_Mismatch
	V_Rectangular = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Buses, sizeof(MKL_Complex16), ALIGN);
	CHECK_MEM(V_Rectangular);
	Ibus = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Buses, sizeof(MKL_Complex16), ALIGN); // Ibus
	CHECK_MEM(Ibus);
	YbusV = (MKL_Complex16 *)mkl_calloc((size_t)Num_of_Buses, sizeof(MKL_Complex16), ALIGN);
	CHECK_MEM(YbusV);

	// Sparse NULL variables
	unsigned long long int Ybus_nnz = 0; // number of non-zero elements in Ybus
	unsigned long long int Jnnz = 0; //number of non-zero elements in J (size of J_IA, J_Val)
	double  *J_val = NULL;
	MKL_INT *J_IA = NULL;
	MKL_INT *J_JA = NULL;
	MKL_INT *J_IA_Array = NULL; // 1:dVa_real 2:dva_imag 3:dvm_real 4:dvm_imag
	MKL_INT *J_IA_Position = NULL; // Location in ds_dVx_vals array
	double *F_dummy = NULL;
	MKL_Complex16 *diagIbus_val = NULL; // Diagonal Matrix
	MKL_Complex16 *diagV_val = NULL; // Diagonal Matrix
	MKL_Complex16 *diagVnorm_val = NULL; // Diagonal Matrix
	MKL_Complex16 *J_Sparse = NULL;
	sparse_matrix_t  Ybus_CSR = NULL, diagV_CSR = NULL, diagVnorm_CSR = NULL, diagIbus_CSR = NULL;
	MKL_INT *IA = NULL;
	MKL_INT *JA = NULL;
	MKL_Complex16 *ds_dVm_vals = NULL, *ds_dVa_vals = NULL;
	MKL_INT *ds_dV_IA = NULL, *ds_dV_JA = NULL;


	if (Sparse_Enable == 1) {
		// ************* Sparse Variables ****************
		// *** Check input data
		CHECK_MEM(Ybus_val);
		CHECK_MEM(Ybus_IA);
		CHECK_MEM(Ybus_JA);

		// *** Number of non-zero elements in Ybus
		Ybus_nnz = Ybus_JA[Num_of_Buses];

		// *** General IA and JA for diagonal matrices of size (Num_of_Buses x Num_of_Buses) ***
		JA = (MKL_INT *)mkl_malloc((size_t)(Num_of_Buses + 1) * sizeof(MKL_INT), ALIGN); // Row pointer
		CHECK_MEM(JA);
		IA = (MKL_INT *)mkl_malloc((size_t)Ybus_nnz * sizeof(MKL_INT), ALIGN); // Column index
		CHECK_MEM(IA);
		Chunk_Size = round(Num_of_Buses / OMP_Cores);	
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Num_of_Buses, IA, JA) if(OMP_Enbale == 1)
		for (i = 0; i < (Num_of_Buses); i++) {
			IA[i] = i;
			JA[i] = i;
		}
		JA[Num_of_Buses] = Num_of_Buses;

		// ** These variable are used in DS_DV_Sparse ***
		diagIbus_val = (MKL_Complex16 *)mkl_malloc((size_t)Num_of_Buses * sizeof(MKL_Complex16), ALIGN); // Diagonal Matrix
		CHECK_MEM(diagIbus_val);
		diagV_val = (MKL_Complex16 *)mkl_malloc((size_t)Num_of_Buses * sizeof(MKL_Complex16), ALIGN); // Diagonal Matrix
		CHECK_MEM(diagV_val);
		diagVnorm_val = (MKL_Complex16 *)mkl_malloc((size_t)Num_of_Buses * sizeof(MKL_Complex16), ALIGN); // Diagonal Matrix
		CHECK_MEM(diagVnorm_val);
		ds_dVm_vals = (MKL_Complex16 *)mkl_malloc((size_t)Ybus_nnz * sizeof(MKL_Complex16), ALIGN); // Values
		CHECK_MEM(ds_dVm_vals);
		ds_dVa_vals = (MKL_Complex16 *)mkl_malloc((size_t)Ybus_nnz * sizeof(MKL_Complex16), ALIGN); // Values
		CHECK_MEM(ds_dVa_vals);
		ds_dV_IA = (MKL_INT *)mkl_malloc((size_t)Ybus_nnz * sizeof(MKL_INT), ALIGN); // Column index
		CHECK_MEM(ds_dV_IA);
		ds_dV_JA = (MKL_INT *)mkl_malloc((size_t)(Num_of_Buses + 1) * sizeof(MKL_INT), ALIGN); // Row pointer
		CHECK_MEM(ds_dV_JA);
		
		// Create sparse handles	
		mkl_sparse_z_create_csr(&Ybus_CSR, SPARSE_INDEX_BASE_ZERO, Num_of_Buses, Num_of_Buses, Ybus_JA, Ybus_JA + 1, Ybus_IA, Ybus_val);
		CHECK_MEM(Ybus_CSR);
		mkl_sparse_z_create_csr(&diagIbus_CSR, SPARSE_INDEX_BASE_ZERO, Num_of_Buses, Num_of_Buses, JA, JA + 1, IA, diagIbus_val);
		CHECK_MEM(diagIbus_CSR);
		mkl_sparse_z_create_csr(&diagV_CSR, SPARSE_INDEX_BASE_ZERO, Num_of_Buses, Num_of_Buses, JA, JA + 1, IA, diagV_val);
		CHECK_MEM(diagV_CSR);
		mkl_sparse_z_create_csr(&diagVnorm_CSR, SPARSE_INDEX_BASE_ZERO, Num_of_Buses, Num_of_Buses, JA, JA + 1, IA, diagVnorm_val);
		CHECK_MEM(diagVnorm_CSR);

		// This is used in Pardiso solver
		F_dummy = malloc((size_t)Jsize * sizeof(double));
		CHECK_MEM(F_dummy);

		// *********** Sparse Variables END *************
	}
	else {	
		// ************* Dense Variables ****************
		CHECK_MEM(Ybus_Row);
		J = (double *)calloc((size_t)(Jsize*Jsize), sizeof(double)); //Jacobian Matrix
		CHECK_MEM(J);
		ds_dVm = (MKL_Complex16 *)mkl_calloc((size_t)(Num_of_Buses*Num_of_Buses), sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(ds_dVm);
		ds_dVa = (MKL_Complex16 *)mkl_calloc((size_t)(Num_of_Buses*Num_of_Buses), sizeof(MKL_Complex16), ALIGN);
		CHECK_MEM(ds_dVa);
		diagVnorm = (MKL_Complex16 *)mkl_calloc((size_t)(Num_of_Buses), sizeof(MKL_Complex16), ALIGN); //Diagonal Vnorm : V ./ abs(V) , Vnorm.real = 1 
		CHECK_MEM(diagVnorm);
		YbusdiagVnorm = (MKL_Complex16 *)mkl_calloc((size_t)(Num_of_Buses*Num_of_Buses), sizeof(MKL_Complex16), ALIGN); //Ybus * diagVnorm	
		CHECK_MEM(YbusdiagVnorm);
		YbusdiagV = (MKL_Complex16 *)mkl_calloc((size_t)(Num_of_Buses*Num_of_Buses), sizeof(MKL_Complex16), ALIGN); // conj(Ybus*diagV)	
		CHECK_MEM(YbusdiagV);
		// ************* Dense Variables END ****************
	}
	timer_stop[1] = dsecnd(); // Memory allocation time
	timer[1] = (timer_stop[1] - timer_start[1]); // Memory allocation time
	


/* **** V, PVPQ, and F(x0) time **** */
	timer_start[2] = dsecnd(); // Initializing variables time

	// Initialize V
	if (FlatStart == 1){
	//  V0 = 1 < 0
		Chunk_Size = round(Num_of_Buses / OMP_Cores);
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Num_of_Buses,V_Mag,V_Angl,V_Rectangular) if(OMP_Enbale == 1)
		for (i=0; i<Num_of_Buses; i++){
			V_Mag[i] = 1;
			V_Angl[i] = 0;
			V_Rectangular[i].real = 1;
			V_Rectangular[i].imag = 0;
		}
	}
	else{
	// V = Vm .* exp(i*pi/180* Va) = Vm .* exp(1j * Va) = Vm * [cos(Va) + j*sin(Va)] = Vm*cos(va) + j*Vm*sin(Va) = Vre + j*Vimg
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,Num_of_Buses,V_Mag,V_Angl,V_Rectangular,V_Mag_Init,V_Angl_Init) if(OMP_Enbale == 1)
		for (i=0;i<Num_of_Buses;i++){
			V_Mag[i] = rec2pol_mag( V_Mag_Init[i]*cos(deg2rad(V_Angl_Init[i])), V_Mag_Init[i]*sin(deg2rad(V_Angl_Init[i])) ); // real part
			V_Angl[i] = rec2pol_ang( V_Mag_Init[i]*cos(deg2rad(V_Angl_Init[i])), V_Mag_Init[i]*sin(deg2rad(V_Angl_Init[i])) ); // imaginary part
			V_Rectangular[i].real = V_Mag[i]*cos(V_Angl[i]);
			V_Rectangular[i].imag = V_Mag[i]*sin(V_Angl[i]); 
		}
	
	// For in-service generators: V0 = [ Vgen[i] / V_Mag(i) ] * V(i)
		Chunk_Size = round(Gen_Buses_Num_Elements / OMP_Cores);
		#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i, k, j, temp, temp1, temp2) shared(Chunk_Size,Gen_Buses_Num_Elements,GenStatus,GenBusNo,PQBuses,PQ_Buses_Num_Elements,Vg,V_Mag,V_Angl,V_Rectangular) if(OMP_Enbale == 1)	
		for(i=0;i<Gen_Buses_Num_Elements;i++){
			if (GenStatus[i]==1){
				temp = 1;
				// Exclude PQ buses from generators list
				for (k = 0; k < PQ_Buses_Num_Elements; k++) {
					if (GenBusNo[i] == PQBuses[k]) {
						temp = 0;
						break;
					}
				}
				if (temp == 1) {
					j = GenBusNo[i] - 1; // generator's bus number            
					temp1 = (Vg[i] / V_Mag[j])* V_Mag[j] * cos(V_Angl[j]); //real
					temp2 = (Vg[i] / V_Mag[j])* V_Mag[j] * sin(V_Angl[j]); // imaginary
					V_Mag[j] = rec2pol_mag(temp1, temp2);
					V_Angl[j] = rec2pol_ang(temp1, temp2);
					V_Rectangular[j].real = temp1;
					V_Rectangular[j].imag = temp2;
				}
			}
		}
	} // end of else


	// Create PVPQ array
	Chunk_Size = round(PVPQsize / OMP_Cores);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(OMP_Cores) private(i) shared(Chunk_Size,PVPQsize,PVPQ,PV_Buses_Num_Elements,PQ_Buses_Num_Elements,PVBuses,PQBuses) if(OMP_Enbale == 1)
	for (i = 0; i < PVPQsize; i++) { // form PVPQ array
		if (i < PV_Buses_Num_Elements)
			PVPQ[i] = PVBuses[i];

		if (i < PQ_Buses_Num_Elements)
			PVPQ[i + PV_Buses_Num_Elements] = PQBuses[i];
	}


	// Evaluate F(x0)
	Sbus(Num_of_Buses, Gen_Buses_Num_Elements, BaseMVA, Pd, Qd, Pg, Qg, GenBusNo, GenStatus, Sbus_Real, Sbus_Imag, OMP_Enbale, OMP_Cores); // Returns vectors of complex bus power injection (Sbus_Real and Sbus_Imag)
	Power_Mimatch(Ybus_Row, Ybus_CSR, V_Rectangular, YbusV, Sbus_Real, Sbus_Imag, Num_of_Buses, Mis_Real, Mis_Imag, OMP_Enbale, OMP_Cores, Sparse_Enable);
	F_Build(Mis_Real, Mis_Imag, Num_of_Buses, PV_Buses_Num_Elements, PQ_Buses_Num_Elements, PVBuses, PQBuses, F, &F_Norm, OMP_Enbale, OMP_Cores); // Build function F and calculate its infinity norm

	if (Sparse_Enable == 0) printf("\n********************* Starting Power Flow Iterations (Dense)  *************************\n");
	if (Sparse_Enable == 1) printf("\n********************* Starting Power Flow Iterations (Sparse) *************************\n");

	printf("\nTolerance at F(x0) = %.9g\n", F_Norm); fflush(stdout);

	if (F_Norm <= Tol){
		printf("Power Flow Converged at F(x0)!\n"); fflush(stdout);
		converged = 0;
		goto memory_free;
	}
	timer_stop[2] = dsecnd(); // Initializing variables time
	timer[2] = (timer_stop[2] - timer_start[2]); // Initializing variables time

	//sparse_matrix_t test_CSR = NULL;
	

// Start Iteration
while(iter < Max_iter && F_Norm > Tol){
	iter++; // Iteration counter

	// Computes partial derivatives of power injection
	timer_start[3] = dsecnd(); // Partial derivates time
	if (Sparse_Enable == 1) {	
		CHECK_OPERATION( dS_dV_Sparse(ds_dVa_vals, ds_dVm_vals, ds_dV_IA, ds_dV_JA, Ybus_nnz, Ybus_CSR, diagIbus_CSR, diagV_CSR, diagVnorm_CSR, V_Rectangular, Ibus, V_Mag, Num_of_Buses, OMP_Enbale, OMP_Cores) );
	}
	else {
		dS_dV(Ybus_Row, V_Rectangular, Ibus, diagVnorm, YbusdiagVnorm, YbusdiagV, V_Mag, Num_of_Buses, ds_dVa, ds_dVm, OMP_Enbale, OMP_Cores);
	}
	timer_stop[3] = dsecnd(); // Partial derivates time
	timer[3] = timer[3] + (timer_stop[3] - timer_start[3]); // Partial derivates time
	
	// Form Jacobian Matrix
	if (iter == 1 && Sparse_Enable == 1) {
		timer_start[4] = dsecnd(); // Sparse Jacobian arrays time
		CHECK_OPERATION(Jacobian_Sparse_Initialize(&J_val, &J_IA, &J_JA, &J_IA_Array, &J_IA_Position, &Jnnz, Ref_Buses, Ref_Buses_Size, Ybus_nnz, ds_dVa_vals, ds_dVm_vals, ds_dV_IA, ds_dV_JA, PVPQ, Num_of_Buses, PV_Buses_Num_Elements, PQ_Buses_Num_Elements, PQBuses, OMP_Enbale, OMP_Cores));
		timer_stop[4] = dsecnd(); // Sparse Jacobian arrays time
		timer[4] = (timer_stop[4] - timer_start[4]); // Sparse Jacobian arrays time
	}
	timer_start[5] = dsecnd(); // Jacobian time
	if (Sparse_Enable == 1) {
		// Initialize Jacobian arrays at 1st iteration		
		Jacobian_Sparse(J_val, J_IA_Position, J_IA_Array, Jnnz, ds_dVa_vals, ds_dVm_vals, OMP_Enbale, OMP_Cores);
	}
	else {	
		Jacobian(ds_dVm, ds_dVa, PVPQ, Num_of_Buses, PV_Buses_Num_Elements, PQ_Buses_Num_Elements, PVBuses, PQBuses, J, OMP_Enbale, OMP_Cores);
	}
	timer_stop[5] = dsecnd(); // Jacobian time
	timer[5] = timer[5] + (timer_stop[5] - timer_start[5]); // Jacobian time


	// Solve system of linear equations
	timer_start[6] = dsecnd(); // Solver time
	/* Dense Solution */
	if (Sparse_Enable == 0) { 
		CHECK_SOLVER( LAPACKE_dgesv(LAPACK_ROW_MAJOR, Jsize, 1, J, Jsize, ipiv, F, 1) ); 
	}
	/* Sparse Solution */
	else { 
		// Call the sparse solver	
		CHECK_SOLVER( pardiso_solver(J_val, J_JA, J_IA, F, F_dummy, Jsize, OMP_Enbale, OMP_Cores) );	
	}
	timer_stop[6] = dsecnd(); // Solver time
	timer[6] = timer[6] + (timer_stop[6] - timer_start[6]); // Solver time


	// Update voltages
	timer_start[7] = dsecnd(); //Voltage updating time
	Voltage_Update(V_Mag, V_Angl, V_Rectangular, Num_of_Buses, PV_Buses_Num_Elements, PVBuses, F, PQ_Buses_Num_Elements, PQBuses, OMP_Enbale, OMP_Cores);
	timer_stop[7] = dsecnd(); //Voltage updating time
	timer[7] = timer[7] + (timer_stop[7] - timer_start[7]); //Voltage updating time

	// calculate power mismatch
	timer_start[8] = dsecnd(); // Power mismatch time
	Power_Mimatch(Ybus_Row, Ybus_CSR, V_Rectangular, YbusV, Sbus_Real, Sbus_Imag, Num_of_Buses, Mis_Real, Mis_Imag, OMP_Enbale, OMP_Cores, Sparse_Enable);
	timer_stop[8] = dsecnd(); // Power mismatch time
	timer[8] = timer[8] + (timer_stop[8] - timer_start[8]); // Power mismatch time

	// Build the norm vector
	timer_start[9] = dsecnd(); // Norm calculation time
	F_Build(Mis_Real, Mis_Imag, Num_of_Buses, PV_Buses_Num_Elements, PQ_Buses_Num_Elements, PVBuses, PQBuses, F, &F_Norm, OMP_Enbale, OMP_Cores); // Build function F and calculate its infinity norm
	timer_stop[9] = dsecnd(); // Norm calculation time
	timer[9] = timer[9] + (timer_stop[9] - timer_start[9]); // Norm calculation time
	printf("Tolerance at iteration [%i] = %.9g\n", iter, F_Norm);	fflush(stdout);

} // end of while loop

	if (F_Norm <= Tol) {
		printf("Power flow converged in [%i] iterations!\n", iter); fflush(stdout);
		converged = 0;
	}
	else{
		printf("Power flow did not converge in [%i] iterations!", iter); fflush(stdout);
		converged = 1;
	}
		

memory_free:
if (Enable_Memory_Stat == 1) {
	*Peak_Mem_Usage = mkl_peak_mem_usage(MKL_PEAK_MEM) / 1024;
	mkl_peak_mem_usage(MKL_PEAK_MEM_DISABLE);
}

printf("\n********************************* End of Power Flow ***********************************\n"); fflush(stdout);


free(Sbus_Real);
free(Sbus_Imag);
free(Mis_Real);
free(Mis_Imag);
free(F);
free(J);
free(PVPQ);
free(F_dummy);
free(J_val);

mkl_free(ipiv);
mkl_free(V_Rectangular);
mkl_free(Ibus);
mkl_free(diagVnorm);
mkl_free(YbusdiagVnorm);
mkl_free(YbusdiagV);
mkl_free(ds_dVa);
mkl_free(ds_dVm);
mkl_free(YbusV);
mkl_free(J_Sparse);
mkl_free(J_IA);
mkl_free(J_JA);
mkl_free(J_IA_Array);
mkl_free(J_IA_Position);
mkl_free(diagIbus_val);
mkl_free(diagV_val);
mkl_free(diagVnorm_val);
mkl_free(IA);
mkl_free(JA);
mkl_free(ds_dVm_vals);
mkl_free(ds_dVa_vals);
mkl_free(ds_dV_IA);
mkl_free(ds_dV_JA);

mkl_sparse_destroy(Ybus_CSR);
mkl_sparse_destroy(diagIbus_CSR);
mkl_sparse_destroy(diagV_CSR);
mkl_sparse_destroy(diagVnorm_CSR);


mkl_free_buffers();

timer_stop[0] = dsecnd(); //Total execution time
timer[0] =  (timer_stop[0] - timer_start[0]); //Total execution time
timer[10] = timer[0] / iter; // average time per iteration


return converged;

} // end of PF_NR_MKL function
