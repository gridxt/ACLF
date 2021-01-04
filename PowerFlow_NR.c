#include "PowerFlow_NR.h"
int Chunk_Size;
int Thread_Num = 1;


// ********************************* Sbus ***************************
// This function calculates complex bus power injection
void Sbus_NR(int NumberOfBuses, int Gen_Num_Buses, double BaseMVA, double *Pd, double *Qd, double *Pg, double *Qg, int *GenBusNo, int *GenStatus, double **Sbus_Real, double **Sbus_Imag){
int i = 0;
int Bus_No = 0;

// Negative of Demand
Chunk_Size = round(NumberOfBuses / Thread_Num);
#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(Thread_Num) private(i) shared(Chunk_Size, NumberOfBuses, Sbus_Real, Sbus_Imag, Pd, Qd, BaseMVA)
for (i=0;i<NumberOfBuses;i++){
    (*Sbus_Real)[i] = (-1)*Pd[i]/BaseMVA;
    (*Sbus_Imag)[i] = (-1)*Qd[i]/BaseMVA;
}

//Online Generation - Demand
Chunk_Size = round(Gen_Num_Buses / Thread_Num);
#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(Thread_Num) private(i, Bus_No) shared(Chunk_Size, Gen_Num_Buses, GenBusNo, Sbus_Real, Sbus_Imag, Pg, Qg, GenStatus, BaseMVA)
for (i=0;i<Gen_Num_Buses;i++){
    Bus_No = GenBusNo[i] - 1;
    (*Sbus_Real)[Bus_No] = (Pg[i]*GenStatus[i]/BaseMVA)+(*Sbus_Real)[Bus_No];
    (*Sbus_Imag)[Bus_No] = (Qg[i]*GenStatus[i]/BaseMVA)+(*Sbus_Imag)[Bus_No];
}

}

// ********************************** Power Mismatch ***************
// This function calculates the power mismatch
void Power_Mimatch_NR(MKL_Complex16 *Ybus_Row_, double *V_mag, double *V_ang, double *Sbus_re, double *Sbus_img, int Num_Buses, double **Mis_re, double **Mis_img){
int i = 0;
int j = 0;
double temp;
// Power Mismatch = V .* conj(Ybus * V) - Sbus(Vm)
// Ybus mxm
// V mx1
// Sbus mx1
// Mis_Real and Mis_Imag mx1

Chunk_Size = round(Num_Buses / Thread_Num);
#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(Thread_Num) private(i, j, temp) shared(Chunk_Size, Num_Buses, Mis_re, Mis_img, Ybus_Row_, V_mag, V_ang, Sbus_re, Sbus_img) 
for (i=0;i<Num_Buses;i++){
(*Mis_re)[i] = 0;
(*Mis_img)[i] = 0;

	for (j=0;j<Num_Buses;j++){			
		//Ybus * V       		
		(*Mis_re)[i] += Mult_re(Ybus_Row_[i*Num_Buses + j].real, Ybus_Row_[i*Num_Buses + j].imag, V_mag[j] * cos(V_ang[j]), V_mag[j] * sin(V_ang[j]));
		(*Mis_img)[i] += Mult_img(Ybus_Row_[i*Num_Buses + j].real, Ybus_Row_[i*Num_Buses + j].imag, V_mag[j]*cos(V_ang[j]), V_mag[j]*sin(V_ang[j]));		
	}

// conj(Ybus * V)
(*Mis_img)[i] = (*Mis_img)[i] * (-1);

// V .* conj(Ybus * V) - Sbus
temp = (*Mis_re)[i]; // save the value of Mis_re before we update it
(*Mis_re)[i] = Mult_re(V_mag[i] * cos(V_ang[i]), V_mag[i] * sin(V_ang[i]), (*Mis_re)[i], (*Mis_img)[i]) - Sbus_re[i];
(*Mis_img)[i] = Mult_img(V_mag[i] * cos(V_ang[i]), V_mag[i] * sin(V_ang[i]), temp, (*Mis_img)[i]) - Sbus_img[i];

}

}

// ***************************** Objective Function ****************
// This function builds norm vector [ Mis_re[PV;PQ] ; Mis_img[PQ] ] and calculate the infinity norm (largest number in the vector)
void F_Build_NR(double *Mis_re, double *Mis_img, int Num_Buses, int PV_Elem, int PQ_Elem, int *PVBuses, int *PQBuses, int Ref_Bus_No, double **F_, double *F_Norm){
int i = 0;
int F_index = 0;
(*F_Norm) = 0;

Chunk_Size = round(PV_Elem / Thread_Num);
#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(Thread_Num) private(i) shared(Chunk_Size, PV_Elem, F_, Mis_re,  PVBuses, F_Norm) 
for (i=0;i<PV_Elem;i++){
    (*F_)[i] = Mis_re[ PVBuses[i] - 1 ] ;	
    if ( *F_Norm < fabs((*F_)[i]) ) *F_Norm = fabs((*F_)[i]); // find the largest number
}

//Chunk_Size = round(PQ_Elem / Thread_Num);
//#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(Thread_Num) private(i, F_index) shared(Chunk_Size, PQ_Elem, PV_Elem, F_, Mis_re, F_Norm, PQBuses, Mis_img) 
for (i=0;i<PQ_Elem;i++){
    F_index = i + PV_Elem ;
    (*F_)[F_index] = Mis_re[ PQBuses[i] - 1 ] ;	
    if ( *F_Norm < fabs((*F_)[F_index]) ) *F_Norm = fabs((*F_)[F_index]);

    F_index = i + PV_Elem + PQ_Elem;
    (*F_)[F_index] = Mis_img[ PQBuses[i] - 1 ] ;	
    if ( *F_Norm < fabs((*F_)[F_index]) ) *F_Norm = fabs((*F_)[F_index]);

}

}

// ***************************** dSbus_dV ***************************
//  Computes partial derivatives of power injection with respect to voltage.
// http://www.pserc.cornell.edu/matpower/TN2-OPF-Derivatives.pdf
// Calculate Ibus (m x 1) = Ybus (m x m) * V (m x 1)
// Calculate dS_dVm_re = diagV * conj(Ybus * diagVnorm) + conj(diagIbus) * diagVnorm
// Calculate dSbus_dVa = conj(diagIbus - Ybus * diagV);

void dS_dV_NR(MKL_Complex16 *Ybus_Row_, double *V_mag, double *V_ang, int Num_Buses, double **dS_dVm_re, double **dS_dVm_img, double **dS_dVa_re, double **dS_dVa_img){
int i = 0;
int j = 0;
double temp1 = 0;
double temp2 = 0;
double *Ibus_re = malloc(Num_Buses * sizeof(double));
double *Ibus_img = malloc(Num_Buses * sizeof(double));
double *VNorm_mag = malloc(Num_Buses * sizeof(double)); // (V_mag ./ abs(V_mag)) = 1
double *VNorm_ang = malloc(Num_Buses * sizeof(double));
double *V_rev_diag_re = malloc(Num_Buses * sizeof(double)); // real(1j * V)
double *V_rev_diag_img = malloc(Num_Buses * sizeof(double)); // imaginary(1j * V)


MKL_Complex16 alpha, beta;
alpha.real = 1;
alpha.imag = 0;
beta.real = 0;
beta.imag = 0;

//MKL_Complex16 *Dummy1, *Dummy2;
//Dummy1 = (MKL_Complex16 *)calloc(Num_Buses*Num_Buses, sizeof(MKL_Complex16));
//Dummy2 = (MKL_Complex16 *)calloc(Num_Buses*Num_Buses, sizeof(MKL_Complex16));


// Calculate norm of V : V./abs(V) & 1j * V
Chunk_Size = round(Num_Buses / Thread_Num);
#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(Thread_Num) private(i, temp1, temp2) shared(Chunk_Size, Num_Buses, VNorm_mag, VNorm_ang,  V_rev_diag_re, V_rev_diag_img, V_mag, V_ang)
for (i=0; i<Num_Buses; i++){
temp1 = ( V_mag[i]*cos(V_ang[i]) )/V_mag[i]; // real V/mag V
temp2 = ( V_mag[i]*sin(V_ang[i]) )/V_mag[i]; // imag V/mag v
// (a+bj)/abs(a+bj) = (a/abs(a+bj)) + j(b/abs(a+bj))
VNorm_mag[i] = 1; // it is equal to one all the times!
VNorm_ang[i] = rec2pol_ang(temp1, temp2);

//1j * diagV  --- 1j*diaV = [digaVimg*(-1) + digaVre*(i)] : same magnitude but different angle
V_rev_diag_re[i] = V_mag[i]*sin(V_ang[i])*(-1);
V_rev_diag_img[i] = V_mag[i]*cos(V_ang[i]);
}


Chunk_Size = round(Num_Buses / Thread_Num);
#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(Thread_Num) private(i, j, temp1, temp2) shared(Chunk_Size, Num_Buses, Ybus_Row_, Ibus_re, Ibus_img, VNorm_mag, VNorm_ang, V_mag, V_ang, dS_dVm_re, dS_dVm_img, dS_dVa_re, dS_dVa_img)
for (i=0; i<Num_Buses; i++){
Ibus_re[i] = 0;
Ibus_img[i] = 0;	
	for (j=0; j<Num_Buses; j++){
        // Ibus = Ybus * V
        Ibus_re[i] = Ibus_re[i] + Mult_re( Ybus_Row_[i*Num_Buses + j].real, Ybus_Row_[i*Num_Buses + j].imag, V_mag[j]*cos(V_ang[j]), V_mag[j]*sin(V_ang[j]) ); // real
        Ibus_img[i] = Ibus_img[i] + Mult_img(Ybus_Row_[i*Num_Buses + j].real, Ybus_Row_[i*Num_Buses + j].imag, V_mag[j]*cos(V_ang[j]), V_mag[j]*sin(V_ang[j]) ); // imaginary

        // diagV * conj(Ybus * diagVnorm)
        temp1 = Mult_re(Ybus_Row_[i*Num_Buses + j].real, Ybus_Row_[i*Num_Buses + j].imag, VNorm_mag[j]*cos(VNorm_ang[j]), VNorm_mag[j]*sin(VNorm_ang[j]) ); // real( Ybus * diagVnorm )
        temp2 = Mult_img(Ybus_Row_[i*Num_Buses + j].real, Ybus_Row_[i*Num_Buses + j].imag, VNorm_mag[j]*cos(VNorm_ang[j]), VNorm_mag[j]*sin(VNorm_ang[j]) ) * (-1); // conjugate of imag( Ybus * diagVnorm )
        dS_dVm_re[i][j] = Mult_re( V_mag[i]*cos(V_ang[i]), V_mag[i]*sin(V_ang[i]), temp1, temp2 ); // real
        dS_dVm_img[i][j] = Mult_img( V_mag[i]*cos(V_ang[i]), V_mag[i]*sin(V_ang[i]), temp1, temp2 ); //imaginary

        // Ybus * digaV
        dS_dVa_re[i][j] = Mult_re(Ybus_Row_[i*Num_Buses + j].real, Ybus_Row_[i*Num_Buses + j].imag, V_mag[j]*cos(V_ang[j]), V_mag[j]*sin(V_ang[j]) ); // real
        dS_dVa_img[i][j] = Mult_img(Ybus_Row_[i*Num_Buses + j].real, Ybus_Row_[i*Num_Buses + j].imag, V_mag[j]*cos(V_ang[j]), V_mag[j]*sin(V_ang[j]) ); // imag
        if (i!=j) dS_dVa_re[i][j] = dS_dVa_re[i][j] * (-1); // to compensate for zero elements in diagIbus when we calculate diagIbus - Ybus * diagV
    }

 // dS_dVm_re = diagV * conj(Ybus * diagVnorm) + conj(diagIbus) * diagVnorm
 dS_dVm_re[i][i] = dS_dVm_re[i][i] + Mult_re( Ibus_re[i], Ibus_img[i]*(-1) , VNorm_mag[i]*cos(VNorm_ang[i]), VNorm_mag[i]*sin(VNorm_ang[i]) );
 dS_dVm_img[i][i] = dS_dVm_img[i][i] + Mult_img( Ibus_re[i], Ibus_img[i]*(-1) , VNorm_mag[i]*cos(VNorm_ang[i]), VNorm_mag[i]*sin(VNorm_ang[i]) );

// diagonal elements = conj(Ibus - Ybus * diagV)
 dS_dVa_re[i][i] = Ibus_re[i] - dS_dVa_re[i][i];
 dS_dVa_img[i][i] = (Ibus_img[i] - dS_dVa_img[i][i]) * (-1);
}


//1j * diagV * conj(diagIbus - Ybus * diagV)
Chunk_Size = round(Num_Buses / Thread_Num);
#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(Thread_Num) collapse(2) private(i, j, temp1, temp2) shared(Chunk_Size, Num_Buses, V_rev_diag_re, V_rev_diag_img, dS_dVa_re, dS_dVa_img)
for (i=0;i<Num_Buses;i++){	
	for (j=0;j<Num_Buses;j++){
        temp1 = Mult_re( V_rev_diag_re[i], V_rev_diag_img[i], dS_dVa_re[i][j], dS_dVa_img[i][j]);
        temp2 = Mult_img( V_rev_diag_re[i], V_rev_diag_img[i], dS_dVa_re[i][j], dS_dVa_img[i][j]);
        dS_dVa_re[i][j] = temp1;
        dS_dVa_img[i][j] = temp2;
     }
}

free(V_rev_diag_re);
free(V_rev_diag_img);
free(Ibus_re);
free(Ibus_img);
free(VNorm_mag);
free(VNorm_ang);
}

// ********************************* Jacobian Matrix *********************
// This function forms the Jacobian matrix
void Jacobian_NR(double **dS_dVm_re, double **dS_dVm_img, double **dS_dVa_re, double **dS_dVa_img, int Num_Buses, int PV_Elem, int PQ_Elem, int *PVBuses, int *PQBuses, int Ref_Bus_No, double *J_){
int i = 0;
int j = 0;
int PVPQsize = PV_Elem + PQ_Elem;
int *PVPQ = malloc( PVPQsize * sizeof(int)); // contains all PV and PQ bus numbers [PV;PQ]
int Jsize = PV_Elem + (PQ_Elem * 2);

Chunk_Size = round(PVPQsize / Thread_Num);
#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(Thread_Num) private(i) shared(Chunk_Size, PVPQsize, PV_Elem, PQ_Elem, PVPQ, PVBuses, PQBuses)
for (i=0;i<PVPQsize;i++){ // form PVPQ array

    if (i<PV_Elem)
        PVPQ[i] = PVBuses[i];

    if (i<PQ_Elem)
        PVPQ[i+PV_Elem] = PQBuses[i];
}


Chunk_Size = round(PVPQsize / Thread_Num);
#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(Thread_Num) collapse(2) private(i, j) shared(Chunk_Size, PVPQsize, J_, Jsize, dS_dVa_re, PVPQ, PQ_Elem, PQBuses, dS_dVa_img, dS_dVm_re, dS_dVm_img)
for (i=0; i<PVPQsize; i++){ // J = [J1 J2; J3 J4]
    for (j=0;j<PVPQsize;j++){
      
		J_[(i*Jsize) + j] = dS_dVa_re[PVPQ[i] - 1][PVPQ[j] - 1];

		if (j < PQ_Elem) 
			J_[(i*Jsize) + (j + PVPQsize)] = dS_dVm_re[PVPQ[i] - 1][PQBuses[j] - 1];           

		if (i < PQ_Elem) 			
			J_[(i+PVPQsize)*Jsize + j] = dS_dVa_img[PQBuses[i] - 1][PVPQ[j] - 1];            

		if ((i < PQ_Elem) && (j < PQ_Elem)) 			
			J_[(i + PVPQsize)*Jsize + (j + PVPQsize)] = dS_dVm_img[PQBuses[i] - 1][PQBuses[j] - 1];		
            
    }
}

free(PVPQ);

}


// ***********************************************************************************************************************************************
// ************************************************************** Main Function ******************************************************************

// This function solves the power flow using full newton method
int PF_NR(MKL_Complex16 *Ybus_Row, double *V_Mag_, double *V_Angl_, int Num_of_Buses,
           int *PVBuses, int *PQBuses,
           int *GenBusNo, int *GenStatus, double BaseMVA, double *Pd, double *Qd, double *Pg, double *Qg, double *Vg,
           int OnlineGen_Num_Elements, int PV_Buses_Num_Elements, int PQ_Buses_Num_Elements, int Ref_Bus,
           int FlatStart, int Max_iter, double Tol, long double **CPU_Execution_TIME){

clock_t start, stop, start1, stop1, start2,  stop2;




int i = 0;
int j = 0;
int iter = 0; // iteration counter
double temp1 = 0; // dummy
double temp2 = 0; // dummy
MKL_INT info; // check singularity of Jacobian

double *V_Mag = malloc(Num_of_Buses * sizeof(double)); // initial V - Magnitude
double *V_Angl = malloc(Num_of_Buses * sizeof(double)); // initial V - Angle
double *Sbus_Real = malloc(Num_of_Buses * sizeof(double)); // real part of complex bus power injection (P)
double *Sbus_Imag = malloc(Num_of_Buses * sizeof(double)); // imaginary part of complex bus power injection (Q)
double *Mis_Real = malloc(Num_of_Buses * sizeof(double)); // real part of power mismatch vector
double *Mis_Imag = malloc(Num_of_Buses * sizeof(double)); // imaginary part of power mismatch vector
double *F = malloc( (PV_Buses_Num_Elements + PQ_Buses_Num_Elements + PQ_Buses_Num_Elements) * sizeof(double)); // update state: imaginary part of power mismatch vector
double F_Norm = 0; // contains the infinity norm of vector F
int Jsize = PV_Buses_Num_Elements + (PQ_Buses_Num_Elements*2); // size of Jacobian matrix
double *J = (double *)calloc(Jsize*Jsize, sizeof(double)); //Jacobian Matrix
MKL_INT *ipiv = malloc(Jsize * sizeof(MKL_INT));

double ** dS_dVm_re = (double **)calloc(Num_of_Buses, sizeof(double *));  // partial derivatives of the complex bus power injections w.r.t voltage magnitude (for all buses)
dS_dVm_re[0] = (double *)calloc(Num_of_Buses*Num_of_Buses, sizeof(double));
    for(i = 0; i < Num_of_Buses; i++)
        dS_dVm_re[i] = (*dS_dVm_re + Num_of_Buses * i);

double ** dS_dVm_img = (double **)calloc(Num_of_Buses, sizeof(double *));
dS_dVm_img[0] = (double *)calloc(Num_of_Buses*Num_of_Buses, sizeof(double));
    for(i = 0; i < Num_of_Buses; i++)
        dS_dVm_img[i] = (*dS_dVm_img + Num_of_Buses * i);

double ** dS_dVa_re = (double **)calloc(Num_of_Buses, sizeof(double *));  // partial derivatives of the complex bus power injections w.r.t voltage angle respectively (for all buses)
dS_dVa_re[0] = (double *)calloc(Num_of_Buses*Num_of_Buses, sizeof(double));
    for(i = 0; i < Num_of_Buses; i++)
        dS_dVa_re[i] = (*dS_dVa_re + Num_of_Buses * i);

double ** dS_dVa_img = (double **)calloc(Num_of_Buses, sizeof(double *));
dS_dVa_img[0] = (double *)calloc(Num_of_Buses*Num_of_Buses, sizeof(double));
    for(i = 0; i < Num_of_Buses; i++)
        dS_dVa_img[i] = (*dS_dVa_img + Num_of_Buses * i);

double *Jinv = (double *)calloc(Jsize*Jsize, sizeof(double *)); // Inverse of Jacobian Matrix
double *dx = (double *)calloc(Jsize, sizeof(double)); // update step
start = clock();



start2 = clock();
// Initialize V
if (FlatStart == 1){
//  V0 = 1 < 0
	Chunk_Size = round(Num_of_Buses / Thread_Num);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(Thread_Num) private(i) shared(Chunk_Size, Num_of_Buses, V_Mag, V_Angl)
	for (i=0; i<Num_of_Buses; i++){
        V_Mag[i] = 1;
        V_Angl[i] = 0;
    }
}
else{
// V = Vm .* exp(i*pi/180* Va) = Vm .* exp(1j * Va) = Vm * [cos(Va) + j*sin(Va)] = Vm*cos(va) + j*Vm*sin(Va) = Vre + j*Vimg
	Chunk_Size = round(Num_of_Buses / Thread_Num);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(Thread_Num) private(i) shared(Chunk_Size, Num_of_Buses, V_Mag, V_Angl, V_Mag_, V_Angl_)
	for (i=0;i<Num_of_Buses;i++){
        V_Mag[i] = rec2pol_mag( V_Mag_[i]*cos(deg2rad(V_Angl_[i])), V_Mag_[i]*sin(deg2rad(V_Angl_[i])) ); // real part
        V_Angl[i] = rec2pol_ang( V_Mag_[i]*cos(deg2rad(V_Angl_[i])), V_Mag_[i]*sin(deg2rad(V_Angl_[i])) ); // imaginary part
    }
// For in-service generators: V0 = [ Vgen[i] / V_Mag(i) ] * V(i)
	Chunk_Size = round(OnlineGen_Num_Elements / Thread_Num);
	#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(Thread_Num) private(i, j, temp1, temp2) shared(Chunk_Size, OnlineGen_Num_Elements, GenStatus, GenBusNo, V_Mag, V_Angl, Vg)
	for(i=0;i<OnlineGen_Num_Elements;i++){
        if (GenStatus[i]==1){
            j = GenBusNo[i] - 1; // generator's bus number            
            temp1 = ( Vg[i]/V_Mag[j] )* V_Mag[j]*cos(V_Angl[j]); //real
            temp2 = ( Vg[i]/V_Mag[j] )* V_Mag[j]*sin(V_Angl[j]); // imaginary
            V_Mag[j] = rec2pol_mag(temp1, temp2);
            V_Angl[j] = rec2pol_ang(temp1, temp2);
        }
    }
} // end of else


// Evaluate F(x0)
Sbus_NR(Num_of_Buses, OnlineGen_Num_Elements, BaseMVA, Pd, Qd, Pg, Qg, GenBusNo, GenStatus, &Sbus_Real, &Sbus_Imag); // Returns vectors of complex bus power injection (Sbus_Real and Sbus_Imag)
Power_Mimatch_NR(Ybus_Row, V_Mag, V_Angl, Sbus_Real, Sbus_Imag, Num_of_Buses, &Mis_Real, &Mis_Imag); // Calculate power mismatch
F_Build_NR(Mis_Real, Mis_Imag, Num_of_Buses, PV_Buses_Num_Elements, PQ_Buses_Num_Elements, PVBuses, PQBuses, Ref_Bus, &F, &F_Norm); // Build function F and calculate its infinity norm
printf("\n********************* Starting Power Flow Iterations *************************\n");
printf("\nTolerance at F(x0) = %.9g\n", F_Norm);
fflush(stdout);

if (F_Norm <= Tol){
    printf("Power Flow Converged at F(x0)!\n");
    exit(1);
}
stop2 = clock();
(*CPU_Execution_TIME)[4] = (((long double)(stop2 - start2)) / CLOCKS_PER_SEC);

// Start Iteration
while(iter < Max_iter && F_Norm >= Tol){
iter++; // Iteration counter

start2 = clock();
dS_dV_NR(Ybus_Row, V_Mag, V_Angl, Num_of_Buses, dS_dVm_re, dS_dVm_img, dS_dVa_re, dS_dVa_img); // Computes partial derivatives of power injection
stop2 = clock();
(*CPU_Execution_TIME)[5] = (*CPU_Execution_TIME)[5] + (((long double)(stop2 - start2)) / CLOCKS_PER_SEC);

start2 = clock();
Jacobian_NR( dS_dVm_re, dS_dVm_img, dS_dVa_re, dS_dVa_img, Num_of_Buses, PV_Buses_Num_Elements, PQ_Buses_Num_Elements, PVBuses, PQBuses, Ref_Bus, J); // Compute Jacobian Matrix
stop2 = clock();
(*CPU_Execution_TIME)[6] = (*CPU_Execution_TIME)[6] + (((long double)(stop2 - start2)) / CLOCKS_PER_SEC);

start1 = clock(); // // calculation time of LU and INV

LU(Jsize, J, Thread_Num); // Calculate the LU of J and store the results in the same matrix J
INVS(Jsize, Jinv, J, Thread_Num); // Calculate the inverse of J and save it in Jinv
// Compute the update state dx using dx = -(J \ F) = - (Jinvs * F)
for (i = 0; i<Jsize; i++) {
	dx[i] = 0;
	for (j = 0; j<Jsize; j++) {
		dx[i] = dx[i] + Jinv[i*Jsize+j] * F[j];		
		//printf("dx[%i] = %f \n", i, dx[i]);
	}
	
}


/*
info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, Jsize, 1, J, Jsize, ipiv, F, 1); // Calculate Jacobian Inverse

if (info > 0) {
		printf("The diagonal element of the triangular factor of A,\n");
		printf("U(%i,%i) is zero, so that A is singular;\n", info, info);
		printf("the solution could not be computed.\n");
		exit(0);
}
*/

stop1 = clock();
(*CPU_Execution_TIME)[7] = (*CPU_Execution_TIME)[7] + ( ((long double) (stop1 - start1)) / CLOCKS_PER_SEC );

start2 = clock();


// Update voltages
Chunk_Size = round(PV_Buses_Num_Elements / Thread_Num);
#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(Thread_Num) private(i, j) shared(Chunk_Size, PV_Buses_Num_Elements, PVBuses, V_Angl, dx)
for (i = 0; i<PV_Buses_Num_Elements; i++) {
	j = PVBuses[i] - 1; // PV bus number
	V_Angl[j] = V_Angl[j] + dx[i]*(-1); // updated angle
}

Chunk_Size = round(PQ_Buses_Num_Elements / Thread_Num);
#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(Thread_Num) private(i, j) shared(Chunk_Size, PQ_Buses_Num_Elements, PQBuses, PV_Buses_Num_Elements, V_Mag, V_Angl, dx)
for (i = 0; i<PQ_Buses_Num_Elements; i++) {
	j = PQBuses[i] - 1; // PQ bus number
	V_Angl[j] = V_Angl[j] + dx[i + PV_Buses_Num_Elements]*(-1); // update angle
	V_Mag[j] = V_Mag[j] + dx[i + PV_Buses_Num_Elements + PQ_Buses_Num_Elements]*(-1); // update magnitude
}

// Update Vm and Va again in case we wrapped around with a negative Vm : V = Vm .* exp(1j * Va) = Vm * [cos(Va) + j*sin(Va)] = Vm*cos(va) + j*Vm*sin(Va) = Vre + j*Vimg
Chunk_Size = round(Num_of_Buses / Thread_Num);
#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(Thread_Num) private(i, temp1, temp2) shared(Chunk_Size, Num_of_Buses, V_Mag, V_Angl)
for (i = 0; i<Num_of_Buses; i++) {
	temp1 = V_Mag[i] * cos(V_Angl[i]); // real part
	temp2 = V_Mag[i] * sin(V_Angl[i]); // imaginary part
	V_Mag[i] = rec2pol_mag(temp1, temp2);
	V_Angl[i] = rec2pol_ang(temp1, temp2);	
}



/*
Chunk_Size = round(PV_Buses_Num_Elements / Thread_Num);
#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(Thread_Num) private(i, j) shared(Chunk_Size, PV_Buses_Num_Elements, PVBuses, V_Angl, F)
for (i=0;i<PV_Buses_Num_Elements;i++){
    j = PVBuses[i] - 1; // PV bus number
    V_Angl[j] = V_Angl[j] + F[i]*(-1); // updated angle
}

Chunk_Size = round(PQ_Buses_Num_Elements / Thread_Num);
#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(Thread_Num) private(i, j) shared(Chunk_Size, PQ_Buses_Num_Elements, PQBuses, PV_Buses_Num_Elements, V_Mag, V_Angl, F)
for (i=0;i<PQ_Buses_Num_Elements;i++){
    j= PQBuses[i] - 1; // PQ bus number
    V_Angl[j] = V_Angl[j] + F[i+PV_Buses_Num_Elements]*(-1); // update angle
    V_Mag[j] = V_Mag[j] + F[i+PV_Buses_Num_Elements+PQ_Buses_Num_Elements]*(-1); // update magnitude
}

// Update Vm and Va again in case we wrapped around with a negative Vm : V = Vm .* exp(1j * Va) = Vm * [cos(Va) + j*sin(Va)] = Vm*cos(va) + j*Vm*sin(Va) = Vre + j*Vimg
Chunk_Size = round(Num_of_Buses / Thread_Num);
#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(Thread_Num) private(i, temp1, temp2) shared(Chunk_Size, Num_of_Buses, V_Mag, V_Angl)
for (i=0;i<Num_of_Buses;i++){
    temp1 = V_Mag[i]*cos(V_Angl[i]); // real part
    temp2 = V_Mag[i]*sin(V_Angl[i]); // imaginary part
    V_Mag[i] = rec2pol_mag(temp1, temp2);
    V_Angl[i] = rec2pol_ang(temp1, temp2);
}
*/
stop2 = clock();
(*CPU_Execution_TIME)[8] = (*CPU_Execution_TIME)[8] + (((long double)(stop2 - start2)) / CLOCKS_PER_SEC);

start2 = clock();
Power_Mimatch_NR(Ybus_Row, V_Mag, V_Angl, Sbus_Real, Sbus_Imag, Num_of_Buses, &Mis_Real, &Mis_Imag); // Calculate power mismatch
stop2 = clock();
(*CPU_Execution_TIME)[9] = (*CPU_Execution_TIME)[9] + (((long double)(stop2 - start2)) / CLOCKS_PER_SEC);

start2 = clock();
F_Build_NR(Mis_Real, Mis_Imag, Num_of_Buses, PV_Buses_Num_Elements, PQ_Buses_Num_Elements, PVBuses, PQBuses, Ref_Bus, &F, &F_Norm); // Build function F and calculate its infinity norm
stop2 = clock();
(*CPU_Execution_TIME)[10] = (*CPU_Execution_TIME)[10] + (((long double)(stop2 - start2)) / CLOCKS_PER_SEC);

printf("Tolerance at iteration [%i] = %.9g\n", iter, F_Norm);	fflush(stdout);

} // end of while loop


stop = clock();
(*CPU_Execution_TIME)[3] = ((long double)(stop - start)) / CLOCKS_PER_SEC;


if (F_Norm <= Tol)
    printf("Power flow converged after [%i] iterations!\n", iter);

if (F_Norm > Tol)
    printf("Power flow did not converge after [%i] iterations!", iter);

printf("\n***************************** End of Power Flow ******************************\n");
fflush(stdout);

free(V_Mag);
free(V_Angl);
free(dS_dVm_re[0]);
free(dS_dVm_re);
free(dS_dVa_re[0]);
free(dS_dVa_re);
free(J);


return 1;

} // end of function
