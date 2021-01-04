#include "YbusBuild.h"


int MakeYbus(double BaseMVA, int Num_of_Buses, int Num_of_Branches, int *FromBus, int *ToBus, int *BranchStatus, double *r, double *x, double *b, double *TapRatio, double *ShiftAngle, double *Gs, double *Bs, MKL_Complex16 *Ybus, double **CPU_Execution_TIME)
{

clock_t start, stop;
start = clock();

long int i = 0;
long int j = 0;
long int k = 0;
long int Fbus = 0;
long int Tbus = 0;
double temp1 = 0;
double temp2 = 0;

double *Ys_Re = malloc(Num_of_Branches * sizeof(double));
double *Ys_Img = malloc(Num_of_Branches * sizeof(double));
double *Tap_Re = malloc(Num_of_Branches * sizeof(double));
double *Tap_Img = malloc(Num_of_Branches * sizeof(double));
double *Ytt_Re = malloc(Num_of_Branches * sizeof(double));
double *Ytt_Img = malloc(Num_of_Branches * sizeof(double));
double *Yff_Re = malloc(Num_of_Branches * sizeof(double));
double *Yff_Img = malloc(Num_of_Branches * sizeof(double));
double *Yft_Re = malloc(Num_of_Branches * sizeof(double));
double *Yft_Img = malloc(Num_of_Branches * sizeof(double));
double *Ytf_Re = malloc(Num_of_Branches * sizeof(double));
double *Ytf_Img = malloc(Num_of_Branches * sizeof(double));


MKL_Complex16 alpha, beta;
alpha.real = 1;
alpha.imag = 0;
beta.real = 0;
beta.imag = 0;


MKL_Complex16 *Yf, *Yt;
Yf = (MKL_Complex16 *)calloc(Num_of_Branches*Num_of_Buses , sizeof(MKL_Complex16)); // Yf = [ Yff Yft]
Yt = (MKL_Complex16 *)calloc(Num_of_Branches*Num_of_Buses , sizeof(MKL_Complex16)); // Yt = [Ytf Ytt]

MKL_Complex16 *CfPrime, *CtPrime;
CfPrime = (MKL_Complex16 *)calloc(Num_of_Buses*Num_of_Branches, sizeof(MKL_Complex16)); // connection matrix for line & from buses (CfPrime actually contains transposed connection matrix to improve computation performance)
CtPrime = (MKL_Complex16 *)calloc(Num_of_Buses*Num_of_Branches, sizeof(MKL_Complex16)); // connection matrix for line & to buses (CtPrime actually contains transposed connection matrix to improve computation performance)

MKL_Complex16 *Dummy1, *Dummy2;
Dummy1 = (MKL_Complex16 *)calloc(Num_of_Buses*Num_of_Buses, sizeof(MKL_Complex16)); 
Dummy2 = (MKL_Complex16 *)calloc(Num_of_Buses*Num_of_Buses, sizeof(MKL_Complex16));


/*
%% for each branch, compute the elements of the branch admittance matrix where
%%
%%      | If |   | Yff  Yft |   | Vf |
%%      |    | = |          | * |    |
%%      | It |   | Ytf  Ytt |   | Vt |
%%
*/


// Initialize data
for (i=0;i<Num_of_Branches;++i){

    if (BranchStatus[i] == 1){ // Series admittance: BranchStatus[i]/(r[i]+x[i]*I)
        temp1 = 1/rec2pol_mag(r[i], x[i]);
        temp2 = (-1)*rec2pol_ang(r[i], x[i]);        
		Ys_Re[i] = temp1 * cos(temp2);
        Ys_Img[i] = temp1 * sin(temp2);
    }
    else{
        Ys_Re[i] = 0;
        Ys_Img[i] = 0;
    }

    Tap_Re[i] = TapRatio[i] * cos(deg2rad(ShiftAngle[i])); // Add phase shifter to tap ratio: TapRatio[i] * cexpl(deg2rad(ShiftAngle[i])*I)
    Tap_Img[i] = TapRatio[i] * sin(deg2rad(ShiftAngle[i]));

    Ytt_Re[i] = Ys_Re[i]; //Ys[i] + (Bc[i]/2)*I
    Ytt_Img[i] = Ys_Img[i] + (b[i]/2) * BranchStatus[i];

    Yff_Re[i] = Ytt_Re[i] / ((Tap_Re[i]*Tap_Re[i]) + (Tap_Img[i]*Tap_Img[i])); // Ytt[i] / (Tap[i] * conjl(Tap[i]))
    Yff_Img[i] = Ytt_Img[i] / ((Tap_Re[i]*Tap_Re[i]) + (Tap_Img[i]*Tap_Img[i]));

    temp1 = rec2pol_mag((-1)*Ys_Re[i],(-1)*Ys_Img[i])/rec2pol_mag(Tap_Re[i],(-1)*Tap_Img[i]);
    temp2 = rec2pol_ang((-1)*Ys_Re[i],(-1)*Ys_Img[i]) - rec2pol_ang(Tap_Re[i],(-1)*Tap_Img[i]);
    Yft_Re[i] = temp1 * cos(temp2); //(Ys[i]*(-1)) / conjl(Tap[i])
    Yft_Img[i] = temp1 * sin(temp2);

    temp1 = rec2pol_mag((-1)*Ys_Re[i],(-1)*Ys_Img[i])/rec2pol_mag(Tap_Re[i], Tap_Img[i]);
    temp2 = rec2pol_ang((-1)*Ys_Re[i],(-1)*Ys_Img[i]) - rec2pol_ang(Tap_Re[i], Tap_Img[i]);
    Ytf_Re[i] = temp1 * cos(temp2); //(Ys[i]*(-1)) / conjl(Tap[i])
    Ytf_Img[i] = temp1 * sin(temp2);

    // Form CfPrime, CtPrime, Yf, and Yt
    Fbus = FromBus[i] - 1;
    Tbus = ToBus[i] - 1;
    
	CfPrime[Fbus*Num_of_Branches + i].real = 1;
	CfPrime[Fbus*Num_of_Branches + i].imag = 0;
	CtPrime[Tbus*Num_of_Branches + i].real = 1;
	CtPrime[Tbus*Num_of_Branches + i].imag = 0;
	
	Yf[i*Num_of_Buses + Fbus].real = Yff_Re[i];
	Yf[i*Num_of_Buses + Fbus].imag = Yff_Img[i];
	Yf[i*Num_of_Buses + Tbus].real = Yft_Re[i];
	Yf[i*Num_of_Buses + Tbus].imag = Yft_Img[i];

	Yt[i*Num_of_Buses + Fbus].real = Ytf_Re[i];
	Yt[i*Num_of_Buses + Fbus].imag = Ytf_Img[i];
	Yt[i*Num_of_Buses + Tbus].real = Ytt_Re[i];
	Yt[i*Num_of_Buses + Tbus].imag = Ytt_Img[i];
	
}


// Calculate Ybus = Cf' * Yf + Ct' * Yt  +.... ::Branch admittances
//                  + Ysh :: Shunt admittance
cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Num_of_Buses, Num_of_Buses, Num_of_Branches, &alpha, CfPrime, Num_of_Branches, Yf, Num_of_Buses, &beta, Dummy1, Num_of_Buses);
cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Num_of_Buses, Num_of_Buses, Num_of_Branches, &alpha, CtPrime, Num_of_Branches, Yt, Num_of_Buses, &beta, Dummy2, Num_of_Buses);
vzAdd(Num_of_Buses*Num_of_Buses, Dummy1, Dummy2, Ybus);


// add Ysh to diagonal elements
for(i = 0; i < Num_of_Buses; ++i)
{	
	Ybus[i*Num_of_Buses + i].real = Ybus[i*Num_of_Buses + i].real + Gs[i] / BaseMVA;
	Ybus[i*Num_of_Buses + i].imag = Ybus[i*Num_of_Buses + i].imag + Bs[i] / BaseMVA;
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
free(Yf);
free(Yt);
free(CfPrime);
free(CtPrime);
free(Dummy1);
free(Dummy2);


  stop = clock();
 (*CPU_Execution_TIME)[2] = ((double) (stop - start)) / CLOCKS_PER_SEC;

return 1;
}
