#include "NumericalFunctions.h"

// This function converts radian to degree
double rad2deg(double radians) {
    return radians * (180.0 / M_PI);
}

// This function converts degree to radian
double deg2rad(double degree) {
    return degree * (M_PI / 180.0);
}

// This function receives two complex numbers in rectangular, multiplies them, returns the real part
double Mult_re(double a, double b, double c, double d) {
    return ( (a*c) - (b*d) );
}

// This function receives two complex numbers in rectangular, multiplies them, returns the imaginary part
double Mult_img(double a, double b, double c, double d) {
    return ( (a*d) + (b*c) );
}

// This function gets the real and imaginary part and returns the respective angle in polar form in RADIAN
double rec2pol_ang(double re, double imag) {
    if(re<0 && imag>=0) return atan(imag/re) + acos(-1);
    if(re<0 && imag<0) return atan(imag/re) - acos(-1);
    if(re==0 && imag>0) return acos(-1)/2;
    if(re==0 && imag<0) return -1*acos(-1)/2;
    if(re>0) return atan(imag/re);
    if(re==0 && imag==0) return 0;
 return -1; // to avoid compiler warning
 }

 // This function gets the real and imaginary part and returns the respective magnitude in polar form
double rec2pol_mag(double re, double imag){
 return sqrt((re*re)+(imag*imag)); // sqrt(re^2 + imag^2)
 }

// This function calculates LU decomposition of matrix ALU
 void LU(int MatSize, double *ALU, int Thread_Num)
{
	int i,j,k,n;
	int Chunk_Size;
	double x;
	n = MatSize - 1;
	
	
	Chunk_Size = round(n-1 / Thread_Num);
	//#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(Thread_Num) private(k, i, j) shared(Chunk_Size, n, x, ALU, MatSize)
    for(k=0;k<=n-1;k++) {		
		for(j=k+1;j<=n;j++) {
        x=ALU[j*MatSize+k]/ALU[k*MatSize + k];
            for(i=k;i<=n;i++) {ALU[j*MatSize + i]=ALU[j*MatSize + i]-x*ALU[k*MatSize + i];}
        ALU[j*MatSize + k]=x;		
	  }    
	}
			
}

// This function takes the inverse of matrix ALU and store it in S
void INVS(int MatSize, double *s, double *ALU, int Thread_Num)
{
int i,j,m,n, k, l;
int xx;
int Chunk_Size;
n = MatSize-1;
double x;
double *d = (double *)calloc(MatSize, sizeof(double));
double *y = (double *)calloc(MatSize, sizeof(double));

Chunk_Size = round(n / Thread_Num);
//#pragma omp parallel for default(none) schedule(guided, Chunk_Size) num_threads(Thread_Num) private(m, i, j) shared(Chunk_Size, x, y, n, d, s, ALU, MatSize)
for (m = 0; m <= n; m++) {
	d[m] = 1.0;
	if (m != 0) d[m - 1] = 0;

	for (i = 0; i <= n; i++) {
		x = 0;
		for (j = 0; j <= i - 1; j++) {
			x = x + ALU[i*MatSize + j] * y[j];
		}
		y[i] = (d[i] - x);
	}

	for (i = n; i >= 0; i--) {
		x = 0.0;
		for (j = i + 1; j <= n; j++) x = x + ALU[i*MatSize + j] * s[j*MatSize + m];
		s[i*MatSize + m] = (y[i] - x) / ALU[i*MatSize + i];
	}

}

free(d);
free(y);

}
