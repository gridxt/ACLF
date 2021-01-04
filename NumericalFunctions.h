#ifndef NUMERICALFUNCTIONS_H_INCLUDED
#define NUMERICALFUNCTIONS_H_INCLUDED

#define _USE_MATH_DEFINES

#include <stdlib.h>
#include <math.h>
#include <omp.h>

// This function converts radian to degree
double rad2deg(double radians);

// This function converts degree to radian
double deg2rad(double degree);

// This function receives two complex numbers in rectangular, multiplies them, returns the real part
double Mult_re(double a, double b, double c, double d);

// This function receives two complex numbers in rectangular, multiplies them, returns the imaginary part
double Mult_img(double a, double b, double c, double d);

// This function gets the real and imaginary part and returns the respective angle in polar form in RADIAN
double rec2pol_ang(double re, double imag);

 // This function gets the real and imaginary part and returns the respective magnitude in polar form
double rec2pol_mag(double re, double imag);

// This function calculates LU decomposition of matrix ALU
 void LU(int MatSize, double *ALU, int Thread_Num);

// This function takes the inverse of matrix ALU and store it in S
void INVS(int MatSize, double *s, double *ALU, int Thread_Num);

#endif // NUMERICALFUNCTIONS_H_INCLUDED
