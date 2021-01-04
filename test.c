#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include "mkl.h"
#include "mkl_types.h"

int main(){

unsigned long int i = 0;
MKL_INT x = 51000;
double *a = NULL;
a = (double *)calloc((size_t)(x*x), sizeof(double));
if (a == NULL) {printf("NULL"); exit(0);}

printf("DBL_MAX     :   %lu\n", (double) DBL_MAX);

for (i = 0; i < x*x; i++)  a[i] = i;

printf("a[%zu] = %f\n", x*x, a[x*x - 1]); 

for (i = 0; i < x*x; i++) {

if (a[i] != i) {
printf("error = %i\n", i);
exit(0);
}


}

return 0;

}