#include <stdio.h>
#include "Complex.h"


int main(int argc, char const *argv[])
{
	
	double complex z1 = 1.0 + 3.0 * I;
	double complex z2 = 1.0 - 3.0 * I;

	double x = 3.0;

	double y = 1.0;

	printf("Working with complex numbers:\n\v");
	printf("Starting values: z1	= %.2f %+.2fi\t z2 = %.2f %+.2fi\n",creal(z1),cimag(z1),creal(z2),cimag(z2));

	double complex sum =  z1 + z2;
	printf("The sum: z1 + z2 = %.2f %+.2fi\n",creal(sum),cimag(sum));	

	double complex diff = z1 - z2;
	printf("The difference: z1 - z2 = %.2f %+.2fi \n",creal(diff),cimag(diff));

	double complex prod = z1 * z2;
	printf("The product: z1 + z2 = %.2f %+.2fi\n",creal(prod),cimag(prod));

	double complex quot = z1 / z2;
	printf("The quotient: z1 / z2 = %.2f %+.2fi\n",creal(quot),cimag(quot));

	double complex conjugate = conj(z1);
	printf("The conjugate of z1 = %.2f %+.2fi\n",creal(conjugate),cimag(conjugate));

	printf("Squareroot of -1 = %.2fi\n", cimag(csqrt(-1.0)));

	//double complex z3 = Complex(x,y);
	//printf("New complex number z3 = %.2f %+.2fi\n", creal(z3),cimag(z3));

	return 0;
}