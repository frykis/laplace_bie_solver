/*
Calculate density mu from the boundary integral formulation for Laplace
eq. with Dirichlet boundary condition (RHS)
Singular integrals calcualted with limit points

*/


#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include "Complex.h"
#include <math.h>
#include "BIELaplace.h"

#include <cblas.h>
#include <clapack.h>


void matvec(double alpha, double * X, double beta, double * Y);
void psolve(double * X, double * Y);
double * A;



void solveDensity(double complex * pzDrops, double complex * pzDropsp, double complex * pzDropspp, double * pwDrops, double * RHS)
{
	int restrt = NBR_POINTS;
	int ldw = NBR_POINTS;
	int ldh = restrt + 1;
	int iter = 1000;
	double resid = 0.00000000001;
	int info;


	A = malloc(NBR_POINTS * NBR_POINTS * sizeof(double));
	for (int i = 0; i < NBR_POINTS; ++i)
		for (int j = 0; j < NBR_POINTS; ++j)
			A[i*NBR_POINTS  + j] = (i == j)?1.0:0.0;

	for (int i = 0; i < NBR_POINTS; ++i)
			for (int j = 0; j < NBR_POINTS; ++j)
				A[i*NBR_POINTS + j] = A[i*NBR_POINTS + j] + (1/M_PI) * pwDrops[j] * cimag(pzDropsp[j]/(pzDrops[j] - pzDrops[i]));


	for (int i = 0; i < NBR_POINTS; ++i)
	{
		for (int j = 0; j < NBR_POINTS; ++j)
		{
			if (i == j)
			{
				A[i*NBR_POINTS  + j] = 0.0;
				A[i*NBR_POINTS  + j] = 1 + (1/M_PI) * pwDrops[i] * cimag(pzDropspp[i]/(2 * pzDropsp[i]));
			}
		}
	}

	//A OK!



	// Initial guess, will contain solution after GMRES.
	double * X = malloc(NBR_POINTS * sizeof(double));
	for (int i = 0; i < NBR_POINTS; ++i)
			X[i] = 1.0;

	double work[ldw][restrt + 4];
	double h[ldh][restrt + 1];

/*
int gmres_(n, b, x, restrt, work, ldw, h, ldh, iter, resid, matvec, psolve, 
           info)
   integer *n, *restrt, *ldw, *ldh, *iter, *info;
   doublereal *b, *x, *work, *h, *resid;
   int (*matvec) (), (*psolve) ();
*/

	gmres_(NBR_POINTS, RHS, X, NBR_POINTS, work, NBR_POINTS, h, ldh, iter, resid, &matvec, &psolve, info);

for (int i = 0; i < NBR_POINTS; ++i)
	printf("mu(%d) = %5.16f\n", i,X[i]);
		
	free(A);
	free(X);
}
/*
 Calculate density mu from the boundary integral formulation for Laplace
 eq. with Dirichlet boundary condition (RHS)
 Singular integrals calcualted with limit points
*/



void matvec(double alpha, double * X, double beta, double * Y){
//Matvec:y := alpha*A*x + beta*y,CALL MATVEC( ALPHA, X, BETA, Y )
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NBR_POINTS, 1, NBR_POINTS, alpha, A, NBR_POINTS, X, 1, beta, Y, 1);

}


/*
Dummy function. Replace M with preconditioner if necessary.
preconditioner solve routine for the linear system

               M*x = b,

          where x and b are vectors, and M a matrix. Vector b must
          remain unchanged.
          The solution is over-written on vector x.
*/
void psolve(double * X, double * Y);








