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
#include <gsl/gsl_test.h>

#include <gsl/gsl_splinalg.h>

//#include <cblas.h> //Conflict with gsl.
#include <clapack.h>

void call_gsl_gmres(double * A, double * RHS,double * X, const double epsrel, const int compress);
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

//	double work[ldw][restrt + 4];
//	double h[ldh][restrt + 1]; 



/*
int gmres_(n, b, x, restrt, work, ldw, h, ldh, iter, resid, matvec, psolve, 
           info)
   integer *n, *restrt, *ldw, *ldh, *iter, *info;
   doublereal *b, *x, *work, *h, *resid;
   int (*matvec) (), (*psolve) ();
*/

	//gmres_(NBR_POINTS, RHS, X, NBR_POINTS, work, NBR_POINTS, h, ldh, iter, resid, &matvec, &psolve, info);
const double epsrel = 1.0e-14;
const int compress = 0;

call_gsl_gmres(A, RHS,X, epsrel, compress);
//const size_t N, const double epsrel, const int compress


for (int i = 0; i < NBR_POINTS; ++i)
	printf("mu(%d) = %5.16f;\n", i+1,X[i]);
		
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



void call_gsl_gmres(double * A, double * RHS,double * X, const double epsrel, const int compress)
{
	const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
	const size_t max_iter = 10;
  	size_t iter = 0;	
  	gsl_spmatrix *A_gmres = gsl_spmatrix_alloc(NBR_POINTS ,NBR_POINTS); /* triplet format */
 	gsl_spmatrix *B;
  	gsl_vector *b = gsl_vector_alloc(NBR_POINTS);        /* right hand side vector */
  	gsl_vector *u = gsl_vector_calloc(NBR_POINTS);       /* solution vector, u0 = 0 */
  	gsl_splinalg_itersolve *w = gsl_splinalg_itersolve_alloc(T, NBR_POINTS, 0);
  	const char *desc = gsl_splinalg_itersolve_name(w);
  	size_t i;
  	int status;
  	 const double tol = 1.0e-14;


	for (int i = 0; i < NBR_POINTS; ++i)
		for (int j = 0; j < NBR_POINTS; ++j)
			gsl_spmatrix_set(A_gmres, i, j, A[i*NBR_POINTS  + j]);
			
	

    for (int i = 0; i < NBR_POINTS; ++i)
	{
		gsl_vector_set(b, i, RHS[i]);
	}

  	if (compress)
    	B = gsl_spmatrix_compcol(A_gmres);
  	else
    	B = A_gmres;


    /* solve the system */
  do
    {
      status = gsl_splinalg_itersolve_iterate(B, b, tol, u, w);
    }
  while (status == GSL_CONTINUE && ++iter < max_iter);

  gsl_test(status, "%s poisson status s=%d N=%zu", desc, status, NBR_POINTS);

    for (i = 0; i < NBR_POINTS; ++i)
    {
      X[i] = gsl_vector_get(u, i);
    }



  gsl_spmatrix_free(A_gmres);
  gsl_vector_free(b);
  gsl_vector_free(u);
  gsl_splinalg_itersolve_free(w);

  if (compress)
    gsl_spmatrix_free(B);

 // test_poisson(7, 1.0e-1, 1);

 //  test_poisson(543, 1.0e-5, 0);
 //  test_poisson(543, 1.0e-5, 1);


}




