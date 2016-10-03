/*
Calculate density mu from the boundary integral formulation for Laplace
eq. with Dirichlet boundary condition (RHS)
Singular integrals calcualted with limit points

*/


#include <stdio.h>
#include <stdlib.h>
#include "Complex.h"
#include <math.h>
#include "BIELaplace.h"
#include <gsl/gsl_test.h>
#include <gsl/gsl_splinalg.h>
#include <clapack.h>

void call_gsl_gmres(gsl_spmatrix *A_gmres, double * RHS,double * X, const double epsrel, const int compress);


void solveDensity(double complex * pzDrops, double complex * pzDropsp, double complex * pzDropspp, double * pwDrops, double * RHS, double * mu)
{
	int restrt = NBR_POINTS;
	int ldw = NBR_POINTS;
	int ldh = restrt + 1;
	int iter = 1000;
	double resid = 0.00000000001;
	int info;

	gsl_spmatrix *A_gmres = gsl_spmatrix_alloc(NBR_POINTS ,NBR_POINTS); /* triplet format */
	

	for (int i = 0; i < NBR_POINTS; ++i)
			for (int j = 0; j < NBR_POINTS; ++j)
					gsl_spmatrix_set(A_gmres, i, j, (i==j)?(1 + (1/M_PI) * pwDrops[i] * cimag(pzDropspp[i]/(2 * pzDropsp[i]))):((1/M_PI) * pwDrops[j] * cimag(pzDropsp[j]/(pzDrops[j] - pzDrops[i]))));


	const double epsrel = 1.0e-14;
	const int compress = 0;

	call_gsl_gmres(A_gmres, RHS, mu, epsrel, compress);


	gsl_spmatrix_free(A_gmres);		
}






void call_gsl_gmres(gsl_spmatrix *A_gmres, double * RHS,double * mu, const double epsrel, const int compress)
{
	const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
	const size_t max_iter = 10;
  	size_t iter = 0;	
  	
 	gsl_spmatrix *B;
  	gsl_vector *b = gsl_vector_alloc(NBR_POINTS);        /* right hand side vector */
  	gsl_vector *u = gsl_vector_calloc(NBR_POINTS);       /* solution vector, u0 = 0 */
  	gsl_splinalg_itersolve *w = gsl_splinalg_itersolve_alloc(T, NBR_POINTS, 0);
  	const char *desc = gsl_splinalg_itersolve_name(w);
  	size_t i;
  	int status;
  	const double tol = 1.0e-14;

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
      mu[i] = gsl_vector_get(u, i);
    }



  
  gsl_vector_free(b);
  gsl_vector_free(u);
  gsl_splinalg_itersolve_free(w);

  if (compress)
    gsl_spmatrix_free(B);


}




