/*
Calculate density mu from the boundary integral formulation for Laplace
eq. with Dirichlet boundary condition (RHS)
Singular integrals calcualted with limit points
The system matrix is symmetric and positive definit.

*/


#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "BIELaplace.h"
//#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_test.h>
//#include <gsl/gsl_splinalg.h>
//#include <clapack.h>

#include <Accelerate/Accelerate.h>



void solveDensity(double complex * pzDrops, double complex * pzDropsp, double complex * pzDropspp, double * pwDrops, double * RHS, double * pmu)
{
  int i,j;
  
  int dim = NBR_PANEL_POINTS;
  int LDA = NBR_PANEL_POINTS;
  int LDB = NBR_PANEL_POINTS;
  int ipiv[NBR_PANEL_POINTS + 1];
  int info;
  int nrhs = 1;

  char trans = 'N';

  double * A = malloc(NBR_PANEL_POINTS * NBR_PANEL_POINTS * sizeof(double));
  
  for (i = 0; i < NBR_PANEL_POINTS; ++i)
    for (j = 0; j < NBR_PANEL_POINTS; ++j)
      A[i * NBR_PANEL_POINTS + j] = (i==j)?(1.0 / 2.0 + (1/(2 * M_PI)) * pwDrops[i] * cimag(pzDropspp[i]/(2 * pzDropsp[i]))):((1 / (2* M_PI)) * pwDrops[i] * cimag(pzDropsp[i]/(pzDrops[i] - pzDrops[j])));

  //Write RHS to pmu. In dgetrs_ is pmu overwritten by the solution.
  for (i = 0; i < NBR_PANEL_POINTS; ++i)
          pmu[i] = RHS[i];

  // Obtain LU-factorization of A. Needed for dgetrs_.
  dgetrf_(&dim, &dim, A, &LDA, ipiv, &info);
  // Solve A*mu = f. The solution mu is in pmu.
  dgetrs_(&trans, &dim, &nrhs, A, &LDA, ipiv, pmu, &LDB, &info);

  


free(A);
}





