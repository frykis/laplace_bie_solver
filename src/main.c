#include <stdio.h>
#include <stdlib.h>
#include "Complex.h"
#include <math.h>
#include "BIELaplace.h"
//#include <gsl/gsl_splinalg.h>

//#include <Accelerate/Accelerate.h> // OSX-specific for LAPACK and BLAS.
//#include <cblas.h> //Conflict with gsl.
#include <clapack.h>



int main(int argc, char const *argv[])
{
int i;
	int N = NBR_POINTS_PER_PANEL * NBR_PANELS;
	
   	double complex pz [NBR_R][NBR_T];

	double complex * ptau = malloc(NBR_T * sizeof(complex double));
	double complex * pzDrops = malloc(NBR_PANELS * NBR_POINTS_PER_PANEL * sizeof(complex double));
	double complex * pzDropsp = malloc(NBR_PANELS * NBR_POINTS_PER_PANEL * sizeof(complex double));
	double complex * pzDropspp = malloc(NBR_PANELS * NBR_POINTS_PER_PANEL * sizeof(complex double));
	
	double * RHS = malloc(NBR_POINTS * sizeof(double));
	double * ppanels = malloc((NBR_PANELS + 1) * sizeof(double));
	double * ptpar = malloc(NBR_POINTS_PER_PANEL * NBR_PANELS * sizeof(double));
	double * pwDrops = malloc(NBR_POINTS_PER_PANEL * NBR_PANELS * sizeof(double));
	

	i = init_domain(pz,ptau,pzDrops,pzDropsp,pzDropspp,ppanels,ptpar,pwDrops);

	double complex zsrc1 = 1.5 + 1.5 * I;
	double complex zsrc2 = -0.25 + 1.5 * I;
	double complex zsrc3 = -0.5 - 1.5 * I;
	

	for (int i = 0; i < NBR_POINTS; ++i)
	{
		RHS[i] = creal( 1 / (pzDrops[i]-zsrc1) + 1 / (pzDrops[i]-zsrc2) + 1 / (pzDrops[i] - zsrc3));
	}

	//pwDrops OK!
	//ptpar OK!
	//pzDrops OK!
	//A OK!


	solveDensity(pzDrops, pzDropsp, pzDropspp, pwDrops, RHS);
	
	

//Free allocated space.
	free(ptau);
	free(pzDrops);
	free(pzDropsp);
	free(pzDropspp);
	free(ppanels);
	free(ptpar);
	free(pwDrops);
	free(RHS);




	return 0;
}

/*function mu_lapl = mubie_lapl(N,zDrops,zpDrops,zppDrops,wDrops,RHS)
% Calculate density mu from the boundary integral formulation for Laplace
% eq. with Dirichlet boundary condition (RHS)
% Singular integrals calcualted with limit points
A = eye(N,N);
for i=1:N
    A(i,:) = A(i,:) + (1/pi*wDrops.*imag(zpDrops./(zDrops-zDrops(i))))';
end
d = 1 + 1/pi*wDrops.*imag(zppDrops./(2*zpDrops));
b = 2*RHS(zDrops);%imag(tauk.^2./(tauk-zp));
A(logical(eye(size(A)))) = 0;
A = A + diag(d,0);
mu_lapl = gmres(A,b,[],1e-13,N);
% mu_lapl = A\b;

end*/