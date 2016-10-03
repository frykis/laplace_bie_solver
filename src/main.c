#include <stdio.h>
#include <stdlib.h>
#include "Complex.h"
#include <math.h>
#include "BIELaplace.h"
#include <clapack.h>



int main(int argc, char const *argv[])
{
	int i;
	int N = NBR_POINTS_PER_PANEL * NBR_PANELS;
	
   	double complex pz [NBR_R][NBR_T];
   	double u [NBR_R][NBR_T];

	double complex * ptau = malloc(NBR_T * sizeof(complex double));
	double complex * pzDrops = malloc(NBR_PANELS * NBR_POINTS_PER_PANEL * sizeof(complex double));
	double complex * pzDropsp = malloc(NBR_PANELS * NBR_POINTS_PER_PANEL * sizeof(complex double));
	double complex * pzDropspp = malloc(NBR_PANELS * NBR_POINTS_PER_PANEL * sizeof(complex double));
	
	double * RHS = malloc(NBR_POINTS * sizeof(double));
	double * ppanels = malloc((NBR_PANELS + 1) * sizeof(double));
	double * ptpar = malloc(NBR_POINTS_PER_PANEL * NBR_PANELS * sizeof(double));
	double * pwDrops = malloc(NBR_POINTS_PER_PANEL * NBR_PANELS * sizeof(double));
	double * mu = malloc(NBR_POINTS * sizeof(double));

	double complex zsrc1 = 1.5 + 1.5 * I;
	double complex zsrc2 = -0.25 + 1.5 * I;
	double complex zsrc3 = -0.5 - 1.5 * I;
	

	
	init_domain(pz,ptau,pzDrops,pzDropsp,pzDropspp,ppanels,ptpar,pwDrops);


	for (int i = 0; i < NBR_POINTS; ++i)
	{
		RHS[i] = creal( 1 / (pzDrops[i]-zsrc1) + 1 / (pzDrops[i]-zsrc2) + 1 / (pzDrops[i] - zsrc3));
	}


	solveDensity(pzDrops, pzDropsp, pzDropspp, pwDrops, RHS, mu);
					
	computeSolution(mu, pz, pwDrops, pzDrops, pzDropsp, u);

	for(i = 0; i < NBR_R; ++i)
		for (int j = 0; j < NBR_T; ++j)
			printf("u(%d,%d) = %5.16f\n", i + 1,j + 1,u[i][j]);

//Free allocated space.
	free(ptau);
	free(pzDrops);
	free(pzDropsp);
	free(pzDropspp);
	free(ppanels);
	free(ptpar);
	free(pwDrops);
	free(RHS);
	free(mu);



	return 0;
}

