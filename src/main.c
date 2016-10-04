#include <stdio.h>
#include <stdlib.h>
#include "Complex.h"
#include <math.h>
#include "BIELaplace.h"
#include <clapack.h>



int main(int argc, char const *argv[])
{
	int i;

	double umax, errormax;

	double * pumax;	

	umax = 0; // Used to obtain relative error.
	pumax = &umax;

	double complex zsrc1 = 1.5 + 1.5 * I;
	double complex zsrc2 = -0.25 + 1.5 * I;
	double complex zsrc3 = -0.5 - 1.5 * I;
   	
   	double complex * pz = malloc(NBR_DOMAIN_POINTS * sizeof(complex double));
	double complex * ptau = malloc(NBR_T * sizeof(complex double));
	double complex * pzDrops = malloc(NBR_PANEL_POINTS * sizeof(complex double));
	double complex * pzDropsp = malloc(NBR_PANEL_POINTS * sizeof(complex double));
	double complex * pzDropspp = malloc(NBR_PANEL_POINTS * sizeof(complex double));
	double complex * ppanels = malloc((NBR_PANELS + 1) * sizeof(complex double));
	
	double * RHS = malloc(NBR_PANEL_POINTS * sizeof(double));
	double * ptpar = malloc(NBR_PANEL_POINTS * sizeof(double));
	double * pwDrops = malloc(NBR_PANEL_POINTS * sizeof(double));
	double * pmu = malloc(NBR_PANEL_POINTS * sizeof(double));
	double * pu = malloc(NBR_DOMAIN_POINTS * sizeof(double));
	double * pu_spec = malloc(NBR_DOMAIN_POINTS * sizeof(double));
	double * pu_ana = malloc(NBR_DOMAIN_POINTS * sizeof(double));
	double * perrorvec = malloc(NBR_DOMAIN_POINTS * sizeof(double));
	
	init_domain(pz,ptau,pzDrops,pzDropsp,pzDropspp,ppanels,ptpar,pwDrops);
	init_function(RHS, pu_ana, pzDrops, pz, pumax);
	solveDensity(pzDrops, pzDropsp, pzDropspp, pwDrops, RHS, pmu);
	computeSolution(pmu, pz, pwDrops, pzDrops, pzDropsp, pu);
	specialquadlapl(pu_spec, pu, pmu, pz, pzDrops, pzDropsp, pwDrops, ppanels);	
	computeError(perrorvec, pu, pu_spec, pu_ana, pumax);

//Free allocated space.
	free(ptau);
	free(pzDrops);
	free(pzDropsp);
	free(pzDropspp);
	free(ppanels);
	free(ptpar);
	free(pwDrops);
	free(RHS);
	free(pmu);
	free(pu);
	free(pu_spec);
	free(pu_ana);
	free(perrorvec);

	return 0;
}

