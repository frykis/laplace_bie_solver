#include <stdio.h>
#include <stdlib.h>
#include "Complex.h"
#include <math.h>
#include "BIELaplace.h"
#include <clapack.h>




void computeSolution(double * mu, double complex pz [NBR_R][NBR_T], double * pwDrops, double complex * pzDrops, double complex * pzDropsp, double u [NBR_R][NBR_T])
{
	//Not working correctly. Getting wrong pz compared to Matlab.
	int i,j,k;
	double t,r;
	for(i = 0; i < NBR_R; ++i)
		for (int j = 0; j < NBR_T; ++j)
			u[i][j] = u[i][j] * 0;


	for(i = 0; i < NBR_R; ++i)
		for (int j = 0; j < NBR_T; ++j)
			for (int k = 0; k < NBR_POINTS; ++k)
				u[i][j] =  u[i][j] + mu[k] * pwDrops[k] * cimag(pzDropsp[k] / (pzDrops[k] - pz[i][j]));


	for(i = 0; i < NBR_R; ++i)
		for (int j = 0; j < NBR_T; ++j)
			u[i][j] = u[i][j] * 1.0 / (2.0 * M_PI);

	
}

