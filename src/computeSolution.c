#include <stdio.h>
#include <stdlib.h>
#include "Complex.h"
#include <math.h>
#include "BIELaplace.h"
#include <clapack.h>




void computeSolution(double * pmu, double complex * pz, double * pwDrops, double complex * pzDrops, double complex * pzDropsp, double * pu)
{
	int i, j, k;

	for(i = 0; i < NBR_R; ++i){
		for (int j = 0; j < NBR_T; ++j){
			pu[i * NBR_T +  j] = pu[i * NBR_T +  j] * 0;
			for (int k = 0; k < NBR_PANEL_POINTS; ++k){
				pu[i * NBR_T +  j] =  pu[i * NBR_T +  j] + pmu[k] * pwDrops[k] * cimag(pzDropsp[k] / (pzDrops[k] - pz[i * NBR_T + j]));
			}
			pu[i * NBR_T +  j] = pu[i * NBR_T +  j] * 1.0 / (2.0 * M_PI);
		}
	}
}

