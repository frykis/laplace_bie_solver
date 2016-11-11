#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "BIELaplace.h"
//#include <clapack.h>


/*
Set up computational doimain. To change boundary, introduce new shape functions.
ptau refers to the pointer corresponding to tau. taup and taupp refers to paramterization of boundary with first and second derivative.
*/

void create_grid(double complex * pz){		
	int i,j;
	double t,r;

	for(i = 0; i < NBR_R; i++)
		for (j = 0; j < NBR_T; j++)
		{
			t = 2.0 * M_PI * j /(NBR_T - 1);
			r = i * 0.999 / (NBR_R - 1);
			pz[i * NBR_T + j] = r * (1.0 + 0.3 *  ccos(5.0 * (t + tStart))) * cexp(I * (t + tStart));
		}
}


void tau(double complex * ptau, double * t, int N)
{
	int i;
	for(i = 0; i<N; i++)
		ptau[i] = (1.0 + 0.3 *  ccos(5.0 * (t[i] + tStart))) * cexp(I * (t[i] + tStart));
}



void taup(double complex * ptaup, double * t, int N)
{
	int i;
	for(i = 0; i<N; i++)
		ptaup[i] = (-1.5 * csin(5.0 * (t[i] + tStart)) + I * (1.0 + 0.3 * ccos(5 * (t[i] + tStart)))) * cexp(I * (t[i] + tStart));
}



void taupp(double complex * ptaupp, double * t, int N)
{
	int i;
	for(i = 0; i<NBR_T; i++)
		ptaupp[i] = cexp(I * (t[i] + tStart)) * (-1.0 - 7.8 * ccos(5.0 * (t[i] + tStart)) - (3.0 * I) * csin(5.0 * (t[i] + tStart)));
}


void glwt(double complex * pzDrops, double complex * pzDropsp, double complex * pzDropspp, double * pwDrops, double * ptpar, double first, double last){
	int i,j;
	int N = NBR_POINTS_PER_PANEL - 1;
	int N1 = N + 1;
	int N2 = N + 2;
	double eps = 1e-14;
	double errormax;
	double xu[N1];
	double y[N1];
	double y0[N1];
	double * L = malloc(N1 * N2 * sizeof(double));
	double * Lp = malloc(N1 * sizeof(double));

	for(i = 0; i < N1; i++){
		xu[i] = -1 + 2.0 * i /(N1 - 1);
		Lp[i] = 0;
		y[i] = 0;
		for (j = 0; j < N2; ++j)
		{
			L[i * N2 + j] = 0;
		}
	}


	// Initial guess
	for (i = 0; i < NBR_POINTS_PER_PANEL; ++i)
	{
		y[i] = cos((2* i +1) * M_PI/( 2 * N + 2)) + (0.27 / N1) * sin(M_PI * xu[i] * N / N2);	
		y0[i] = 2;
	}

	
	
	errormax = 2;

	while(errormax > eps)
	{
		for (i = 0; i < N1; ++i)
		{
			L[i*N2] = 1;
			L[i*N2 + 1] = y[i];
			
			for (j = 2; j < N2; ++j)
			{
				L[i*N2 + j ]=( (2*j-1)*y[i]*L[i * N2 + j - 1] - (j-1) * L[i * N2 + j-2] )/j;
			}
			Lp[i] = N2 * (L[i*N2 + N1 - 1] -y[i] * L[i * N2 + N1]) / (1 - pow(y[i],2));
			y0[i] = y[i];
			y[i] = y0[i] - L[i * N2 + N1] / Lp[i];
		}

		



		errormax = 0;
		for (i = 0; i < N1; i++) {
			if (errormax < fabs(y0[i] - y[i]))
				errormax = fabs(y0[i] - y[i]);
		}
	
	} 

	for (i = 0; i < NBR_POINTS_PER_PANEL; ++i){
		ptpar[i] = (first * (1.0 - y[NBR_POINTS_PER_PANEL - 1 - i]) + last * (1.0 + y[NBR_POINTS_PER_PANEL - 1 - i])) * 0.5;
		pwDrops[i] = 1.12890625 * (last - first) / ((1 - pow(y[i],2)) * pow(Lp[i],2));
		tau(&pzDrops[i],&ptpar[i],1);
		taup(&pzDropsp[i],&ptpar[i],1);
		taupp(&pzDropspp[i],&ptpar[i],1);
	} 

	free(L);
	free(Lp);
}

void gl16(double complex * pzDrops, double complex * pzDropsp, double complex * pzDropspp, double * pwDrops, double * ptpar, double complex * ppanels)
{
	int i;
	ppanels[0] = 0;

	for(i = 0; i < NBR_PANELS; ++i)
	{
		ppanels[i + 1] = (i + 1) * 2 * M_PI / NBR_PANELS;		
		//gaussleg requires clapack
		//gaussleg(&pzDrops[i*NBR_POINTS_PER_PANEL],&pzDropsp[i*NBR_POINTS_PER_PANEL],&pzDropspp[i*NBR_POINTS_PER_PANEL],&pwDrops[i*NBR_POINTS_PER_PANEL],&ptpar[i*NBR_POINTS_PER_PANEL],ppanels[i],ppanels[i+1],ptau);
		
		
		glwt(&pzDrops[i*NBR_POINTS_PER_PANEL],&pzDropsp[i*NBR_POINTS_PER_PANEL],&pzDropspp[i*NBR_POINTS_PER_PANEL],&pwDrops[i*NBR_POINTS_PER_PANEL],&ptpar[i*NBR_POINTS_PER_PANEL],ppanels[i],ppanels[i+1]);
		
	}
}




void init_domain(double complex * pz, double complex * pzDrops, double complex * pzDropsp, double complex * pzDropspp, double complex * ppanels,double * ptpar,double * pwDrops)
{

	int i;
	double t[NBR_T];	
	double * pt; 
	for (i = 0; i < NBR_T; i++)
		t[i] = 2.0 * M_PI * i /(NBR_T - 1);

	create_grid(pz);
	

//Obtain GL-16 nodes and weights
	gl16(pzDrops,pzDropsp,pzDropspp,pwDrops,ptpar,ppanels);

	//Transform ppanels to complex points instead of angular paramterisation
	for(i = 0; i<(NBR_PANELS + 1); i++)
		ppanels[i] = (1.0 + 0.3 *  ccos(5.0 * (ppanels[i] + tStart))) * cexp(I * (ppanels[i] + tStart));
	//glwt();
}



