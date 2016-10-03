#include <stdio.h>
#include <stdlib.h>
#include "Complex.h"
#include <math.h>
#include "BIELaplace.h"
#include <clapack.h>


/*
Set up computational doimain. To change boundary, introduce new shape functions.
ptau refers to the pointer corresponding to tau. taup and taupp refers to paramterization of boundary with first and second derivatine.
*/
void tau(double complex *, double * t,int N); 
void taup(double complex *, double * t,int N);
void taupp(double complex *, double * t,int N);
void create_grid(double complex[NBR_R][NBR_T]);
void gl16(double complex * pzDrops, double complex * pzDropsp, double complex * pzDropspp, double * pwDrops, double * ptpar, double * ppanels, double complex * ptau);
void gaussleg(double complex * pzDrops, double complex * pzDropsp, double complex * pzDropspp, double * pwDrops, double * ptpar, double,double,double complex * ptau);
void transposeMatrix(double *, double *,int,int);
void printArrayD(double * A, int N, int M);
void printArrayC(double complex * A, int N, int M);
//double tStart =0.0;


void init_domain(double complex pz [NBR_R][NBR_T], double complex * ptau, double complex * pzDrops, double complex * pzDropsp, double complex * pzDropspp, double * ppanels,double * ptpar,double * pwDrops)
{

	int i;
	double t[NBR_T];	
	for (int i = 0; i < NBR_T; i++)
			t[i] = 2.0 * M_PI * i /(NBR_T - 1);

	create_grid(pz);
	tau(ptau,t,NBR_T);

//Obtain GL-16 nodes and weights


	
	gl16(pzDrops,pzDropsp,pzDropspp,pwDrops,ptpar,ppanels,ptau);

}


void gl16(double complex * pzDrops, double complex * pzDropsp, double complex * pzDropspp, double * pwDrops, double * ptpar, double * ppanels, double complex * ptau)
{
	int i;

	for(i = 0; i< NBR_PANELS + 1; ++i){
		ppanels[i] = i * 2 * M_PI / NBR_PANELS;
	}

	for(i = 0; i < NBR_PANELS; ++i){
		gaussleg(&pzDrops[i*NBR_POINTS_PER_PANEL],&pzDropsp[i*NBR_POINTS_PER_PANEL],&pzDropspp[i*NBR_POINTS_PER_PANEL],&pwDrops[i*NBR_POINTS_PER_PANEL],&ptpar[i*NBR_POINTS_PER_PANEL],ppanels[i],ppanels[i+1],ptau);
	}

}


void gaussleg(double complex * pzDrops, double complex * pzDropsp, double complex * pzDropspp, double * pwDrops, double * ptpar, double first, double last,double complex * ptau){
	
	int i,j,INFO;
	double beta[NBR_POINTS_PER_PANEL-1];
	double diag[NBR_POINTS_PER_PANEL];
	double Z[NBR_POINTS_PER_PANEL*NBR_POINTS_PER_PANEL];
	double WORK[2*NBR_POINTS_PER_PANEL-2];
	// Parameters for dstev_. c4 gives that we obtain eigenvalues and eigenvectors. N is needed as the value must be passed by reference.
	char dstevParameter='V';
	int N = NBR_POINTS_PER_PANEL;
	double temp;


	for(i = 1; i < NBR_POINTS_PER_PANEL;++i)
		beta[i-1] = 0.5 * pow(1.0 - pow(2.0 * i,-2.0),-0.5);


	for (int i = 0; i < NBR_POINTS_PER_PANEL; ++i)
	{
		diag[i] = 0.0;
		for (int j = 0; j < NBR_POINTS_PER_PANEL; ++j)
			Z[i*NBR_POINTS_PER_PANEL + j] = (i == j)?1.0:0.0;
	}



//Obtain eigenvalues and eigenvectors for symmetric tridiagonal real matrix.
dstev_(&dstevParameter,&N,diag,beta,Z,&N,WORK,&INFO);



for (int i = 0; i < NBR_POINTS_PER_PANEL; ++i)
{
	ptpar[i] = (first * (1 - diag[i]) + last * (1 + diag[i])) * 0.5;
	pwDrops[i] = (last - first) * 0.5 * 2 * pow( Z[i*NBR_POINTS_PER_PANEL],2);
	tau(&pzDrops[i],&ptpar[i],1);
	taup(&pzDropsp[i],&ptpar[i],1);
	taupp(&pzDropspp[i],&ptpar[i],1);

}

}


void create_grid(double complex pz[NBR_R][NBR_T]){		
	int i,j;
	double t,r;

	for(i = 0; i < NBR_R; i++){
		for (int j = 0; j < NBR_T; j++)
		{
			t = 2.0 * M_PI * j /(NBR_T - 1);
			r = i * 0.999 / (NBR_R - 1);
			pz[i][j] = r * (1.0 + 0.3 *  ccos(5.0 * (t + tStart))) * cexp(I * (t + tStart));
		}
	}
}


void tau(double complex * ptau, double * t, int N)
{
	int i;
	for(int i = 0; i<N; i++){
		ptau[i] = (1.0 + 0.3 *  ccos(5.0 * (t[i] + tStart))) * cexp(I * (t[i] + tStart));

}
}



void taup(double complex * ptaup, double * t, int N)
{
	int i;
	for(int i = 0; i<N; i++){
		ptaup[i] = (-1.5 * csin(5.0 * (t[i] + tStart)) + I * (1.0 + 0.3 * ccos(5 * (t[i] + tStart)))) * cexp(I * (t[i] + tStart));
	}	
}



void taupp(double complex * ptaupp, double * t, int N)
{
	int i;
	for(int i = 0; i<NBR_T; i++){
		ptaupp[i] = cexp(I * (t[i] + tStart)) * (-1.0 - 7.8 * ccos(5.0 * (t[i] + tStart)) - (3.0 * I) * csin(5.0 * (t[i] + tStart)));
	}
}


void transposeMatrix(double * A, double * AT, int N, int M)
{
	int i,j;
	for(i=0; i<N; ++i)
    	for(j=0; j<M; ++j)
    		AT[j*N + i]=A[i * N + j];
}

void printArrayD(double * A, int N, int M)
{
	int i,j;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			printf("A[%d][%d] = %5.4f\t",i,j,A[i*M + j]);
		}
		printf("\n");
	}

}

void printArrayC(double complex * A, int N, int M)
{
	int i,j;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			printf("A[%d][%d] = %5.4f %+5.4fi \t",i,j,creal(A[i*M + j]),cimag(A[i*M + j]));
		}
		printf("\n");
	}

}
