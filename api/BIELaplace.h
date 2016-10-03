#ifndef BIELAPLACE_H_INCLUDED
#define BIELAPLACE_H_INCLUDED

#define tStart 0.0
#define NBR_R 4 // Resolution in radial direction
#define NBR_T 20 // Resolution in boundary parameterisation
#define NBR_PANELS 5 // Number of GL-panels
#define NBR_POINTS_PER_PANEL 16 // Number of points per GL-panel
#define NBR_POINTS  NBR_POINTS_PER_PANEL *  NBR_PANELS // Useful abbrevation

void init_domain(double complex pz [NBR_R][NBR_T], double complex * ptau, double complex * pzDrops, double complex * pzDropsp, double complex * pzDropspp, double * ppanels,double * ptpar,double * pwDrops);
void solveDensity(double complex * pzDrops, double complex * pzDropsp, double complex * pzDropspp, double * pwDrops,double * RHS,double * mu);
void computeSolution(double * mu, double complex pz [NBR_R][NBR_T], double * pwDrops, double complex * pzDrops, double complex * pzDropsp, double u [NBR_R][NBR_T]);

#endif