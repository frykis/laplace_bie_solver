#ifndef BIELAPLACE_H_INCLUDED
#define BIELAPLACE_H_INCLUDED

#define S 0.0
#define NBR_R 2 // Resolution in radial direction
#define NBR_T 10 // Resolution in boundary parameterisation
#define NBR_PANELS 2 // Number of GL-panels
#define NBR_POINTS_PER_PANEL 16 // Number of points per GL-panel
#define NBR_POINTS  NBR_POINTS_PER_PANEL *  NBR_PANELS // Useful abbrevation

int init_domain(double complex pz [NBR_R][NBR_T], double complex * ptau, double complex * pzDrops, double complex * pzDropsp, double complex * pzDropspp, double * ppanels,double * ptpar,double * pwDrops);
void solveDensity(double complex * pzDrops, double complex * pzDropsp, double complex * pzDropspp, double * pwDrops,double * RHS);
#endif