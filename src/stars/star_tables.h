#ifndef STAR_TABLES_H
#define STAR_TABLES_H

#include "../main/allvars.h"


extern int Z_COUNT;
extern int M_COUNT;

extern double *Z_VALUES;
extern double *M_VALUES;

extern double *logZ_VALUES;
extern double *logM_VALUES;

extern int **N;

extern double ***Age;
extern double ***Radius;
extern double ***Temperature;

#ifdef WINDS
extern double ***MassLossRate;
#if GRACKLE_CHEMISTRY >= 1
extern double ***HLossRate;
extern double ***HeLossRate;
#endif
#ifdef METALS
extern double ***MetalsLossRate;
#endif
extern double ***WindVelocity;
#endif

#ifdef STAR_RADIATION_ACTIVE
extern WavebandData ***Flux[WAVEBANDS];
#endif

#ifdef SUPERNOVAE
extern double **SN_MassLoss; 
#if GRACKLE_CHEMISTRY >= 1
extern double **SN_HLoss; 
extern double **SN_HeLoss; 
#endif
#ifdef METALS
extern double **SN_MetalsLoss; 
#endif 
#endif

#ifdef AGB 
extern double **AGB_MassLoss; 
#ifdef METALS
extern double **AGB_MetalsLoss; 
#endif 
#endif

#endif