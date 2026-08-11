#ifndef STAR_PARTICLE_H
#define STAR_PARTICLE_H

#include "../main/allvars.h"

#define NBINS 114
#define MMIN 0.10
#define MMAX 120.0

#define N_CDF_BINS 10000

enum IMFType
{
    POPII = 0,
#ifdef POPIII_SF
    POPIII,
#endif
    N_IMF_TYPES  /* 2 #ifdef POPIII_SF, else 1 */
};

extern double cdf_masses[N_IMF_TYPES][N_CDF_BINS + 1];   
extern double cdf_values[N_IMF_TYPES][N_CDF_BINS + 1];   

#if defined(STAR_PARTICLES) && STAR_PARTICLES < 2
extern double StarMassBins[NBINS + 1];
extern double StarMeanMassInBins[N_IMF_TYPES][NBINS];
#endif 

#if STAR_PARTICLES == 0

#include <gsl/gsl_rng.h>

extern gsl_rng *rng;

extern double bin_imf[N_IMF_TYPES][NBINS];
#endif

#endif