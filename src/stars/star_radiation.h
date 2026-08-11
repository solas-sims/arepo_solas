#ifndef STAR_RADIATION_H
#define STAR_RADIATION_H

#include <stdint.h>


#define TAG_RAD 18

#define RAD_TRUNC_FRAC 0.01

/* Safety cap on the number of Voronoi cells a single ray may cross in one call to raytrace_voronoi() /* 
/* Only trips on degenerate geometry */
#define RAY_MAX_CELL_STEPS (1 << 22)

/* Relative tolerance (in units of the host cell radius) used to break ties
   between coincident bisectors at cell edges/vertices */
#define RAY_TOL 1.0e-10

/* Healpix: Multiples of 2! */
#define NSIDE_MIN 1
#define NSIDE_MAX 128

/* Starting number of rays */
#define NRays (12 * NSIDE_MIN * NSIDE_MIN)

/* Dissociation of H2 */
#define SIGMA_DISS 2.47e-18 /* cm^2, dissociation-weighted eff. cross
                               section (DB96, Baczynski+15) */

#define F_DISS 0.15 /* dissociation branching per absorption */
#define SIGMA_PUMP (SIGMA_DISS / F_DISS) /* total line absorption */

#define H2_SHIELD_B5 3.0 /* Doppler b / (km/s); fixed */

/* Shielding log table parameters */
#define H2TAB_N 1024
#define H2TAB_LOGNMIN 11.0 /* log10 N_H2 [cm^-2]: f_sh = 1 below */
#define H2TAB_LOGNMAX 24.0 /* f_sh negligible above */

/* Change at init_rays */
#define ALL_BANDS_ACTIVE ((uint8_t)((1u << WAVEBANDS) - 1u))
#define NO_IR_ACTIVE ((uint8_t)(ALL_BANDS_ACTIVE & ~(1u << INFRARED)))
#define NO_IONIZING_ACTIVE ((uint8_t)(ALL_BANDS_ACTIVE & ~(1u << IONIZING_HI) & ~(1u << IONIZING_HeI) & ~(1u << IONIZING_HeII)))
#define ONLY_IONIZING_ACTIVE ((uint8_t)(ALL_BANDS_ACTIVE & ((1u << IONIZING_HI) | (1u << IONIZING_HeI) | (1u << IONIZING_HeII))))

/*
 * INFRARED: inf A - 12398.4 A (0 eV - 1 eV)
 *
 * OPTICAL: 12398.4 A - 2066.4 A (1 eV - 6 eV)
 *
 * ULTRAVIOLET: 2066.4 A - 1107.0 A (6 eV - 11.2 eV)
 *
 * LYMAN_WERNER: 1107.0 A - 911.6 A (11.2 eV - 13.6 eV)
 *
 * IONIZING_HI: 911.6 A - 504.0 A (13.6 eV - 24.6 eV)
 *
 * IONIZING_HeI: 504.0 A - 227.9 A (24.6 eV - 54.4 eV)
 *
 * IONIZING_HeII: 227.9 A - 0 A (54.4 eV - inf eV)
 */
typedef enum
{ INFRARED = 0,
  OPTICAL,
  ULTRAVIOLET,
  LYMAN_WERNER,
  IONIZING_HI,
  IONIZING_HeI,
  IONIZING_HeII,
  WAVEBANDS
} Waveband;

#if WAVEBANDS > 8
#error "Active_bands is uint8_t but WAVEBANDS > 8 - use uint16_t instead"
#endif

typedef struct WavebandData
{
  double Photons;
  double Energy;
} WavebandData;

extern double Kappa_E[WAVEBANDS];
extern double Kappa_N[WAVEBANDS];

extern double ReradiatedFraction[WAVEBANDS];

/*
 * RayPacket
 * ---------
 * A ray is fully described by the cell it currently occupies plus its
 * entry point into that cell. Because the Voronoi cell is exactly the
 * intersection of the half spaces
 *
 *     (x - m_ij) . d_ij <= 0 ,   d_ij = s_j - s_i ,  m_ij = (s_i + s_j)/2
 *
 * over the face-defining neighbours listed in DC(i), the exit point and the
 * next cell follow from a single min-reduction over that list.
 */
typedef struct RayPacket
{
  MyIDType star_id;
    
  /* Cell currently occupied: local SphP index on the owning task */
  int cell;

  double pos[3]; /* Position relative to the generator of `cell` */
  double dir[3]; /* Unit propagation direction (global frame) */

  double t; /* Path length accumulated since the star */
  double t_maximum; /* Hard stop */

  int nside; /* Current HEALPix nside level */
  int healpix_pixel; /* Current HEALPix pixel (NESTED) */
  unsigned long long rotation_seed; /* HEALPix basis rotation seed (per star) */
  uint8_t needs_locate; /* Locate flag for split rays */

  /* Bitmask: bit w is SET while band w is still alive */
  /* Cleared when Radiated[w] < RAD_TRUNC_FRAC * Radiated_Init[w] */
  /* When active_bands == 0 the ray is fully absorbed - return immediately */
  uint8_t active_bands;

  /* Ray energy and photons */
  WavebandData Radiated[WAVEBANDS];
  WavebandData Radiated_Init[WAVEBANDS];

  /* Accumulated H2 column since source */
  double N_H2;

  /* Ray bookkeeping */
  int ray_id;
  int home_task;
} RayPacket;

typedef struct RayWorkStack
{
  long long n; /* Number of rays on this stack */
  long long capacity; /* Allocated capacity */
  RayPacket *rays;
} RayWorkStack;

typedef struct RayExportBuffer
{
  long long n; /* Number of rays to export */
  long long capacity; /* Allocated capacity */
  int *task; /* Where to send each ray */
  RayPacket *rays;
} RayExportBuffer;

#endif