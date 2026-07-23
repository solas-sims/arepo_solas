#ifndef STAR_RADIATION_H
#define STAR_RADIATION_H

#include <stdint.h>


#define TAG_RAD 18
#define RAD_TRUNC_FRAC 0.01 
#define MAX_NUM_RAYS 12288 
#define RAY_STACK_SIZE 64

#define NSIDE_MIN 1
#define NSIDE_MAX 32  

/* Dissociation of H2 */
#define SIGMA_DISS 2.47e-18 /* cm^2, dissociation-weighted eff. cross 
                               section (Rollig+07 rate / DB96 flux;
                               Baczynski+15) */

#define F_DISS 0.15 /* dissociation branching per absorption */
#define SIGMA_PUMP (SIGMA_DISS / F_DISS) /* total line absorption   */
 
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

extern double HealpixDirs[MAX_NUM_RAYS][3];
/* 12*NSIDE^2 */
extern int NRays;

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

typedef struct StackEntry
{
  double t_enter;
  double t_exit;
  int node;
} StackEntry;

typedef struct RayPacket
{
  double pos[3];
  double dir[3];
  double t;
  double t_exit;
  double t_maximum;

  /* Bitmask: bit w is SET while band w is still alive */
  /* Cleared when Radiated[w] < RAD_TRUNC_FRAC * Radiated_Init[w] */
  /* When active_bands == 0 the ray is fully absorbed - return immediately */
  uint8_t  active_bands;

  WavebandData Radiated[WAVEBANDS];
  WavebandData Radiated_Init[WAVEBANDS];

  double N_H2; /* accumulated H2 column since source */

  int ray_id;
  int home_task;
  
  /* Pending top-level nodes still to traverse after current domain */
  StackEntry pending[RAY_STACK_SIZE];
  int n_pending;
  int target_node;
  int is_paused;
  
  int nside; /* Current HEALPix nside level */
  int healpix_pixel; /* Pixel index in nested scheme */
} RayPacket;

typedef struct RayWorkStack
{
  long long n; /* Number of rays on this stack */
  long long capacity; /* Allocated capacity */
  RayPacket *rays; /* Ray information */
} RayWorkStack;

typedef struct RayExportBuffer
{
  long long n; /* Number of rays to export */
  long long capacity; /* Allocated capacity */
  int *task; /* Which task to send each ray to */
  RayPacket *rays; /* Ray information */
} RayExportBuffer;

extern struct rad_resultsactiveimported_data
{
  int index; /* Local SphP index on home task */
  WavebandData Radiated[WAVEBANDS];
  double StarMomentumFeed[3];
} *Rad_ResultsActiveImported;

#endif