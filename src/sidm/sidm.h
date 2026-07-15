#ifndef SIDM_H
#define SIDM_H

#include "../main/allvars.h"

#ifdef SIDM

/*! \brief Per-DM-particle SIDM data, in a dedicated side array rather
 *  than on the shared particle_data struct (P[]).
 *
 *  Mirrors Star_Particle_Data / SP[] (src/stars/star.h) exactly in
 *  mechanism: PIndex is a back-reference into P[], kept in sync with
 *  P[i].SIDMID (the forward reference) by the same swap-compaction
 *  logic domain_exchange.c already uses for STARS/BLACKHOLES -- see
 *  the Type==1 branch there.
 *
 *  This existed as direct P[] fields in earlier versions of this
 *  module; moved to a side array because DM is typically a large
 *  fraction of total particles in a mixed baryon+DM run, so a
 *  dedicated array avoids every gas/star/BH particle carrying fields
 *  it never uses. In a pure N-body (100% DM) run this costs slightly
 *  more than the old direct-P[]-fields approach (the PIndex
 *  back-reference is pure overhead when there's no non-DM population
 *  to exclude) -- a deliberate, documented trade-off, not an oversight.
 */
typedef struct DM_Particle_Data
{
  MyIDType PIndex; /*!< back-reference into P[] -- see PIndex naming note:
                    * chosen over the star/BH convention's "PID" name
                    * specifically because "PID" reads as "particle ID"
                    * (a persistent identifier) when it's actually just
                    * an array index, which is confusing to a reader
                    * meeting this pattern for the first time. Same
                    * mechanism as Star_Particle_Data.PID, different name. */

  MyFloat SidmDensity;     /*!< local DM density estimate, elastic v1 */
  MyFloat SidmHsml;        /*!< adaptive smoothing length / search radius */
  MyFloat SidmVelDisp;     /*!< local 1D DM velocity dispersion, mean-subtracted */
  int     SidmNumNgb;      /*!< neighbours found within SidmHsml, for Hsml convergence */
  integertime SidmLastScatterTime; /*!< diagnostic / rate-limiting */
  int     SidmScatterFlag;         /*!< set by Monte Carlo step, consumed by kick routine */
  int     SidmScatterCount;        /*!< cumulative count of scatters this particle has undergone --
                                     * unlike SidmScatterFlag (reset-able, single-step), this only
                                     * ever increments, so differencing it between two snapshots
                                     * gives a genuine empirical scattering rate for validation
                                     * against the analytic Gamma = rho*(sigma/m)*v_rel prediction. */
} DM_Particle_Data;

extern DM_Particle_Data *DMSP;
extern int               NumDM;

/*! DMPS(i): given a P[]-index i, get its DM_Particle_Data (forward,
 *  mirrors SPP(i) = SP[P[i].SID]).
 *  PDMS(i): given a DMSP[]-index i, get its particle_data (backward,
 *  mirrors PPS(i) = P[SP[i].PID]). */
#define DMPS(i) DMSP[P[i].SIDMID]
#define PDMS(i) P[DMSP[i].PIndex]

void sidm_density(void);
void sidm_scatter(int timebin);

/*! Safety factor bounding per-step scattering probability (rate*dt) for
 *  the SIDM timestep criterion in get_timestep_gravity() (timestep.c).
 *  Hardcoded for v1 rather than a param.txt entry -- promote to a
 *  runtime parameter later if tuning turns out to be needed. */
#define SIDM_TIMESTEP_SAFETY_FACTOR 0.1

#endif /* #ifdef SIDM */

#endif /* #ifndef SIDM_H */
