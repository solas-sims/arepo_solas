#include <math.h>

#include "dust.h"
#include "dust_proto.h"

/*! \brief Grain growth (accretion) rate, Dwek (1998) timescale form.
 *
 *  tau_growth = DUST_GROWTH_TAU_REF_GYR * (n_ref / n_H) * sqrt(T_ref / T)
 *  dM/dt = (M_dust / tau_growth) * (1 - M_dust / M_metal)
 *
 *  The (1 - M_dust/M_metal) factor saturates growth as the dust mass
 *  approaches the cell's total metal budget, so callers do not need a
 *  separate clamp against M_metal for this term alone (dust_cell() still
 *  clamps the combined update against GasMetals as a safety net).
 *
 *  Pure function -- M_dust and M_metal must share a consistent mass unit;
 *  the return value is in that same unit, per Gyr. n_H_cgs/temp_K must be
 *  physical (proper) values in cm^-3 / K.
 *
 *  \param[in] species  dust species index (Phase 1: always DUST_PHASE1_SPECIES)
 *  \param[in] M_dust   current dust mass in the cell
 *  \param[in] M_metal  current total (gas-phase + dust-phase) metal mass in the cell
 *  \param[in] n_H_cgs  hydrogen number density, cm^-3
 *  \param[in] temp_K   gas temperature, K
 *
 *  \return dM/dt due to growth, in [mass unit of M_dust] / Gyr
 */
double dust_growth_rate(int species, double M_dust, double M_metal, double n_H_cgs, double temp_K)
{
  (void)species; /* Phase 1: single species; kept for the Phase 2b API */

  if(M_dust <= 0 || M_metal <= 0 || n_H_cgs <= 0 || temp_K <= 0)
    return 0.0;

  double remaining_fraction = 1.0 - M_dust / M_metal;
  if(remaining_fraction <= 0)
    return 0.0;

  double tau_growth_gyr = DUST_GROWTH_TAU_REF_GYR * (DUST_GROWTH_REF_NH_CGS / n_H_cgs) *
                           sqrt(DUST_GROWTH_REF_TEMP_K / temp_K);

  return (M_dust / tau_growth_gyr) * remaining_fraction;
}
