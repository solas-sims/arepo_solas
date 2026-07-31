#include <math.h>

#include "dust.h"
#include "dust_proto.h"

/*! \brief Thermal sputtering rate, Tsai & Mathews (1995) timescale form.
 *
 *  tau_sp = DUST_SPUTTER_TAU_REF_GYR * (n_ref / n_H) * [(T_ref/T)^omega + 1]
 *  dM/dt = -M_dust / tau_sp
 *
 *  Pure function -- M_dust and the return value share a consistent mass
 *  unit, per Gyr. n_H_cgs/temp_K must be physical (proper) values in
 *  cm^-3 / K.
 *
 *  \param[in] species  dust species index (Phase 1: always DUST_PHASE1_SPECIES)
 *  \param[in] M_dust   current dust mass in the cell
 *  \param[in] n_H_cgs  hydrogen number density, cm^-3
 *  \param[in] temp_K   gas temperature, K
 *
 *  \return dM/dt due to sputtering (<= 0), in [mass unit of M_dust] / Gyr
 */
double dust_sputtering_rate(int species, double M_dust, double n_H_cgs, double temp_K)
{
  (void)species; /* Phase 1: single species; kept for the Phase 2b API */

  if(M_dust <= 0 || n_H_cgs <= 0 || temp_K <= 0)
    return 0.0;

  double tau_sp_gyr = DUST_SPUTTER_TAU_REF_GYR * (DUST_SPUTTER_REF_NH_CGS / n_H_cgs) *
                       (pow(DUST_SPUTTER_REF_TEMP_K / temp_K, DUST_SPUTTER_OMEGA) + 1.0);

  return -M_dust / tau_sp_gyr;
}
