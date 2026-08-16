#include <math.h>
#include <mpi.h>
#include <stdlib.h>

#include "../main/allvars.h"
#include "../main/proto.h"
#include "dust.h"
#include "dust_proto.h"

#ifdef DUST

/*! \brief Physical hydrogen number density and gas temperature for a cell.
 *
 *  Mirrors the unit handling in DoCooling()/cool_cell() (src/cooling/cooling.c)
 *  so the dust rate laws see the same physical density/temperature cooling
 *  itself uses.
 */
static void dust_get_nH_and_temp(int i, double *n_H_cgs, double *temp_K)
{
  double dens_code = SphP[i].Density * All.cf_a3inv;
  double dens_cgs   = dens_code * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;

  *n_H_cgs = HYDROGEN_MASSFRAC * dens_cgs / PROTONMASS;

#ifdef USE_GRACKLE
  *temp_K = CallGrackle(SphP[i].Utherm, dens_code, 0.0, i, 2);
#else
  double ne    = SphP[i].Ne;
  double u_cgs = SphP[i].Utherm * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
  *temp_K      = convert_u_to_temp(u_cgs, dens_cgs, &ne);
#endif
}

/*! \brief Residual for the implicit backward-Euler dust growth/sputtering
 *  update: M - M_old - dt * (growth_rate(M) + sputtering_rate(M)).
 */
static double dust_ode_residual(double M, double M_old, double M_metal, double n_H_cgs, double temp_K,
                                 double Z_gas_massfrac, double dt_gyr)
{
  double rate = dust_growth_rate(DUST_PHASE1_SPECIES, M, M_metal, n_H_cgs, temp_K, Z_gas_massfrac) +
                dust_sputtering_rate(DUST_PHASE1_SPECIES, M, n_H_cgs, temp_K);

  return M - M_old - dt_gyr * rate;
}

/*! \brief Apply the combined growth + thermal sputtering rate-ODE to a
 *  single gas cell's dust mass, over its own hydro timestep.
 *
 *  Solved via bisection, bracketed by the cell's physical dust-mass bounds
 *  [0, GasMetals] (dust can never exceed the metal budget it condensed
 *  from) -- the same bracket+bisect spirit as DoCooling() in
 *  src/cooling/cooling.c, but bounded by physics rather than an expanding
 *  search, since M_dust (unlike internal energy) has known finite bounds.
 *
 *  \param[in] i index of the gas cell
 */
void dust_cell(int i)
{
  double dt, dtime, dt_gyr;

  dt    = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
  dtime = All.cf_atime * dt / All.cf_time_hubble_a;

  if(dtime <= 0)
    return;

  dt_gyr = (dtime * All.UnitTime_in_s / All.HubbleParam) / SEC_PER_GIGAYEAR;

  double M_old   = SphP[i].GasDustMass;
  double M_metal = SphP[i].GasMetals;

  if(M_metal <= 0)
    {
      SphP[i].GasDustMass = 0;
      sync_primitive_from_conserved(i, DUST_INDEX);
      return;
    }

  double n_H_cgs, temp_K;
  dust_get_nH_and_temp(i, &n_H_cgs, &temp_K);

  double Z_gas_massfrac = SphP[i].GasMetallicity;

  double lo = 0.0, hi = M_metal;
  double f_lo = dust_ode_residual(lo, M_old, M_metal, n_H_cgs, temp_K, Z_gas_massfrac, dt_gyr);
  double f_hi = dust_ode_residual(hi, M_old, M_metal, n_H_cgs, temp_K, Z_gas_massfrac, dt_gyr);

  double M_new;

  if(f_lo >= 0)
    {
      M_new = lo; /* even M=0 overshoots: sputtering removes all dust this step */
    }
  else if(f_hi <= 0)
    {
      M_new = hi; /* even M=M_metal undershoots: growth saturates at the metal budget */
    }
  else
    {
      int iter = 0;
      double dM;

      do
        {
          M_new     = 0.5 * (lo + hi);
          double f  = dust_ode_residual(M_new, M_old, M_metal, n_H_cgs, temp_K, Z_gas_massfrac, dt_gyr);

          if(f > 0)
            hi = M_new;
          else
            lo = M_new;

          dM = hi - lo;
          iter++;
        }
      while(dM > DUST_ODE_TOLERANCE * fmax(hi, 1e-300) && iter < DUST_ODE_MAXITER);

      M_new = 0.5 * (lo + hi);
    }

  if(M_new < 0)
    M_new = 0;
  if(M_new > M_metal)
    M_new = M_metal;

  SphP[i].GasDustMass = M_new;
  sync_primitive_from_conserved(i, DUST_INDEX);
}

/*! \brief Apply dust growth/sputtering to all active gas cells.
 *
 *  Mirrors cooling_only() in src/cooling/cooling.c, and is called
 *  immediately after the COOLING block in
 *  calculate_non_standard_physics_end_of_step() (src/main/run.c) so the
 *  rate laws see post-cooling temperature.
 */
void dust_processes(void)
{
  int idx, i;

  CPU_Step[CPU_MISC] += measure_time();

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i >= 0)
        {
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue; /* skip cells that have been swallowed or eliminated */

          dust_cell(i);
        }
    }

  CPU_Step[CPU_COOLINGSFR] += measure_time();
}

#endif /* DUST */
