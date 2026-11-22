/*!
 * \copyright   This file is part of the public version of the AREPO code.
 * \copyright   Copyright (C) 2009-2019, Max-Planck Institute for Astrophysics
 * \copyright   Developed by Volker Springel (vspringel@MPA-Garching.MPG.DE) and
 *              contributing authors.
 * \copyright   Arepo is free software: you can redistribute it and/or modify
 *              it under the terms of the GNU General Public License as published by
 *              the Free Software Foundation, either version 3 of the License, or
 *              (at your option) any later version.
 *
 *              Arepo is distributed in the hope that it will be useful,
 *              but WITHOUT ANY WARRANTY; without even the implied warranty of
 *              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *              GNU General Public License for more details.
 *
 *              A copy of the GNU General Public License is available under
 *              LICENSE as part of this program.  See also
 *              <https://www.gnu.org/licenses/>.
 *
 * \file        src/star_formation/sfr_AGORA.c
 * \date        05/2018
 * \brief       Star formation rate routines for the effective multi-phase
 *              model.
 * \details     contains functions:
 *                void cooling_and_starformation(void)
 *                double get_starformation_rate(int i)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../gravity/forcetree.h"

#ifdef AGORA_SF

#ifdef JEANS_ENERGY_FLOOR
/* Jeans energy floor (Robertson & Kravtsov 2008, their eqs 7-10 and A2; Truelove et al. 1997 in grid codes).
 *
 * A cell whose Jeans mass is below JEANS_FLOOR_NCELLS cell masses has its specific internal energy raised until the Jeans mass reaches that
 * value. Since m_Jeans = pi^(5/2) c_s^3 / (6 G^(3/2) rho^(1/2)) and c_s^2 = gamma (gamma-1) u, this gives
 *     u_floor = (6 N M_cell / pi^(5/2))^(2/3) G rho^(1/3) / (gamma (gamma-1)),
 * which scales as rho^(1/3) at fixed cell mass. All quantities are physical, in code units.
 *
 * JEANS_FLOOR_NCELLS: required Jeans mass in units of the cell mass (default 33: Truelove et al.'s lambda_J >= 4 dx; the SPH criterion
 *   of Robertson & Kravtsov, 2 N_Jeans = 30 kernel masses, corresponds to roughly 420 if a kernel is the cell plus ~13 neighbours).
 * JEANS_FLOOR_SF_GATE: 0 (default) the temperature gate of the star formation criteria sees the floored energy;
 *   1 cells whose energy was raised by the floor pass the temperature gate (the gas would be cold without the artificial pressure). */
#ifndef JEANS_FLOOR_NCELLS
#define JEANS_FLOOR_NCELLS 33.0
#endif
#ifndef JEANS_FLOOR_SF_GATE
#define JEANS_FLOOR_SF_GATE 0
#endif

static double JFloorCalls = 0, JFloorActive = 0, JFloorHit = 0, JFloorRatio = 0, JFloorEnergy = 0, JFloorEnergyTot = 0;

static void apply_jeans_energy_floor(void)
{
  const double coef = pow(6.0 * JEANS_FLOOR_NCELLS / pow(M_PI, 2.5), 2.0 / 3.0) * All.G / (GAMMA * GAMMA_MINUS1);

  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0 || P[i].Mass == 0)
        continue;

      double rho_phys = SphP[i].Density * All.cf_a3inv;
      double u_floor  = coef * pow(P[i].Mass, 2.0 / 3.0) * cbrt(rho_phys);

      JFloorActive += 1;
      SphP[i].JeansFloored = 0;

      if(SphP[i].Utherm < u_floor)
        {
          double du = u_floor - SphP[i].Utherm;

          JFloorHit += 1;
          JFloorRatio += u_floor / fmax(SphP[i].Utherm, 1e-300);
          JFloorEnergy += du * P[i].Mass;
          JFloorEnergyTot += SphP[i].Energy;

          SphP[i].Utherm += du;
          SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;
          SphP[i].JeansFloored = 1;
          set_pressure_of_cell(i);
        }
    }

  /* report every 200 calls: fraction of active cells floored, mean u_floor/u before, energy added in erg */
  if(++JFloorCalls >= 200)
    {
      double loc[4] = {JFloorActive, JFloorHit, JFloorRatio, JFloorEnergy}, glob[4];
      MPI_Allreduce(loc, glob, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      mpi_printf("JEANS_FLOOR: a=%g active_cells=%.0f floored=%.2f%% mean(u_floor/u)=%.3g energy_added=%.4g erg (over last %d steps)\n", All.Time,
                 glob[0], 100.0 * glob[1] / fmax(glob[0], 1.0), glob[2] / fmax(glob[1], 1.0),
                 glob[3] * All.cf_UnitEnergy_in_cgs, (int)JFloorCalls);
      JFloorCalls = JFloorActive = JFloorHit = JFloorRatio = JFloorEnergy = JFloorEnergyTot = 0;
    }
}
#endif /* #ifdef JEANS_ENERGY_FLOOR */

/* Function that checks whether a cell i satisfies star formation criteria*/
static int sf_criteria(int i)
{
  /* comoving-overdensity floor, shared with EEOS_SF/JEANS_SF (see set_overdens_thresh() in
   * starformation.c) -- without this, diffuse, cosmologically-unvirialized gas can
   * transiently satisfy the local density/temperature gate below well before any real
   * structure has collapsed, especially with less conservative threshold parameters */
  if(All.ComovingIntegrationOn)
    if(SphP[i].Density < All.OverDensThresh)
      return 0;

#ifdef USE_GRACKLE
#if defined(POPIII_SF) && (GRACKLE_CHEMISTRY >= 2)
  if(SphP[i].GasMetallicity < All.PopIIIMetallicityThreshold &&
     get_H2_fraction(i) < All.PopIIIH2FractionThreshold)
    return 0;
#endif
  double mu = compute_mu(i);
#else // To be replaced 
  double mu = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));
#endif

  double number_dens = (SphP[i].Density * All.cf_UnitDensity_in_cgs) / mu / PROTONMASS;
  double temp = (SphP[i].Utherm * All.cf_UnitVelocity_in_cm_per_s * All.cf_UnitVelocity_in_cm_per_s)
  * mu * PROTONMASS * GAMMA_MINUS1 / BOLTZMANN;

  double number_dens_threshold = All.NumberDensThreshold;
#ifdef SF_THRESHOLD_HALO_MASS_DEPENDENT
  number_dens_threshold *= sf_threshold_halo_mass_factor(i);
#endif /* #ifdef SF_THRESHOLD_HALO_MASS_DEPENDENT */

  if(number_dens < number_dens_threshold)
    return 0;
  
  if(temp > All.TemperatureThreshold)
    {
#if defined(JEANS_ENERGY_FLOOR) && JEANS_FLOOR_SF_GATE
      /* the heat of a cell raised by the Jeans floor is artificial pressure, not gas that is too hot to form stars */
      if(!SphP[i].JeansFloored)
        return 0;
#else
      return 0;
#endif
    }

#ifdef DIVVEL
  if(SphP[i].DivVel >= 0)
    return 0;
#endif

  return 1;
}

/*! \brief Main driver for star formation and gas cooling.
 *
 *  This function loops over all the active gas cells. If a given cell
 *  meets the criteria for star formation to be active the multi-phase
 *  model is activated, the properties of the cell are updated according to
 *  the latter and the star formation rate computed. In the other case, the
 *  standard isochoric cooling is applied to the gas cell by calling the
 *  function cool_cell() and the star formation rate is set to 0.
 *
 *  \return void
 */
void cooling_and_starformation(void)
{
  TIMER_START(CPU_COOLINGSFR);

  int idx, i, bin;
  double unew, du;
    
  /* clear the SFR stored in the active timebins */
  for(bin = 0; bin < TIMEBINS; bin++)
    if(TimeBinSynchronized[bin])
      TimeBinSfr[bin] = 0;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      
      /* skip cells that have been swallowed or eliminated */
      if(P[i].Mass == 0)
        continue; 
       
      /* apply the temperature floor */
      unew = dmax(All.MinEgySpec, SphP[i].Utherm);

      if(unew < 0)
        terminate("Invalid Temperature: Task=%d i=%d unew=%g\n", ThisTask, i, unew);

      du = unew - SphP[i].Utherm;
      SphP[i].Utherm += du;
      SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;

    } /* end of temperature-floor loop */

  /* cool all active cells in one go (with Grackle: array-based calls, one per time bin) */
  cool_active_cells();

#ifdef JEANS_ENERGY_FLOOR
  /* applied after cooling so that the floor is the last word on the internal energy of each active cell */
  apply_jeans_energy_floor();
#endif

  /* sf_criteria() only depends on the state of the cell itself, so it can be evaluated after the cooling pass */
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0 || P[i].Mass == 0)
        continue;

      if(sf_criteria(i))
        {
          SphP[i].Sfr = get_starformation_rate(i);
          TimeBinSfr[P[i].TimeBinHydro] += SphP[i].Sfr;
        }
      else
        SphP[i].Sfr = 0;

    } /* end of star formation loop over active particles */

  TIMER_STOP(CPU_COOLINGSFR);
}

/*! \brief Return the star formation rate associated with the gas cell i.
 *
 *  \param[in] i the index of the gas cell.
 *
 *  \return star formation rate in solar masses / yr.
 */
double get_starformation_rate(int i)
{
  if(RestartFlag == 3)
    return SphP[i].Sfr;

  double rateOfSF, t_freefall;

  if(!sf_criteria(i))
    return 0;
    
  /* freefall time in code units; Density is stored comoving, so convert to physical
   * density (* cf_a3inv) before using it in a physical free-fall time */
  t_freefall = sqrt(3. * M_PI / 32 / All.G / (SphP[i].Density * All.cf_a3inv));

  rateOfSF = All.StarFormationEfficiency * P[i].Mass / t_freefall;

  /* convert to solar masses per yr */
  rateOfSF *= (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

  return rateOfSF;
}

#endif /* #ifdef AGORA_SF */
