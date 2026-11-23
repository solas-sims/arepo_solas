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
/* Jeans energy floor / stiffened equation of state for gas whose Jeans mass is not resolved.
 *
 * Mode 0 and 1 (Robertson & Kravtsov 2008, eqs 7-10 and A2; Truelove et al. 1997 in grid codes): a cell whose Jeans mass is below
 * JEANS_FLOOR_NCELLS cell masses has its specific internal energy raised until the Jeans mass reaches that value. Since
 * m_Jeans = pi^(5/2) c_s^3 / (6 G^(3/2) rho^(1/2)) and c_s^2 = gamma (gamma-1) u, this gives
 *     u_floor = (6 N M_cell / pi^(5/2))^(2/3) G rho^(1/3) / (gamma (gamma-1)),
 * which scales as rho^(1/3) at fixed cell mass, so the floor temperature scales as M_cell^(2/3): at higher resolution, with the same N,
 * it falls below the physical gas temperatures and the floor removes itself. All quantities are physical, in code units.
 *
 * Mode 2 (Schaye & Dalla Vecchia 2008, as used in EAGLE): a fixed polytropic floor T = T0 (n_H/n0)^(gamma_eos-1) above n0, independent
 * of the cell mass.
 *
 * Options (compile time):
 *   JEANS_FLOOR_NCELLS   required Jeans mass in cell masses (default 33: Truelove et al.'s lambda_J >= 4 dx; the SPH criterion of
 *                        Robertson & Kravtsov, 2 N_Jeans = 30 kernel masses, is roughly 420 if a kernel is the cell plus ~13 neighbours)
 *   JEANS_FLOOR_MODE     0: raise the energy after the cooling step
 *                        1 (default): as 0, and in addition cells already at the floor are not passed to the cooling solver, so cooling
 *                           can never remove energy from gas that sits on the floor (and the solver is not called for it)
 *                        2: fixed EAGLE-style equation of state (JEANS_FLOOR_EOS_GAMMA, _NH0, _T0), applied as in mode 1
 *   JEANS_FLOOR_NH_MIN   [cm^-3] no floor below this hydrogen number density (default 0: none)
 *   JEANS_FLOOR_SF_GATE  0 (default): the star formation temperature gate sees the floored energy
 *                        1: cells at the floor pass the temperature gate
 *                        2: gas whose energy is within JEANS_FLOOR_SF_DEX dex of the floor passes the gate (EAGLE-like) */
#ifndef JEANS_FLOOR_NCELLS
#define JEANS_FLOOR_NCELLS 33.0
#endif
#ifndef JEANS_FLOOR_MODE
#define JEANS_FLOOR_MODE 1
#endif
#ifndef JEANS_FLOOR_NH_MIN
#define JEANS_FLOOR_NH_MIN 0.0
#endif
#ifndef JEANS_FLOOR_SF_GATE
#define JEANS_FLOOR_SF_GATE 0
#endif
#ifndef JEANS_FLOOR_SF_DEX
#define JEANS_FLOOR_SF_DEX 0.5
#endif
#ifndef JEANS_FLOOR_EOS_GAMMA
#define JEANS_FLOOR_EOS_GAMMA (4.0 / 3.0)
#endif
#ifndef JEANS_FLOOR_EOS_NH0
#define JEANS_FLOOR_EOS_NH0 0.1
#endif
#ifndef JEANS_FLOOR_EOS_T0
#define JEANS_FLOOR_EOS_T0 8000.0
#endif

static double JFloorCalls = 0, JFloorActive = 0, JFloorHit = 0, JFloorRatio = 0, JFloorEnergy = 0, JFloorSkipped = 0;

/* specific internal energy floor of gas cell i (code units); 0 where no floor applies */
static double jeans_floor_u(int i)
{
  double rho_phys = SphP[i].Density * All.cf_a3inv;                                  /* code units */
  double nH       = HYDROGEN_MASSFRAC * SphP[i].Density * All.cf_UnitDensity_in_cgs / PROTONMASS; /* cm^-3, includes the h^2/a^3 factor */

  if(nH < JEANS_FLOOR_NH_MIN)
    return 0;

#if JEANS_FLOOR_MODE == 2
  if(nH < JEANS_FLOOR_EOS_NH0)
    return 0;

  double T = JEANS_FLOOR_EOS_T0 * pow(nH / JEANS_FLOOR_EOS_NH0, JEANS_FLOOR_EOS_GAMMA - 1.0);
#ifdef USE_GRACKLE
  double mu = compute_mu(i);
#else
  double mu = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));
#endif
  /* inverse of the temperature definition used in sf_criteria() */
  return T * BOLTZMANN / (mu * PROTONMASS * GAMMA_MINUS1) / (All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);
#else
  const double coef = pow(6.0 * JEANS_FLOOR_NCELLS / pow(M_PI, 2.5), 2.0 / 3.0) * All.G / (GAMMA * GAMMA_MINUS1);
  return coef * pow(P[i].Mass, 2.0 / 3.0) * cbrt(rho_phys);
#endif
}

/* Cool the active cells; in mode 1 and 2 the cells that already sit on the floor are not passed to the cooling solver. */
static void jeans_floor_cool_active_cells(void)
{
#if JEANS_FLOOR_MODE >= 1
  int *list = (int *)mymalloc("jeans_floor_list", (TimeBinsHydro.NActiveParticles + 1) * sizeof(int));
  int n = 0;

  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0 || P[i].Mass == 0)
        continue;

      double uf = jeans_floor_u(i);
      if(uf > 0 && SphP[i].Utherm <= uf)
        JFloorSkipped += 1;
      else
        list[n++] = i;
    }

  cool_cells(n, list);
  myfree(list);
#else
  cool_active_cells();
#endif
}

static void apply_jeans_energy_floor(void)
{
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0 || P[i].Mass == 0)
        continue;

      double u_floor = jeans_floor_u(i);

      JFloorActive += 1;
      SphP[i].JeansFloored = 0;

      if(u_floor > 0 && SphP[i].Utherm <= u_floor)
        {
          double du = u_floor - SphP[i].Utherm;

          if(du > 0)
            {
              JFloorHit += 1;
              JFloorRatio += u_floor / fmax(SphP[i].Utherm, 1e-300);
              JFloorEnergy += du * P[i].Mass;

              SphP[i].Utherm += du;
              SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;
              set_pressure_of_cell(i);
            }
          SphP[i].JeansFloored = 1;
        }
    }

  /* report every 200 calls: fraction of active cells raised, mean u_floor/u before, energy added after cooling (in erg), and the number of
   * cells that were kept out of the cooling solver */
  if(++JFloorCalls >= 200)
    {
      double loc[5] = {JFloorActive, JFloorHit, JFloorRatio, JFloorEnergy, JFloorSkipped}, glob[5];
      MPI_Allreduce(loc, glob, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      mpi_printf("JEANS_FLOOR: a=%g active_cells=%.0f raised=%.2f%% skipped_cooling=%.2f%% mean(u_floor/u)=%.3g energy_added=%.4g erg (over last %d steps)\n",
                 All.Time, glob[0], 100.0 * glob[1] / fmax(glob[0], 1.0), 100.0 * glob[4] / fmax(glob[0], 1.0), glob[2] / fmax(glob[1], 1.0),
                 glob[3] * All.cf_UnitEnergy_in_cgs, (int)JFloorCalls);
      JFloorCalls = JFloorActive = JFloorHit = JFloorRatio = JFloorEnergy = JFloorSkipped = 0;
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
#if defined(JEANS_ENERGY_FLOOR) && JEANS_FLOOR_SF_GATE == 1
      /* the heat of a cell at the Jeans floor is artificial pressure, not gas that is too hot to form stars */
      if(!SphP[i].JeansFloored)
        return 0;
#elif defined(JEANS_ENERGY_FLOOR) && JEANS_FLOOR_SF_GATE == 2
      /* EAGLE-like: gas within JEANS_FLOOR_SF_DEX dex of the floor is treated as star-forming gas on the equation of state */
      {
        double uf = jeans_floor_u(i);
        if(!(uf > 0 && SphP[i].Utherm <= uf * pow(10.0, JEANS_FLOOR_SF_DEX)))
          return 0;
      }
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
#ifdef JEANS_ENERGY_FLOOR
  jeans_floor_cool_active_cells();

  /* applied after cooling so that the floor is the last word on the internal energy of each active cell */
  apply_jeans_energy_floor();
#else
  cool_active_cells();
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
