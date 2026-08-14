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

/* Function that checks whether a cell i satisfies star formation criteria*/
static int sf_criteria(int i)
{
#ifdef USE_GRACKLE
#if defined(POPIII_SF) && (GRACKLE_CHEMISTRY >= 2)
  if(SphP[i].GasMetallicity < All.PopIIIMetallicityThreshold &&
     get_H2_fraction(i) < All.PopIIIH2FractionThreshold)
    return 0;
#endif
  double mu = compute_mu(i);
#else
  double mu = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));
#endif

  double number_dens = (SphP[i].Density * All.cf_UnitDensity_in_cgs) / mu / PROTONMASS;
  double temp = (SphP[i].Utherm * All.cf_UnitVelocity_in_cm_per_s * All.cf_UnitVelocity_in_cm_per_s) 
  * mu * PROTONMASS * GAMMA_MINUS1 / BOLTZMANN;

  if(number_dens < All.NumberDensThreshold)
    return 0;
  
  if(temp > All.TemperatureThreshold)
    return 0;

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

      cool_cell(i);

      if(sf_criteria(i))
        {
          SphP[i].Sfr = get_starformation_rate(i);
          TimeBinSfr[P[i].TimeBinHydro] += SphP[i].Sfr;
        }
      else
        SphP[i].Sfr = 0;

    } /* end of main loop over active particles */

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
    
  /* freefall time in code units */
  t_freefall = sqrt(3. * M_PI / 32 / All.G / SphP[i].Density); 

  rateOfSF = All.StarFormationEfficiency * P[i].Mass / t_freefall;

  /* convert to solar masses per yr */
  rateOfSF *= (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

  return rateOfSF;
}