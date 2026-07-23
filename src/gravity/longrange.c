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
 * \file        src/gravity/longrange.c
 * \date        05/2018
 * \brief       Driver routines for computation of long-range gravitational
 *              PM force
 * \details     contains functions:
 *                void long_range_init(void)
 *                void long_range_init_regionsize(void)
 *                void long_range_force(void)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 06.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef PMGRID
/*! \brief Driver routine to call initialization of periodic or/and
 *         non-periodic FFT routines.
 *
 *  \return void
 */
void long_range_init(void)
{
#ifndef GRAVITY_NOT_PERIODIC
  pm_init_periodic();
#ifdef TWODIMS
  pm2d_init_periodic();
#endif /* #ifdef TWODIMS */
#ifdef PLACEHIGHRESREGION
  pm_init_nonperiodic();
#endif /* #ifdef PLACEHIGHRESREGION */
#else  /* #ifndef GRAVITY_NOT_PERIODIC */
  pm_init_nonperiodic();
#endif /* #ifndef GRAVITY_NOT_PERIODIC #else */
}

/*! \brief Driver routine to determine the extend of the non-
 *         periodic or high resolution region.
 *
 *  The initialization is done by pm_init_regionsize(). Afterwards
 *  the convolution kernels are computed by pm_setup_nonperiodic_kernel().
 *
 *  \return void
 */
void long_range_init_regionsize(void)
{
#ifndef GRAVITY_NOT_PERIODIC
#ifdef PLACEHIGHRESREGION
  if(RestartFlag != 1)
    pm_init_regionsize();
  pm_setup_nonperiodic_kernel();
#endif /* #ifdef PLACEHIGHRESREGION */

#else  /* #ifndef GRAVITY_NOT_PERIODIC */
  if(RestartFlag != 1)
    pm_init_regionsize();
  pm_setup_nonperiodic_kernel();
#endif /* #ifndef GRAVITY_NOT_PERIODIC #else */
}

/*! \brief This function computes the long-range PM force for all particles.
 *
 *  In case of a periodic grid the force is calculated by pmforce_periodic()
 *  otherwise by pmforce_nonperiodic(). If a high resolution region is
 *  specified for the PM force, pmforce_nonperiodic() calculates that force in
 *  both cases.
 *
 *  \return void
 */
void long_range_force(void)
{
  int i;

#ifdef GRAVITY_NOT_PERIODIC
  int j;
  double fac;
#endif /* #ifdef GRAVITY_NOT_PERIODIC */

  for(i = 0; i < NumPart; i++)
    {
      P[i].GravPM[0] = P[i].GravPM[1] = P[i].GravPM[2] = 0;
#ifdef EVALPOTENTIAL
      P[i].PM_Potential = 0;
#endif /* #ifdef EVALPOTENTIAL */
    }

  /* Moved here (before the SELFGRAVITY early-return just below), not
   * left at this function's own end where it originally lived -- a
   * real, confirmed bug: this function is the ONLY place
   * All.DtDisplacement ever gets computed, and get_timestep_gravity()
   * (timestep.c) unconditionally clamps every particle's timestep to
   * All.DtDisplacement whenever PMGRID is compiled in, regardless of
   * SELFGRAVITY. FDM requires PMGRID for its own FFT infrastructure
   * without necessarily needing SELFGRAVITY's own PM force computation
   * -- exactly this project's own SELFGRAVITY-off configuration -- so
   * leaving this call after the early return meant All.DtDisplacement
   * stayed at its zero-initialized default forever, clamping every
   * timestep to zero regardless of how good the actual force
   * computation was (confirmed directly: GravAccel was genuinely
   * correct and substantial, dt was still 0). Calling it here instead
   * is safe with SELFGRAVITY off too: find_long_range_step_constraint()
   * computes its own acceleration from GravPM directly (zeroed just
   * above, for every particle, unconditionally), floors that at
   * MIN_FLOAT_NUMBER, and the resulting (very large) dt gets clamped
   * to All.MaxSizeTimestep by this same function's own bounds-checking
   * -- exactly the correct behaviour when there's no PM force to
   * constrain the displacement timestep by. */
  find_long_range_step_constraint();

#ifndef SELFGRAVITY
  return;
#endif /* #ifndef SELFGRAVITY */

  /* TIMER_START moved to here (after the SELFGRAVITY early-return
   * check above), not at this function's own top -- without this,
   * the timer starts unconditionally, but the early return above skips
   * the matching TIMER_STOP at this function's own end whenever
   * SELFGRAVITY isn't compiled in, leaving CPU_PM_GRAVITY stuck
   * "running" on the timer stack -- the NEXT call to this function
   * (there are two calls per sync point: once from the one-time
   * startup bootstrap, once from the main loop) then fails its own
   * TIMER_START with "already running". PMGRID has historically only
   * ever been used alongside SELFGRAVITY (its own tree-PM force
   * split); FDM also requires PMGRID (for the FFT/fft_plan
   * infrastructure it reuses) without necessarily needing
   * SELFGRAVITY's own PM force computation at all -- exactly this
   * project's own SELFGRAVITY-off configuration. */
  TIMER_START(CPU_PM_GRAVITY);

#ifndef GRAVITY_NOT_PERIODIC

#ifdef TWODIMS
  pm2d_force_periodic(0);
#else  /* #ifdef TWODIMS */
  pmforce_periodic(0, NULL);
#endif /* #ifdef TWODIMS #else */

#ifdef PLACEHIGHRESREGION
  i = pmforce_nonperiodic(1);

  if(i == 1) /* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmforce_nonperiodic(1); /* try again */
    }
  if(i == 1)
    terminate("despite we tried to increase the region, we still don't fit all particles in it");
#endif /* #ifdef PLACEHIGHRESREGION */

#else /* #ifndef GRAVITY_NOT_PERIODIC */
  i = pmforce_nonperiodic(0);

  if(i == 1) /* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmforce_nonperiodic(0); /* try again */
    }
  if(i == 1)
    terminate("despite we tried to increase the region, somehow we still don't fit all particles in it");
#ifdef PLACEHIGHRESREGION
  i = pmforce_nonperiodic(1);

  if(i == 1) /* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();

      /* try again */

      for(i = 0; i < NumPart; i++)
        P[i].GravPM[0] = P[i].GravPM[1] = P[i].GravPM[2] = 0;

      i = pmforce_nonperiodic(0) + pmforce_nonperiodic(1);
    }
  if(i != 0)
    terminate("despite we tried to increase the region, somehow we still don't fit all particles in it");
#endif /* #ifdef PLACEHIGHRESREGION */
#endif /* #ifndef GRAVITY_NOT_PERIODIC #else */

#ifdef GRAVITY_NOT_PERIODIC
  if(All.ComovingIntegrationOn)
    {
      fac = 0.5 * All.Hubble * All.Hubble * All.Omega0;

      for(i = 0; i < NumPart; i++)
        for(j = 0; j < 3; j++)
          P[i].GravPM[j] += fac * P[i].Pos[j];
    }

  /* Finally, the following factor allows a computation of cosmological simulation
     with vacuum energy in physical coordinates */
  if(All.ComovingIntegrationOn == 0)
    {
      fac = All.OmegaLambda * All.Hubble * All.Hubble;

      for(i = 0; i < NumPart; i++)
        for(j = 0; j < 3; j++)
          P[i].GravPM[j] += fac * P[i].Pos[j];
    }
#endif /* #ifdef GRAVITY_NOT_PERIODIC */

  TIMER_STOP(CPU_PM_GRAVITY);
}
#endif /* #ifdef PMGRID */
