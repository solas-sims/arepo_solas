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
 * \file        src/gravity/accel.c
 * \date        05/2018
 * \brief       Routines to carry out gravity force computation.
 * \details     contains functions:
 *                void compute_grav_accelerations(int timebin, int fullflag)
 *                void gravity(int timebin, int fullflag)
 *                void gravity_force_finalize(int timebin)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 03.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"
#ifdef FDM
#include "../fdm/fdm.h"
#endif /* #ifdef FDM */

/*! \brief Computes the gravitational accelerations for all active particles.
 *
 *  If the particle mesh is used and the current time step
 *  requires a PM force computation, new long range forces are
 *  computed by long_range_force(). Then the short-range tree forces
 *  are computed by gravity(). The force tree is rebuild every time step.
 *
 *  \param[in] timebin Current timebin for which gravity is calculated
 *             (positive integer).
 *  \param[in] fullflag Flag whether this is a global timestep
 *             (Flag_Full_Tree, Flag_Partial_Tree).
 *
 *  \return void
 */
void compute_grav_accelerations(int timebin, int fullflag)
{
  if(TimeBinsGravity.GlobalNActiveParticles > 0)
    {
      if(All.TypeOfOpeningCriterion == 1 && All.Ti_Current == 0 && All.ErrTolTheta > 0)
        {
          /* For the first timestep, we do one gravity calculation up front
           * with the Barnes & Hut Criterion to allow usage of relative opening
           * criterion with consistent accuracy.
           */
#ifdef PMGRID
          long_range_force();
#endif /* #ifdef PMGRID */
          gravity(timebin, fullflag);
        }

      gravity(timebin, fullflag); /* computes (short-range) gravity accel. */

#ifdef FORCETEST
      gravity_forcetest();
#endif /* #ifdef FORCETEST */
    }
}

/*! \brief Main routine for tree force calculation.
 *
 *  This routine handles the tree force calculation. First it builds a new
 *  force tree calling force_treebuild() at every timestep. This tree is then
 *  used to calculate a new tree force for every active particle by calling
 *  gravity_tree().
 *
 *  \param[in] timebin Current timebin for which gravity is calculated.
 *  \param[in] fullflag Flag whether this is a global timestep.
 *
 *  \return void
 */
#ifdef FDM
/*! \brief Adds FDM_StarResult's force to P[i].GravAccel for every
 *  ACTIVE star (Type==4) particle -- factored out as its own function
 *  (rather than left inline inside gravity_force_finalize(), where it
 *  originally lived) purely for readability.
 *
 *  Called from TWO places: gravity_force_finalize() itself (the
 *  SELFGRAVITY-on path), and gravity()'s own "self-gravity is switched
 *  off" branch. The second call site was added, then deliberately
 *  reverted earlier in this project (as speculative, unvalidated
 *  scaffolding for a SELFGRAVITY-off configuration that wasn't fully
 *  working at the time), then reinstated once SELFGRAVITY-off became
 *  genuinely needed -- as a diagnostic, to isolate FDM-only dynamics
 *  from stellar self-gravity when investigating whether some observed
 *  behaviour is real physics or an artifact of the two being coupled
 *  together. Without this call in the SELFGRAVITY-off branch, stars
 *  get literally zero gravity of any kind there (no self-gravity,
 *  since that's what's off; no FDM force either, since this function
 *  is otherwise only reached via gravity_force_finalize(), which that
 *  branch never calls) -- producing an unconstrained, meaningless
 *  single timestep spanning the entire run, exactly the symptom that
 *  led back here.
 *
 *  FDM_StarResult's own force is -grad(Phi) where Phi is a potential
 *  per unit mass (this project's consistent convention throughout) --
 *  already an acceleration, already including G (baked into Phi's own
 *  normalization in fdm_poisson.c) -- must be added AFTER any *=All.G
 *  step (SELFGRAVITY-on path) or the zeroing loop (SELFGRAVITY-off
 *  path), not before, or G would be double-applied to this term
 *  specifically.
 *
 *  Requires FDM_StarResult[] to have already been populated for the
 *  current sync point by fdm_interpolate_to_stars() (fdm_advance_to_time(),
 *  fdm_integrator.c) -- see run.c's own two call sites for that function:
 *  the one-time startup bootstrap sequence AND the main loop both call
 *  do_gravity_step_second_half() (which reaches this function), and both
 *  needed their own fdm_advance_to_time() call immediately before it.
 *  Missing the bootstrap one specifically left FDM_StarResult[] as
 *  genuinely uninitialized memory the first time this function ever
 *  ran -- found from a real crash (an enormous, nonsensical FDM force on
 *  the very first sync point) via instrumentation showing the
 *  "populate" side had never executed at all before the "consume" side
 *  read it, not a staleness/reordering issue as originally suspected.
 */
static void fdm_add_force_to_active_stars(void)
{
  for(int idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      int i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type == 4)
        {
          P[i].GravAccel[0] += FDM_StarResult[i].ForceX;
          P[i].GravAccel[1] += FDM_StarResult[i].ForceY;
          P[i].GravAccel[2] += FDM_StarResult[i].ForceZ;
        }
    }
}
#endif /* #ifdef FDM */

void gravity(int timebin, int fullflag)
{
  double tstart = second();

#if defined(SELFGRAVITY)
  /* set new softening lengths on global steps to take into account possible cosmological time variation */
  if(timebin == All.HighestOccupiedGravTimeBin)
    set_softenings();

#ifdef ALLOW_DIRECT_SUMMATION
  if(TimeBinsGravity.GlobalNActiveParticles < DIRECT_SUMMATION_THRESHOLD)
    {
      gravity_direct(timebin);

#ifndef ONEDIMS_SPHERICAL
      gravity_force_finalize(timebin);
#endif /* #ifndef ONEDIMS_SPHERICAL */

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
      calc_exact_gravity_for_particle_type();
#endif /* #ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE */

#ifdef EXTERNALGRAVITY
      gravity_external();
#endif /* #ifdef EXTERNALGRAVITY */
    }
  else
#endif /* #ifdef ALLOW_DIRECT_SUMMATION */
    {
#ifdef ONEDIMS_SPHERICAL
      gravity_monopole_1d_spherical();
#else /* #ifdef ONEDIMS_SPHERICAL */

    if(TimeBinsGravity.GlobalNActiveParticles >= 10 * NTask)
      construct_forcetree(0, 1, 0, timebin); /* build force tree with all particles */
    else
      construct_forcetree(0, 0, 0, timebin); /* build force tree with all particles */

    gravity_tree(timebin);

    gravity_force_finalize(timebin);

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
    calc_exact_gravity_for_particle_type();
#endif /* #ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE */

#ifdef EXTERNALGRAVITY
    gravity_external();
#endif /* #ifdef EXTERNALGRAVITY */

    /* note: we here moved 'gravity_force_finalize' in front of the non-standard physics;
     * reminder: restart flag 18: post-processing calculation potential without running simulation
     */
    if(fullflag == FLAG_FULL_TREE && RestartFlag != 18)
      calculate_non_standard_physics_with_valid_gravity_tree();

    /* this is for runs which have the full tree at each time step; no HIERARCHICAL_GRAVITY */
    calculate_non_standard_physics_with_valid_gravity_tree_always();

    myfree(Father);
    myfree(Nextnode);
    myfree(Tree_Points);
    force_treefree();
#endif /* #ifdef ONEDIMS_SPHERICAL #else */
    }

#else /* defined(SELFGRAVITY) */

  /* self-gravity is switched off */
  int idx, i, j;
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];

      if(i < 0)
        continue;

#ifdef EVALPOTENTIAL
      P[i].Potential = 0;
#endif /* #ifdef EVALPOTENTIAL */

      for(j = 0; j < 3; j++)
        P[i].GravAccel[j] = 0;
    }

#ifdef FDM
  fdm_add_force_to_active_stars(); /* without this here, SELFGRAVITY-off runs give stars literally zero
                                     * gravity of any kind (no self-gravity, since that's what's off; no
                                     * FDM force either, since gravity_force_finalize() -- where this call
                                     * normally lives -- is never reached from this branch), producing an
                                     * unconstrained, meaningless single timestep spanning the whole run. */
#endif /* #ifdef FDM */

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
  calc_exact_gravity_for_particle_type();
#endif /* #ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE */

#ifdef EXTERNALGRAVITY
  gravity_external();
#endif /* #ifdef EXTERNALGRAVITY */

#endif /* defined(SELFGRAVITY) #else */

  double tend = second();
  mpi_printf("GRAVITY: done for timebin %d,  %lld particles  (took %g sec)\n", timebin, TimeBinsGravity.GlobalNActiveParticles,
             timediff(tstart, tend));
}

/*! \brief Adds individual gravity contribution and appropriate factors.
 *
 *  Routine combines accelerations of particle mesh and tree and applies
 *  the required physical constants and scaling factors e.g. for a cosmological
 *  simulation with nonperiodic gravity.
 *
 *  \param[in] timebin Current timebin for which gravity is calculated.
 *
 *  \return void
 */
void gravity_force_finalize(int timebin)
{
  int i, j, idx;
  double ax, ay, az;

  TIMER_START(CPU_TREE);

  /* now add things for comoving integration */
#ifdef GRAVITY_NOT_PERIODIC
#ifndef PMGRID
  if(All.ComovingIntegrationOn)
    {
      double fac = 0.5 * All.Hubble * All.Hubble * All.Omega0 / All.G;

      for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
        {
          i = TimeBinsGravity.ActiveParticleList[idx];
          if(i < 0)
            continue;

          for(j = 0; j < 3; j++)
            P[i].GravAccel[j] += fac * P[i].Pos[j];
        }
    }
#endif /* #ifndef PMGRID */
#endif /* #ifdef GRAVITY_NOT_PERIODIC */

#ifdef HIERARCHICAL_GRAVITY
  if(timebin == All.HighestOccupiedGravTimeBin)
#endif /* #ifdef HIERARCHICAL_GRAVITY */
    {
      mpi_printf("GRAVTREE: Setting OldAcc!\n");

      for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
        {
          i = TimeBinsGravity.ActiveParticleList[idx];
          if(i < 0)
            continue;

#ifdef PMGRID
          ax = P[i].GravAccel[0] + P[i].GravPM[0] / All.G;
          ay = P[i].GravAccel[1] + P[i].GravPM[1] / All.G;
          az = P[i].GravAccel[2] + P[i].GravPM[2] / All.G;
#else  /* #ifdef PMGRID */
        ax = P[i].GravAccel[0];
        ay = P[i].GravAccel[1];
        az = P[i].GravAccel[2];
#endif /* #ifdef PMGRID #else */

          P[i].OldAcc = sqrt(ax * ax + ay * ay + az * az);
        }
    }

  /*  muliply by G */
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      for(j = 0; j < 3; j++)
        P[i].GravAccel[j] *= All.G;

#ifdef EVALPOTENTIAL

#if defined(PMGRID) && !defined(GRAVITY_NOT_PERIODIC)
      P[i].Potential += All.MassPMregions[0] * M_PI / (All.Asmth[0] * All.Asmth[0] * boxSize_X * boxSize_Y * boxSize_Z);
#ifdef PLACEHIGHRESREGION
      P[i].Potential += All.MassPMregions[1] * M_PI / (All.Asmth[1] * All.Asmth[1] * boxSize_X * boxSize_Y * boxSize_Z);
#endif /* #ifdef PLACEHIGHRESREGION */
#endif /* #if defined(PMGRID) && !defined(GRAVITY_NOT_PERIODIC) */

      /* It's better to not remove the self-potential here to get a smooth potential field for co-spatial particles with varying mass
       * or softening. For calculating the binding energy of a particle, the self-energy should then be removed as
       *
       *  P[i].Potential += P[i].Mass / (All.ForceSoftening[P[i].SofteningType] / 2.8);
       */

      P[i].Potential *= All.G;

#ifdef PMGRID
#ifndef FORCETEST_TESTFORCELAW
      P[i].Potential += P[i].PM_Potential; /* add in long-range potential */
#endif                                     /* #ifndef FORCETEST_TESTFORCELAW */
#endif                                     /* #ifdef PMGRID */
#endif                                     /* #ifdef EVALPOTENTIAL */
      if(All.ComovingIntegrationOn)
        {
#ifdef GRAVITY_NOT_PERIODIC
          double fac, r2;
          int k;

          fac = -0.5 * All.Omega0 * All.Hubble * All.Hubble;

          for(k = 0, r2 = 0; k < 3; k++)
            r2 += P[i].Pos[k] * P[i].Pos[k];

#ifdef EVALPOTENTIAL
          P[i].Potential += fac * r2;
#endif /* #ifdef EVALPOTENTIAL */
#endif /* #ifdef GRAVITY_NOT_PERIODIC */
        }
      else
        {
          double fac, r2;
          int k;

          fac = -0.5 * All.OmegaLambda * All.Hubble * All.Hubble;

          if(fac != 0)
            {
              for(k = 0, r2 = 0; k < 3; k++)
                r2 += P[i].Pos[k] * P[i].Pos[k];
#ifdef EVALPOTENTIAL
              P[i].Potential += fac * r2;
#endif /* #ifdef EVALPOTENTIAL */
            }
        }
    }

    /* Finally, the following factor allows a computation of a cosmological
     * simulation with vacuum energy in physical coordinates
     */
#ifdef GRAVITY_NOT_PERIODIC
#ifndef PMGRID
  if(All.ComovingIntegrationOn == 0)
    {
      double fac = All.OmegaLambda * All.Hubble * All.Hubble;

      for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
        {
          i = TimeBinsGravity.ActiveParticleList[idx];
          if(i < 0)
            continue;

          for(j = 0; j < 3; j++)
            P[i].GravAccel[j] += fac * P[i].Pos[j];
        }
    }
#endif /* #ifndef PMGRID */
#endif /* #ifdef GRAVITY_NOT_PERIODIC */

#ifdef FDM
  fdm_add_force_to_active_stars();
#endif /* #ifdef FDM */

  TIMER_STOP(CPU_TREE);
}
