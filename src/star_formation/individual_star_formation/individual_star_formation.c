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
 * \file        src/star_formation/starformation.c
 * \date        05/2018
 * \brief       Generic creation routines for star particles.
 * \details     Star formation rates are calculated in sfr_eEOS for the
 *              multiphase model.
 *              contains functions:
 *                void sfr_init()
 *                void sfr_create_star_particles(void)
 *                void convert_cell_into_star(int i, double birthtime)
 *                void spawn_star_from_cell(int igas, double birthtime, int
 *                  istar, MyDouble mass_of_star)
 *                void make_star(int idx, int i, double prob, MyDouble
 *                  mass_of_star, double *sum_mass_stars)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 07.06.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../main/allvars.h"
#include "../../main/proto.h"

static void spawn_heavy(int igas, double birthtime, int istar, MyDouble mass_of_star);
static void spawn_light(int igas, double birthtime, int istar, MyDouble mass_of_star);
static void make_individual_star(int i, MyDouble mass_of_star, double *local_stars_mass);

static int stars_spawned;      /*!< local number of star particles spawned in the time step */
static int tot_stars_spawned;  /*!< global number of star paricles spawned in the time step */


/*! \brief This routine creates star particles according to their
 *         respective rates.
 *
 *  This function loops over all the active gas cells. If in a given cell the
 *  SFR is greater than zero, the probability of forming a star is computed
 *  and the corresponding particle is created stichastically according to the
 *  model in Springel & Hernquist (2003, MNRAS). It also saves information
 *  about the formed stellar mass and the star formation rate in the file
 *  FdSfr.
 *
 *  \return void
 */
void individual_starbystar_formation(void)
{
  TIMER_START(CPU_COOLINGSFR);

  int idx, i;
  double dt, dtff, u;
  MyDouble mass_of_star;
  double p = 0, prob, p_decide;
  double rate, local_stars_mass, global_stars_mass;

  stars_spawned = local_stars_mass = 0;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      
      /* skip cells that have been swallowed or eliminated */
      if(P[i].Mass == 0 && P[i].ID == 0)
        continue; 

      dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
      dt /= All.cf_hubble_a;

      dtff = sqrt(3.0 * M_PI / 32 / All.G / SphP[i].Density);
          
      p = 0;
      prob = 0;
      mass_of_star = 0;
    
      if(SphP[i].Sfr > 0)
        {
          p = All.StarFormationEfficiency * dt / dtff;
          prob = (1 - exp(-p));
        }

      if(prob == 0)
        continue;

      if(prob < 0)
        terminate("SFR: prob < 0");

      if(prob > 1)
        {
          terminate(
          "Individual Star Formation: need to make a heavier star than desired. Task=%d prob=%g P[i].Mass=%g mass_of_star=%g",
          ThisTask, prob, P[i].Mass, mass_of_star);
        }

      /* decide what process to consider */
      p_decide = get_random_number();

      if(p_decide < prob) /* ok, it is decided to consider star formation */
        {
          u = get_random_number_aux();
          mass_of_star = sample_imf(u);
          mass_of_star /= All.cf_UnitMass_in_Msun;
          make_individual_star(i, mass_of_star, &local_stars_mass);
        }
    } /* end of main loop over active gas particles */

  /* Check if we are overflowing the stars array based on stars actually formed this step */
  int local_star_load = NumStars;
  int global_star_load;

  MPI_Allreduce(&local_star_load, &global_star_load, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if(global_star_load != 0 && global_star_load > (1.0 - ALLOC_TOLERANCE) * All.MaxPartStars)
    {
      int old_alloc = All.MaxPartStars;
      All.MaxPartStars = global_star_load / (1.0 - 2 * ALLOC_TOLERANCE);

      if(All.MaxPartStars < ALLOC_STAR_ROOM)
        All.MaxPartStars = ALLOC_STAR_ROOM;

      if(All.MaxPartStars != old_alloc)
        reallocate_memory_maxpartstars();
    }

  MPI_Allreduce(&stars_spawned, &tot_stars_spawned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(tot_stars_spawned > 0)
    mpi_printf("Individual Star Formation: spawned %d stars\n", tot_stars_spawned);

  if(tot_stars_spawned)
    {
      /* need to assign new unique IDs to the spawned stars */

      int *list;

      if(All.MaxID == 0) /* MaxID not calculated yet */
        calculate_maxid();

      list = mymalloc("list", NTask * sizeof(int));

      MPI_Allgather(&stars_spawned, 1, MPI_INT, list, 1, MPI_INT, MPI_COMM_WORLD);

      MyIDType newid = All.MaxID + 1;

      for(i = 0; i < ThisTask; i++)
        newid += list[i];

      myfree(list);

      for(i = 0; i < stars_spawned; i++)
        {
          P[NumPart + i].ID = newid;

          newid++;
        }

      All.MaxID += tot_stars_spawned;
    }

  /* Note: New tree construction can be avoided because of  `force_add_star_to_tree()' */
  if(tot_stars_spawned > 0)
    {
      All.TotNumPart += tot_stars_spawned;
      NumPart += stars_spawned;

      All.TotNumStars += tot_stars_spawned;
    }

  MPI_Reduce(&local_stars_mass, &global_stars_mass, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      if(All.TimeStep > 0)
        rate = global_stars_mass / All.TimeStep;
      else
        rate = 0;

      /* convert to solar masses per yr */
      rate *= All.cf_UnitMass_in_Msun / All.cf_UnitTime_in_yr;

      fprintf(FdSfr, "%14e %14e %14e\n", All.Time, global_stars_mass, rate);
      myflush(FdSfr);
    }

  sf_starbystar();
  sf_massdrain();

  /* Apply drain and finalize heavy stars */
  for(i = 0; i < NumStars; i++)
    {
      /* Heavy star */
      if(PPS(i).Mass == 0 && SP[i].MassOfStar > 0) 
        {
          PPS(i).Mass = SP[i].MassOfStar;
          SP[i].Metallicity /= PPS(i).Mass;
        }
    }
      
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(SphP[i].StarMassDrain > 0)
        {
          if(P[i].Mass - SphP[i].StarMassDrain < 0.1*P[i].Mass)
            {
              terminate("STAR FORMATION DRAIN ERROR!");
            }
          else
            {
              double factor = (P[i].Mass - SphP[i].StarMassDrain) / P[i].Mass;
                  
              P[i].Mass -= SphP[i].StarMassDrain;
                    
              // Update total energy 
              SphP[i].Energy *= factor;
                    
              // Update momentum 
              SphP[i].Momentum[0] *= factor;
              SphP[i].Momentum[1] *= factor;
              SphP[i].Momentum[2] *= factor;
#ifdef MAXSCALARS
              for(int s = 0; s < N_Scalar; s++)
              *(MyFloat *)(((char *)(&SphP[i])) + scalar_elements[s].offset_mass) *= factor;
#endif
            }
          SphP[i].StarMassDrain = 0;
        }
    }

  TIMER_STOP(CPU_COOLINGSFR);
}

/*! \brief Spawn a star particle from a gas cell.
 *
 *  This function spawns a star particle from an active star-forming
 *  cell. The particle information of the gas cell is copied to the
 *  location istar and the fields necessary for the creation of the star
 *  particle are initialized. The conserved variables of the gas cell
 *  are then updated according to the mass ratio between the two components
 *  to ensure conservation.
 *
 *  \param[in] igas Index of the gas cell from which the star is spawned.
 *  \param[in] birthtime Time of birth (in code units) of the stellar particle.
 *  \param[in] istar Index of the spawned stellar particle.
 *  \param[in] mass_of_star The mass of the spawned stellar particle.
 *
 *  \return void
 */
static void spawn_heavy(int igas, double birthtime, int istar, MyDouble mass_of_star)
{  
  P[istar] = P[igas];
  P[istar].Type = 4;
  P[istar].Mass = 0;

  P[istar].SofteningType = All.SofteningTypeOfPartType[P[istar].Type];
  
#ifdef INDIVIDUAL_GRAVITY_SOFTENING
  if(((1 << P[istar].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
    P[istar].SofteningType = get_softening_type_from_mass(P[istar].Mass);
#endif /* #ifdef INDIVIDUAL_GRAVITY_SOFTENING */

  timebin_add_particle(&TimeBinsGravity, istar, igas, P[istar].TimeBinGrav, TimeBinSynchronized[P[istar].TimeBinGrav]);
  
  /* Zero star struct */
  memset(&SP[NumStars], 0, sizeof(Star_Particle_Data));
  /* Assign star_ids */
  P[istar].SID = NumStars;
  SP[NumStars].PID = istar;
  
  /* Prepare for star forming loop */
  SP[NumStars].MassOfStar = mass_of_star;
  SP[NumStars].Hsml = get_cell_radius(igas);
  
#ifdef STAR_FEEDBACK_ACTIVE
  /* Set timebin */
  SP[NumStars].Active = 0;
  SP[NumStars].WithFeedback = 1;
  SP[NumStars].HostHydroBin = P[igas].TimeBinHydro;
  timebin_add_particle(&TimeBinsStar, NumStars, -1, 0, 1);
#endif

  NumStars++;

  return;
}

/*! \brief Spawn a star particle from a gas cell.
 *
 *  This function spawns a star particle from an active star-forming
 *  cell. The particle information of the gas cell is copied to the
 *  location istar and the fields necessary for the creation of the star
 *  particle are initialized. The conserved variables of the gas cell
 *  are then updated according to the mass ratio between the two components
 *  to ensure conservation.
 *
 *  \param[in] igas Index of the gas cell from which the star is spawned.
 *  \param[in] birthtime Time of birth (in code units) of the stellar particle.
 *  \param[in] istar Index of the spawned stellar particle.
 *  \param[in] mass_of_star The mass of the spawned stellar particle.
 *
 *  \return void
 */
static void spawn_light(int igas, double birthtime, int istar, MyDouble mass_of_star)
{
  P[istar] = P[igas];
  P[istar].Type = 4;
  P[istar].Mass = mass_of_star;

  P[istar].SofteningType = All.SofteningTypeOfPartType[P[istar].Type];

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
  if(((1 << P[istar].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
    P[istar].SofteningType = get_softening_type_from_mass(P[istar].Mass);
#endif /* #ifdef INDIVIDUAL_GRAVITY_SOFTENING */

  timebin_add_particle(&TimeBinsGravity, istar, igas, P[istar].TimeBinGrav, TimeBinSynchronized[P[istar].TimeBinGrav]);

  /* now change the conserved quantities in the cell in proportion */
  double fac = (P[igas].Mass - P[istar].Mass) / P[igas].Mass;

//#ifdef MHD
//  double Emag = 0.5 * (SphP[igas].B[0] * SphP[igas].B[0] + SphP[igas].B[1] * SphP[igas].B[1] + SphP[igas].B[2] * SphP[igas].B[2]) *
//                SphP[igas].Volume * All.cf_atime;
//  SphP[igas].Energy -= Emag;
//#endif /* #ifdef MHD */

  P[igas].Mass *= fac;
  SphP[igas].Energy *= fac;
  SphP[igas].Momentum[0] *= fac;
  SphP[igas].Momentum[1] *= fac;
  SphP[igas].Momentum[2] *= fac;

//#ifdef MHD
//  SphP[igas].Energy += Emag;
//#endif /* #ifdef MHD */

#ifdef MAXSCALARS
  for(int s = 0; s < N_Scalar; s++) 
    *(MyFloat *)(((char *)(&SphP[igas])) + scalar_elements[s].offset_mass) *= fac;
#endif 

  /* Zero star struct */
  memset(&SP[NumStars], 0, sizeof(Star_Particle_Data));
  /* Assign star_ids */
  P[istar].SID = NumStars;
  SP[NumStars].PID = istar;

#ifdef METALS 
  SP[NumStars].Metallicity = SphP[igas].GasMetals / P[igas].Mass;
#endif

#ifdef STAR_FEEDBACK_ACTIVE
  /* Assign density loop properties */
  SP[NumStars].Hsml = get_cell_radius(igas);
  
  /* Set timebin */
  SP[NumStars].Active = 0;
  SP[NumStars].WithFeedback = 1;
  SP[NumStars].HostHydroBin = P[igas].TimeBinHydro;
  timebin_add_particle(&TimeBinsStar, NumStars, -1, 0, 1);

  /* This is needed for lower res star by star simulations
     Give star small random displacement */
  if(mass_of_star * All.cf_UnitMass_in_Msun > 2)
    {
      double cell_size = get_cell_radius(igas);

      double rx = (rand()/RAND_MAX - 0.5) * cell_size / 50;
      double ry = (rand()/RAND_MAX - 0.5) * cell_size / 50; 
      double rz = (rand()/RAND_MAX - 0.5) * cell_size / 50;

      P[istar].Pos[0] += rx;
      P[istar].Pos[1] += ry;
      P[istar].Pos[2] += rz;
    }
#endif

  NumStars++;

  return;
}

/*! \brief Make a star particle from a gas cell.
 *
 *  Given a gas cell where star formation is active and the probability
 *  of forming a star, this function selectes either to convert the gas
 *  cell into a star particle or to spawn a star depending on the
 *  target mass for the star.
 *
 *  \param[in] idx Index of the gas cell in the hydro list of active cells.
 *  \param[in] i Index of the gas cell.
 *  \param[in] prob Probability of making a star.
 *  \param[in] mass_of_star Desired mass of the star particle.
 *  \param[in, out] sum_mass_stars Holds the mass of all the stars created at the
 *             current time-step (for the local task)
 *
 *  \return void
 */
static void make_individual_star(int i, MyDouble mass_of_star, double *local_stars_mass)
{
  if(NumPart + stars_spawned >= All.MaxPart)
    terminate("NumPart=%d spwawn %d particles no space left (All.MaxPart=%d)\n", NumPart, stars_spawned, All.MaxPart);

  int j = NumPart + stars_spawned;
  *local_stars_mass += mass_of_star;
  stars_spawned++;

  if(mass_of_star < P[i].Mass)
    spawn_light(i, All.Time, j, mass_of_star);
  else 
    spawn_heavy(i, All.Time, j, mass_of_star);
}