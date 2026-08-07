/*!
 * \file        src/blackholes/bh_seed.c
 * \brief       Seed black holes in FOF halos identified on-the-fly.
 * \details     contains functions:
 *                void seed_black_holes_from_events(HaloSeedEvent *events, int n_events)
 *                static void spawn_black_hole_from_cell(int igas, int ibh, double seed_mass)
 *
 *              Seeding mirrors the star-formation spawn path (make_star() /
 *              spawn_star_from_cell()): a new Type-5 particle is spawned from
 *              the densest gas cell of the halo, removing the seed mass (and
 *              conserved quantities in proportion) from the donor cell rather
 *              than converting it, so that mass is conserved and no mesh cell
 *              is destroyed.
 */

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"

#include "../fof/fof.h"
#include "../fof/fof_seeding.h"

#if defined(HALO_SEEDING) && defined(BLACKHOLE_SEEDING)

static int bhs_spawned;            /*!< local number of black holes spawned at this seeding event */
static int tot_bhs_spawned;        /*!< global number of black holes spawned at this seeding event */
static double cum_mass_bhs = 0.0;  /*!< cumulative seeded black hole mass (global, task 0 log only) */

static void spawn_black_hole_from_cell(int igas, int ibh, double seed_mass);

/*! \brief Spawn black holes for the seed events of this FOF pass.
 *
 *  Collective call: all tasks must enter with the same global event list, as
 *  produced by fof_seeding_list(). Each task spawns the black holes whose
 *  donor gas cells it owns; new particle IDs are then assigned collectively
 *  (same scheme as star formation).
 *
 *  Must be called at the same synchronisation point as the FOF pass that
 *  produced the events (donor indices become stale after the next domain
 *  decomposition).
 *
 *  \param[in] events   Global list of seed events (identical on all tasks).
 *  \param[in] n_events Number of events (identical on all tasks).
 *
 *  \return void
 */
void seed_black_holes_from_events(HaloSeedEvent *events, int n_events)
{
  int i;

  if(n_events == 0)
    return; /* consistent on all tasks */

  bhs_spawned = 0;
  double sum_mass_seeded = 0;

  for(i = 0; i < n_events; i++)
    {
      if(events[i].DonorTask != ThisTask)
        continue;

      int igas = events[i].DonorIndex;

      /* validate the donor cell; indices come from this timestep's FOF pass */
      if(igas < 0 || igas >= NumGas || P[igas].Type != 0 || P[igas].ID != events[i].DonorID ||
         (P[igas].Mass == 0 && P[igas].ID == 0))
        {
          printf("FOF_SEEDING: WARNING Task %d: stale donor cell (index=%d) for halo MinID=%llu, skipping spawn!\n", ThisTask, igas,
                 (unsigned long long)events[i].HaloMinID);
          continue;
        }

      if(NumPart + bhs_spawned >= All.MaxPart)
        terminate("FOF_SEEDING: NumPart=%d spawn %d particles no space left (All.MaxPart=%d)", NumPart, bhs_spawned, All.MaxPart);

      double seed_mass = All.BlackHoleSeedMass;

      if(seed_mass >= 0.9 * P[igas].Mass)
        {
          printf("FOF_SEEDING: WARNING Task %d: BlackHoleSeedMass=%g >= 0.9 x donor cell mass=%g, capping seed mass.\n", ThisTask,
                 All.BlackHoleSeedMass, P[igas].Mass);
          seed_mass = 0.9 * P[igas].Mass;
        }

      int ibh = NumPart + bhs_spawned; /* index of the new black hole */

      spawn_black_hole_from_cell(igas, ibh, seed_mass);

      printf("FOF_SEEDING: Task %d: seeded BH (mass=%g) from cell ID=%llu in halo MinID=%llu (M=%g)\n", ThisTask, seed_mass,
             (unsigned long long)events[i].DonorID, (unsigned long long)events[i].HaloMinID, events[i].HaloMass);

      sum_mass_seeded += seed_mass;
      bhs_spawned++;
    }

  MPI_Allreduce(&bhs_spawned, &tot_bhs_spawned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(tot_bhs_spawned > 0)
    {
      /* assign new unique IDs to the spawned black holes (same scheme as SFR) */
      int *list;

      if(All.MaxID == 0) /* MaxID not calculated yet */
        calculate_maxid();

      list = mymalloc("list", NTask * sizeof(int));

      MPI_Allgather(&bhs_spawned, 1, MPI_INT, list, 1, MPI_INT, MPI_COMM_WORLD);

      MyIDType newid = All.MaxID + 1;

      for(i = 0; i < ThisTask; i++)
        newid += list[i];

      myfree(list);

      for(i = 0; i < bhs_spawned; i++)
        {
          P[NumPart + i].ID = newid;
          newid++;
        }

      All.MaxID += tot_bhs_spawned;

      All.TotNumPart += tot_bhs_spawned;
      NumPart += bhs_spawned;
    }

  double tot_mass_seeded;
  MPI_Reduce(&sum_mass_seeded, &tot_mass_seeded, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    cum_mass_bhs += tot_mass_seeded;

  mpi_printf("FOF_SEEDING: spawned %d black holes out of %d seed events (cumulative seeded mass %g)\n", tot_bhs_spawned, n_events,
             cum_mass_bhs);
}

/*! \brief Spawn a black hole particle from a gas cell.
 *
 *  Mirror of spawn_star_from_cell(): the new Type-5 particle takes the seed
 *  mass from the donor cell, whose conserved quantities are reduced in
 *  proportion. The donor cell remains in the mesh.
 *
 *  \param[in] igas      Index of the donor gas cell.
 *  \param[in] ibh       Index for the new black hole particle.
 *  \param[in] seed_mass Mass of the new black hole (must be < cell mass).
 *
 *  \return void
 */
static void spawn_black_hole_from_cell(int igas, int ibh, double seed_mass)
{
  P[ibh]               = P[igas];
  P[ibh].Type          = 5;
  P[ibh].SofteningType = All.SofteningTypeOfPartType[P[ibh].Type];
  P[ibh].Mass          = seed_mass;

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
  if(((1 << P[ibh].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
    P[ibh].SofteningType = get_softening_type_from_mass(P[ibh].Mass);
#endif /* #ifdef INDIVIDUAL_GRAVITY_SOFTENING */

  timebin_add_particle(&TimeBinsGravity, ibh, igas, P[ibh].TimeBinGrav, TimeBinSynchronized[P[ibh].TimeBinGrav]);

  /* now change the conserved quantities in the cell in proportion */
  double fac = (P[igas].Mass - seed_mass) / P[igas].Mass;

  P[igas].Mass *= fac;
  SphP[igas].Energy *= fac;
  SphP[igas].Momentum[0] *= fac;
  SphP[igas].Momentum[1] *= fac;
  SphP[igas].Momentum[2] *= fac;

#ifdef METALS
  SphP[igas].Metals *= fac;
#endif /* #ifdef METALS */

#ifdef MAXSCALARS
  for(int s = 0; s < N_Scalar; s++) /* Note, the changes in MATERIALS, HIGHRESGASMASS, etc., are treated as part of the Scalars */
    *(MyFloat *)(((char *)(&SphP[igas])) + scalar_elements[s].offset_mass) *= fac;
#endif /* #ifdef MAXSCALARS */

#ifdef BLACKHOLES
  /* register the new black hole in the BhP array */
  if(NumBhs >= All.MaxPartBhs)
    {
      All.MaxPartBhs = (int)(1.25 * All.MaxPartBhs) + 1;
      reallocate_memory_maxpartbhs();
    }

  P[ibh].BhID     = NumBhs;
  BhP[NumBhs].PID = ibh;

  /* initial guess for the smoothing length: radius of the donor cell */
  BhP[NumBhs].Hsml        = cbrt((3.0 * SphP[igas].Volume) / (4.0 * M_PI));
  BhP[NumBhs].Density     = 0;
  BhP[NumBhs].NgbMass     = 0;
  BhP[NumBhs].NgbMassFeed = 0;
  BhP[NumBhs].NgbMinStep  = 0;
  BhP[NumBhs].DensityFlag = 0;
  BhP[NumBhs].TimeBinBh   = 0;

  /* timebin_add_particle(&TimeBinsBh, NumBhs, -1, 0, 1); */

  NumBhs++;
#endif /* #ifdef BLACKHOLES */

  return;
}

#endif /* #if defined(HALO_SEEDING) && defined(BLACKHOLE_SEEDING) */
