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

/*! \brief Spawn a black hole particle from the densest gas cell in a FOF group.
 *
 *  Finds the densest gas cell belonging to this group across all MPI tasks,
 *  then the owning task spawns a Type 5 BH particle from it, conserving mass.
 *
 *  \param[in] grp_index  Index into Group[] for this group.
 *  \return void
 */

// #ifdef BLACKHOLE_SEEDING

static int bhs_spawned;             /*!< local number of black hole particles spawned in the time step */
static int tot_bhs_spawned;         /*!< global number of black hole particles spawned in the time step */
static int bhs_converted;           /*!< local number of gas cells converted into black holes in the time step */
static int tot_bhs_converted;       /*!< global number of gas cells converted into black holes in the time step */
static int altogether_spawned;      /*!< local number of particles spawned in the time step */
static int tot_altogether_spawned;  /*!< global number of particles spawned in the time step */
static double cum_mass_bhs = 0.0;   /*!< cumulative mass of black holes created in the time step (global value) */

void seed_black_hole_in_group(int grp_index, int *n_seeded)
{
  /* grp_index == -1 means this task has no group to seed but must
   * still participate in the MPI_Allreduce collective */
  
  MyIDType target_minid = (grp_index >= 0) ? Group[grp_index].MinID : 0;

  double local_max_density = -1.0;
  int    local_best_index  = -1;

  if(grp_index >= 0)
    {
      for(int i = 0; i < NumPart; i++)
        {
          if(P[i].Type != 0) continue;
          if(P[i].Mass == 0 && P[i].ID == 0) continue;
          if(FOF_PList[i].MinID != target_minid) continue;
          if(SphP[i].Density > local_max_density)
            {
              local_max_density = SphP[i].Density;
              local_best_index  = i;
            }
        }
    }

  struct { double density; int task; } local_val, global_val;
  local_val.density = (local_best_index >= 0) ? local_max_density : -1.0;
  local_val.task    = ThisTask;

  MPI_Allreduce(&local_val, &global_val, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

  /* Rest of spawning logic unchanged - only winning task acts */
 if(global_val.density < 0)
    {
      mpi_printf("FOF_SEEDING: WARNING - no gas cell found for group MinID=%llu, skipping.\n",
                 (unsigned long long)target_minid);
      return;  /* all tasks return — no winner */
    }

  /* Non-winning tasks and dummy tasks exit here — winning task continues */
  if(grp_index < 0 || ThisTask != global_val.task)
    return;

  /* --- Phase 3: winning task spawns the BH --- */
  if(ThisTask != global_val.task)
    return;  /* nothing to do on non-winning tasks */

  /* Check we have room — same pattern as make_star() */
  if(NumPart >= All.MaxPart)
    terminate("FOF_SEEDING: no space to spawn BH (NumPart=%d MaxPart=%d)", NumPart, All.MaxPart);

#ifdef BLACKHOLES
  /* Check BhP array has room */
  if(NumBhs >= All.MaxPartBhs)
    {
      All.MaxPartBhs = (int)(1.25 * All.MaxPartBhs) + 1;
      reallocate_memory_maxpartbhs();
    }
#endif

  int igas = local_best_index;

/* store original mass */
double gas_mass = P[igas].Mass;

/* convert gas → BH */
P[igas].Type = 5;
P[igas].SofteningType = All.SofteningTypeOfPartType[5];

/* assign BH mass */
double actual_seed_mass = All.BlackHoleSeedMass;
if(actual_seed_mass >= gas_mass)
  actual_seed_mass = 0.1 * gas_mass;

P[igas].Mass = actual_seed_mass;

/* DO NOT touch SphP extensively — just neutralize */
SphP[igas].Density = 0;
SphP[igas].Volume  = 0;

/* initialise BH structure */
// P[igas].BhID = NumBhs;
// BhP[NumBhs].PID = igas;
// BhP[NumBhs].Hsml = cbrt((3.0 * SphP[igas].Volume) / (4.0 * M_PI));
// BhP[NumBhs].Density = 0;

// NumBhs++;
/* Safety check — mass must not go to zero */
  if(P[igas].Mass <= 0)
    terminate("FOF_SEEDING: gas cell mass went to zero after BH seeding");
  
  mpi_printf("FOF_SEEDING: donor gas cell mass after seeding: %g\n", P[igas].Mass);

  mpi_printf("FOF_SEEDING: Seeded BH (ID=%llu) in group MinID=%llu, "
             "mass=%g at task %d\n",
             (unsigned long long)P[igas].ID,
             (unsigned long long)target_minid,
             All.BlackHoleSeedMass, ThisTask);


  (*n_seeded)++;
}
// #endif /* #ifdef BLACKHOLE_SEEDING */
