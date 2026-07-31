/*!
 * \file        src/blackholes/bh_merger.c
 * \brief       Merge gravitationally-bound, close black hole pairs.
 * \details     contains functions:
 *                void bh_merger(void)
 *
 *              There is no on-the-fly BH-BH interaction search elsewhere in the code:
 *              the accretion/swallow neighbour search (bh_accretion.c/bh_swallow.c) only
 *              ever considers Type-0 gas neighbours (ngb_treefind_variable_threads() is
 *              gas-only by construction). FOF_SEEDING already refuses to seed a halo that
 *              already hosts a black hole (see fof_seeding.c), but nothing merges two black
 *              holes that end up in the same halo after their respective progenitor halos
 *              merge -- they simply coexist indefinitely, each with its own accretion
 *              neighbourhood, which can lead to several black holes competing to drain the
 *              same shared gas cells (see bh_swallow.c/bh_update.c's BhMassDrain handling).
 *
 *              Since the global black hole count is small compared to the gas particle
 *              count, this does not attempt a tree-based neighbour search: every task
 *              MPI_Allgathers a snapshot of every live black hole's position/velocity/mass/
 *              Hsml (mirroring the global event-list idiom in fof_seeding.c/bh_seed.c), and
 *              then every task independently runs an identical O(N_bh^2) scan over that
 *              snapshot to decide which pairs merge. Because the snapshot and the scan are
 *              both identical on every task, no further communication is needed to agree on
 *              the outcome -- each task only has to apply the decisions that involve one of
 *              its own local particles.
 */

#include <math.h>
#include <mpi.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#if defined(BH_MERGER)

/*! \brief Snapshot of one live black hole, as broadcast to every task. */
typedef struct
{
  MyIDType ID;
  int Task;
  int Index; /*!< this BH's index into the owning task's local BhP[] array */
  MyDouble Mass;
  MyDouble Pos[3];
  MyDouble Vel[3];
  MyDouble Hsml;
} BhMergerCandidate;

/*! \brief Merges gravitationally-bound black hole pairs closer than
 *         All.BhMergerRadiusFactor x max(Hsml_i, Hsml_j).
 *
 *  Collective call: every task must enter this function together (it does two
 *  collective MPI calls). Should be called once per BH-active step, after
 *  bh_perform_end_of_step_physics() has already applied this step's accretion
 *  onto every currently active black hole.
 *
 *  \return void
 */
void bh_merger(void)
{
  int i, k;

  /* Step A: collect this task's live black holes and share them with every task */

  BhMergerCandidate *local_bhs =
      (BhMergerCandidate *)mymalloc("bh_merger_local", (NumBhs > 0 ? NumBhs : 1) * sizeof(BhMergerCandidate));

  int n_local = 0;

  for(i = 0; i < NumBhs; i++)
    {
      if(PPB(i).Mass == 0 && PPB(i).ID == 0) /* already killed (e.g. by an earlier merge this run) */
        continue;

      local_bhs[n_local].ID    = PPB(i).ID;
      local_bhs[n_local].Task  = ThisTask;
      local_bhs[n_local].Index = i;
      local_bhs[n_local].Mass  = PPB(i).Mass;
      local_bhs[n_local].Hsml  = BhP[i].Hsml;
      for(k = 0; k < 3; k++)
        {
          local_bhs[n_local].Pos[k] = PPB(i).Pos[k];
          local_bhs[n_local].Vel[k] = PPB(i).Vel[k];
        }

      n_local++;
    }

  int *counts = (int *)mymalloc("bh_merger_counts", NTask * sizeof(int));

  MPI_Allgather(&n_local, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);

  int *bcounts = (int *)mymalloc("bh_merger_bcounts", NTask * sizeof(int));
  int *bdispls = (int *)mymalloc("bh_merger_bdispls", NTask * sizeof(int));

  int n_global = 0;
  for(i = 0; i < NTask; i++)
    {
      bcounts[i] = counts[i] * sizeof(BhMergerCandidate);
      bdispls[i] = n_global * sizeof(BhMergerCandidate);
      n_global += counts[i];
    }

  if(n_global < 2) /* nothing that could possibly pair up */
    {
      myfree(bdispls);
      myfree(bcounts);
      myfree(counts);
      myfree(local_bhs);
      return;
    }

  BhMergerCandidate *global_bhs = (BhMergerCandidate *)mymalloc("bh_merger_global", n_global * sizeof(BhMergerCandidate));

  MPI_Allgatherv(local_bhs, n_local * sizeof(BhMergerCandidate), MPI_BYTE, global_bhs, bcounts, bdispls, MPI_BYTE,
                 MPI_COMM_WORLD);

  myfree(bdispls);
  myfree(bcounts);
  myfree(counts);
  myfree(local_bhs);

  /* Step B: every task independently computes the identical accepted-pair list.
   * global_bhs[] is byte-identical on every task (same Allgatherv, same task-then-local
   * ordering), so this scan produces the same result everywhere without further MPI. */

  int *best_partner = (int *)mymalloc("bh_merger_best_partner", n_global * sizeof(int));
  double *best_r2    = (double *)mymalloc("bh_merger_best_r2", n_global * sizeof(double));

  for(i = 0; i < n_global; i++)
    {
      best_partner[i] = -1;
      best_r2[i]      = 0;
    }

  MyDouble xtmp, ytmp, ztmp;

  for(i = 0; i < n_global; i++)
    {
      for(int j = i + 1; j < n_global; j++)
        {
          double dx = NEAREST_X(global_bhs[i].Pos[0] - global_bhs[j].Pos[0]);
          double dy = NEAREST_Y(global_bhs[i].Pos[1] - global_bhs[j].Pos[1]);
          double dz = NEAREST_Z(global_bhs[i].Pos[2] - global_bhs[j].Pos[2]);
          double r2 = dx * dx + dy * dy + dz * dz;

          if(r2 <= 0) /* degenerate/identical positions: nothing sane to compute */
            continue;

          double hmax = All.BhMergerRadiusFactor * fmax(global_bhs[i].Hsml, global_bhs[j].Hsml);

          if(r2 >= hmax * hmax)
            continue;

          double r = sqrt(r2);

          double dvx = global_bhs[i].Vel[0] - global_bhs[j].Vel[0];
          double dvy = global_bhs[i].Vel[1] - global_bhs[j].Vel[1];
          double dvz = global_bhs[i].Vel[2] - global_bhs[j].Vel[2];
          double v2  = dvx * dvx + dvy * dvy + dvz * dvz;

          /* gravitationally-bound (two-body escape velocity) check */
          double vesc2 = 2.0 * All.G * (global_bhs[i].Mass + global_bhs[j].Mass) / r;

          if(v2 >= vesc2)
            continue;

          if(best_partner[i] < 0 || r2 < best_r2[i])
            {
              best_partner[i] = j;
              best_r2[i]      = r2;
            }
          if(best_partner[j] < 0 || r2 < best_r2[j])
            {
              best_partner[j] = i;
              best_r2[j]      = r2;
            }
        }
    }

  /* Step C: deterministic apply. Accept a pair only where both sides picked each other
   * (mutual nearest neighbour), and process each accepted pair exactly once. */

  for(i = 0; i < n_global; i++)
    {
      int j = best_partner[i];

      if(j < 0 || j < i) /* no candidate, or this pair already handled from j's side */
        continue;

      if(best_partner[j] != i) /* not mutual */
        continue;

      int survivor, loser;
      if(global_bhs[i].Mass > global_bhs[j].Mass || (global_bhs[i].Mass == global_bhs[j].Mass && global_bhs[i].ID < global_bhs[j].ID))
        {
          survivor = i;
          loser    = j;
        }
      else
        {
          survivor = j;
          loser    = i;
        }

      double mass_new = global_bhs[survivor].Mass + global_bhs[loser].Mass;
      MyDouble vel_new[3];
      for(k = 0; k < 3; k++)
        vel_new[k] = (global_bhs[survivor].Mass * global_bhs[survivor].Vel[k] +
                      global_bhs[loser].Mass * global_bhs[loser].Vel[k]) /
                     mass_new;

      if(global_bhs[survivor].Task == ThisTask)
        {
          int idx = global_bhs[survivor].Index;

          if(PPB(idx).ID == global_bhs[survivor].ID) /* re-validate: still the same particle we saw in Step A */
            {
              PPB(idx).Mass = mass_new;
              for(k = 0; k < 3; k++)
                PPB(idx).Vel[k] = vel_new[k];
              BhP[idx].NumMergers++;

              mpi_printf("BH_MERGER: Task %d: BH ID=%llu (mass=%g) absorbed BH ID=%llu (mass=%g) -> new mass=%g\n", ThisTask,
                         (unsigned long long)global_bhs[survivor].ID, global_bhs[survivor].Mass,
                         (unsigned long long)global_bhs[loser].ID, global_bhs[loser].Mass, mass_new);
            }
        }

      if(global_bhs[loser].Task == ThisTask)
        {
          int idx = global_bhs[loser].Index;

          if(PPB(idx).ID == global_bhs[loser].ID) /* re-validate: still the same particle we saw in Step A */
            {
              /* standard "this particle is dead" signal (mirrors voronoi_derefinement.c's cell
               * dissolution): swept up by the next domain decomposition, which already keeps
               * BhP[] in sync with relocated/eliminated particles. bh_reconstruct_timebins()
               * skips Mass==0/ID==0 entries so this BH drops out of TimeBinsBh from the next
               * rebuild onward; TimeBinsGravity relies on the same generic Mass==0/ID==0 skip
               * that every other particle-kill path in this code already depends on. */
              PPB(idx).Mass   = 0;
              PPB(idx).ID     = 0;
              PPB(idx).Vel[0] = 0;
              PPB(idx).Vel[1] = 0;
              PPB(idx).Vel[2] = 0;
            }
        }
    }

  myfree(best_r2);
  myfree(best_partner);
  myfree(global_bhs);
}

#endif /* #if defined(BH_MERGER) */
