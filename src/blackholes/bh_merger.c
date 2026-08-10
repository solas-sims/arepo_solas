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
 *              Hsml/softening (mirroring the global event-list idiom in fof_seeding.c/
 *              bh_seed.c), and then every task independently runs an identical O(N_bh^2)
 *              scan over that snapshot to decide which pairs merge. Because the snapshot
 *              and the scan are both identical on every task, no further communication is
 *              needed to agree on the outcome -- each task only has to apply the decisions
 *              that involve one of its own local particles.
 *
 *              The "close enough to merge" radius is configurable (see bh_merger_radius()
 *              and All.BhMergerRadiusCriterion / enum bh_merger_radius_criterion in bh.h):
 *              by default it is factor*max(Hsml_i,Hsml_j), same as originally, but Hsml is
 *              set by an unrelated criterion (target gas neighbour count) and can diverge
 *              arbitrarily from the gravitational softening length that actually governs
 *              the timestep criterion for a close BH pair, which can let a bound pair sit
 *              at a separation small enough to stall the timestep indefinitely while never
 *              satisfying an Hsml-only merger radius test.
 *
 *              Pair selection is a greedy global match by ascending separation (see Step B/
 *              C below), not a mutual-nearest-neighbour rule: mutual-NN can permanently
 *              deadlock compact 3+ body subsystems (if A's closest bound neighbour is B, but
 *              B's closest bound neighbour is a closer C, A-B never merges even though it
 *              independently qualifies on its own), which is a real, observed failure mode,
 *              not a hypothetical one.
 */

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

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
  MyDouble Softening;
} BhMergerCandidate;

/*! \brief One qualifying (bound, within-radius) candidate pair, indices into global_bhs[]. */
typedef struct
{
  int i, j;
  double r2;
} BhMergerPair;

/*! \brief Gravitational force-softening length for the BH at local BhP[]/PPB() index idx.
 *  See declaration in bh_proto.h for why this isn't static.
 */
double bh_softening_for_index(int idx)
{
  return All.ForceSoftening[PPB(idx).SofteningType];
}

/*! \brief Computes the "close enough to merge" radius for a candidate pair, per
 *         All.BhMergerRadiusCriterion (enum bh_merger_radius_criterion in bh.h).
 *
 *  BH_MERGER_RADIUS_HSML (the default) reproduces the original factor*max(Hsml_i,Hsml_j)
 *  behaviour exactly, so a parameter file that sets BhMergerRadiusCriterion=HSML gets
 *  identical merger behaviour to before this option existed.
 *
 *  \param[in] a First candidate.
 *  \param[in] b Second candidate.
 *
 *  \return The merger radius for this pair (All.BhMergerRadiusFactor already applied).
 */
static double bh_merger_radius(const BhMergerCandidate *a, const BhMergerCandidate *b)
{
  double h_hsml = fmax(a->Hsml, b->Hsml);
  double h_soft = fmax(a->Softening, b->Softening);

  double h;
  switch(All.BhMergerRadiusCriterion)
    {
      case BH_MERGER_RADIUS_HSML:
        h = h_hsml;
        break;
      case BH_MERGER_RADIUS_SOFTENING:
        h = h_soft;
        break;
      case BH_MERGER_RADIUS_MAX_HSML_SOFTENING:
        h = fmax(h_hsml, h_soft);
        break;
      case BH_MERGER_RADIUS_MIN_HSML_SOFTENING:
        h = fmin(h_hsml, h_soft);
        break;
      case BH_MERGER_RADIUS_OTHER_PHYSICAL:
      default:
        /* stub for a future physically-motivated criterion (e.g. mutual Bondi radius);
         * check_parameters() already rejects unrecognised strings, so reaching this means
         * BH_MERGER_RADIUS_OTHER_PHYSICAL was explicitly selected -- not implemented yet. */
        terminate(
            "BH_MERGER: BhMergerRadiusCriterion selects an unimplemented criterion "
            "(BH_MERGER_RADIUS_OTHER_PHYSICAL). Use HSML, SOFTENING, MAX_HSML_SOFTENING, or "
            "MIN_HSML_SOFTENING.\n");
        h = 0; /* unreachable; terminate() aborts */
    }

  return All.BhMergerRadiusFactor * h;
}

/*! \brief qsort() comparator for BhMergerPair, sorting by ascending r2.
 *
 *  global_bhs[] is byte-identical and identically ordered on every task, so every task's
 *  candidate list is identical too -- but qsort() is not guaranteed stable, and different
 *  qsort implementations (or even the same implementation on different inputs) can order
 *  exact ties differently. Break ties deterministically by (i,j) index so every task's
 *  sorted list is byte-identical regardless of qsort implementation details.
 */
static int bh_merger_pair_compare(const void *va, const void *vb)
{
  const BhMergerPair *a = (const BhMergerPair *)va;
  const BhMergerPair *b = (const BhMergerPair *)vb;

  if(a->r2 < b->r2)
    return -1;
  if(a->r2 > b->r2)
    return +1;

  if(a->i != b->i)
    return (a->i < b->i) ? -1 : +1;

  if(a->j != b->j)
    return (a->j < b->j) ? -1 : +1;

  return 0;
}

/*! \brief Merges gravitationally-bound black hole pairs closer than the radius
 *         returned by bh_merger_radius() for that pair (see All.BhMergerRadiusCriterion).
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

      local_bhs[n_local].ID        = PPB(i).ID;
      local_bhs[n_local].Task      = ThisTask;
      local_bhs[n_local].Index     = i;
      local_bhs[n_local].Mass      = PPB(i).Mass;
      local_bhs[n_local].Hsml      = BhP[i].Hsml;
      local_bhs[n_local].Softening = bh_softening_for_index(i);
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

  /* mymalloc()/myfree() is a strict LIFO stack: counts/bcounts/bdispls/local_bhs can't be
   * freed yet even though Allgatherv is their last use, because global_bhs (allocated after
   * them) is still on top of the stack and still needed through Step B/C below. They're
   * freed at the very end, in exact reverse allocation order alongside everything else. */

  /* Step B: every task independently computes the identical qualifying-pair list.
   * global_bhs[] is byte-identical on every task (same Allgatherv, same task-then-local
   * ordering), so this scan produces the same result everywhere without further MPI. */

  int max_pairs = n_global * (n_global - 1) / 2;
  BhMergerPair *pairs = (BhMergerPair *)mymalloc("bh_merger_pairs", (max_pairs > 0 ? max_pairs : 1) * sizeof(BhMergerPair));
  int n_pairs = 0;

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

          double hmax = bh_merger_radius(&global_bhs[i], &global_bhs[j]);

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

          pairs[n_pairs].i  = i;
          pairs[n_pairs].j  = j;
          pairs[n_pairs].r2 = r2;
          n_pairs++;
        }
    }

  qsort(pairs, n_pairs, sizeof(BhMergerPair), bh_merger_pair_compare);

  /* Step C: greedy global match by ascending separation, then deterministic apply.
   *
   * A mutual-nearest-neighbour rule (only merge a pair if each side's single closest bound
   * candidate is the other) can permanently deadlock compact 3+ body subsystems: e.g. if A's
   * closest bound neighbour is B, but B's closest bound neighbour is a closer C, A-B never
   * merges even though the A-B pair independently satisfies the bound + within-radius test
   * on its own. This is a real, observed failure mode (runs stalling indefinitely at the
   * minimum timebin with exactly 3-4 active BH particles), not a hypothetical one.
   *
   * Greedy-by-separation instead walks the sorted qualifying-pair list once, accepting the
   * globally closest pair first, then the next-closest pair whose both members are still
   * unconsumed, and so on. Every accepted pair still independently satisfies the bound +
   * within-radius test, and no BH merges more than once per call, but a BH is no longer
   * required to be any particular *other* BH's single closest neighbour to merge. Do not
   * "simplify" this back to mutual-NN bookkeeping -- that's the deadlock this replaced. */

  char *consumed = (char *)mymalloc("bh_merger_consumed", (n_global > 0 ? n_global : 1) * sizeof(char));
  memset(consumed, 0, n_global * sizeof(char));

  for(int p = 0; p < n_pairs; p++)
    {
      i = pairs[p].i;
      int j = pairs[p].j;

      if(consumed[i] || consumed[j]) /* one or both already merged this call */
        continue;

      consumed[i] = consumed[j] = 1;

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

  myfree(consumed);
  myfree(pairs);
  myfree(global_bhs);
  myfree(bdispls);
  myfree(bcounts);
  myfree(counts);
  myfree(local_bhs);
}

#endif /* #if defined(BH_MERGER) */
