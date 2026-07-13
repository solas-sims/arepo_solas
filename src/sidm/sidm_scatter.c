/*!
 * \file        src/sidm/sidm_scatter.c
 * \brief       v1 elastic SIDM: Monte Carlo pairwise scattering/kick,
 *              following Vogelsberger, Zavala & Loeb (2012)'s algorithm
 *              (as described in e.g. Valdarnini 2023, arXiv:2309.10374,
 *              Eq. 18-20 -- not derived independently).
 *
 *              v1 restrictions, both deliberate and documented:
 *              - LOCAL PARTNERS ONLY: a scattering partner must be a
 *                genuinely local P[] particle. Cross-task partners
 *                (reached only via the pseudo-particle branch of the
 *                tree walk) are skipped. This means particles near a
 *                task boundary have a slightly undercounted candidate
 *                pool -- a real, bounded approximation, not a bug,
 *                deferred to a later phase alongside DMSP[] and the
 *                VelDisp mean-subtraction fix.
 *              - LOWER-ID-INITIATES: for a pair (i,j), only the
 *                particle with the smaller ID evaluates and (if
 *                triggered) applies the scatter. This deliberately
 *                REPLACES Vogelsberger's P_i = sum(P_ij)/2 convention
 *                (the /2 there corrects for each pair being evaluated
 *                twice, once from each side) -- since lower-ID-
 *                initiates structurally prevents double-evaluation in
 *                the first place, there is nothing left for a /2 to
 *                correct, so it is dropped. This is a derived
 *                consequence of the parallel-implementation choice,
 *                not a verbatim copy of the source algorithm.
 *
 *              Called once per particle at EXACTLY its own native
 *              timebin (P[i].TimeBinGrav == timebin), not from the
 *              cumulative hierarchical active-list used by gravity's
 *              own kick loop -- that list includes a particle on every
 *              coarser timebin's pass too, which would re-roll its
 *              scatter dice multiple times per actual step. See the
 *              call site in do_gravity_hydro.c for where this matters.
 */

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"
#include "../time_integration/timestep.h"
#include "sidm.h"
#include "sidm_tree.h"

#ifdef SIDM

/*! \brief One candidate scattering partner found during the local walk. */
typedef struct
{
  int index;     /*!< P[] index of the candidate */
  double r;      /*!< distance from the target particle */
  double p_ij;   /*!< this pair's scattering probability contribution */
} sidm_candidate;

static int sidm_candidate_compare(const void *a, const void *b)
{
  const sidm_candidate *ca = (const sidm_candidate *)a;
  const sidm_candidate *cb = (const sidm_candidate *)b;
  if(ca->r < cb->r)
    return -1;
  if(ca->r > cb->r)
    return 1;
  return 0;
}

/*! \brief Physical elapsed time for a full step at the given hierarchical
 *  timebin, in code time units.
 *
 *  Needed because All.Timebase_interval steps are linear in ln(a), not
 *  in physical time, for a comoving run -- (tend-tstart)*Timebase_interval
 *  is Delta(ln a), NOT a time interval, so it cannot be used directly
 *  in the P_ij formula. get_time_difference_in_Gyr() does the actual
 *  cosmological-expansion-aware conversion; the Gyr->code-units step
 *  afterwards is ordinary dimensional bookkeeping (SEC_PER_YEAR,
 *  All.UnitTime_in_s), not something requiring a literature check.
 */
static double sidm_dt_code_units(int timebin)
{
  integertime ti_step = timebin ? (((integertime)1) << timebin) : 0;
  integertime tend     = All.Ti_begstep[timebin];
  integertime tstart    = tend - ti_step;

  double dt_code;

  if(All.ComovingIntegrationOn)
    {
      double a0 = All.TimeBegin * exp(tstart * All.Timebase_interval);
      double a1 = All.TimeBegin * exp(tend * All.Timebase_interval);

      double dt_gyr = get_time_difference_in_Gyr(a0, a1);
      dt_code       = dt_gyr * 1.0e9 * SEC_PER_YEAR / All.UnitTime_in_s;
    }
  else
    {
      dt_code = (tend - tstart) * All.Timebase_interval;
    }

  return dt_code;
}

/*! \brief Draws a uniformly-distributed random unit vector. */
static void sidm_random_unit_vector(double e[3])
{
  double costheta = 2.0 * get_random_number_aux() - 1.0;
  double sintheta = sqrt(fmax(0.0, 1.0 - costheta * costheta));
  double phi      = 2.0 * M_PI * get_random_number_aux();

  e[0] = sintheta * cos(phi);
  e[1] = sintheta * sin(phi);
  e[2] = costheta;
}

/*! \brief Applies the isotropic elastic two-body kick to particles i, j.
 *
 *  Conserves momentum and kinetic energy exactly by construction
 *  (Vogelsberger et al. 2012, Eq. 20) -- verified below on every real
 *  scatter event as a live diagnostic, rather than via a separate
 *  isolated test.
 *
 *  NOTE: assumes m_i == m_j (matches the source paper's own stated
 *  assumption of equal-mass DM particles, and this codebase's
 *  equal-mass DM box). V = (v_i+v_j)/2 and the symmetric 50/50 velocity
 *  split are both mass-weighting shortcuts that would need generalizing
 *  (mass-weighted CM velocity, asymmetric split by mass fraction) if
 *  unequal-mass particles were ever involved -- not needed for v1.
 */
static long long sidm_scatter_count_local     = 0;
static double     sidm_max_momentum_error     = 0.0;
static double     sidm_max_energy_error       = 0.0;

static void sidm_apply_kick(int i, int j)
{
  double vij[3], V[3], e[3];

  /* TEMPORARY diagnostic: momentum/energy before the kick. */
  double p_before[3], ke_before;
  for(int k = 0; k < 3; k++)
    p_before[k] = P[i].Mass * P[i].Vel[k] + P[j].Mass * P[j].Vel[k];
  ke_before = 0.5 * P[i].Mass *
                  (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]) +
              0.5 * P[j].Mass *
                  (P[j].Vel[0] * P[j].Vel[0] + P[j].Vel[1] * P[j].Vel[1] + P[j].Vel[2] * P[j].Vel[2]);

  for(int k = 0; k < 3; k++)
    {
      vij[k] = P[i].Vel[k] - P[j].Vel[k];
      V[k]   = 0.5 * (P[i].Vel[k] + P[j].Vel[k]);
    }

  double v_rel = sqrt(vij[0] * vij[0] + vij[1] * vij[1] + vij[2] * vij[2]);

  sidm_random_unit_vector(e);

  for(int k = 0; k < 3; k++)
    {
      P[i].Vel[k] = V[k] + 0.5 * v_rel * e[k];
      P[j].Vel[k] = V[k] - 0.5 * v_rel * e[k];
    }

  /* TEMPORARY diagnostic: momentum/energy after the kick, compare. */
  {
    double p_after[3], ke_after;
    for(int k = 0; k < 3; k++)
      p_after[k] = P[i].Mass * P[i].Vel[k] + P[j].Mass * P[j].Vel[k];
    ke_after = 0.5 * P[i].Mass *
                   (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]) +
               0.5 * P[j].Mass *
                   (P[j].Vel[0] * P[j].Vel[0] + P[j].Vel[1] * P[j].Vel[1] + P[j].Vel[2] * P[j].Vel[2]);

    double p_mag_before = sqrt(p_before[0] * p_before[0] + p_before[1] * p_before[1] + p_before[2] * p_before[2]);
    double dp[3]         = {p_after[0] - p_before[0], p_after[1] - p_before[1], p_after[2] - p_before[2]};
    double dp_mag        = sqrt(dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2]);
    double p_rel_err     = (p_mag_before > 0) ? dp_mag / p_mag_before : dp_mag;
    double ke_rel_err     = (ke_before > 0) ? fabs(ke_after - ke_before) / ke_before : fabs(ke_after - ke_before);

    if(p_rel_err > sidm_max_momentum_error)
      sidm_max_momentum_error = p_rel_err;
    if(ke_rel_err > sidm_max_energy_error)
      sidm_max_energy_error = ke_rel_err;

    if(p_rel_err > 1e-10 || ke_rel_err > 1e-10)
      printf("SIDM_SCATTER conservation WARNING: task=%d IDs=(%lld,%lld) p_rel_err=%.3e ke_rel_err=%.3e\n", ThisTask,
             (long long)P[i].ID, (long long)P[j].ID, p_rel_err, ke_rel_err);

    sidm_scatter_count_local++;
  }

  P[i].SidmLastScatterTime = All.Ti_Current;
  P[j].SidmLastScatterTime = All.Ti_Current;
  P[i].SidmScatterFlag     = 1;
  P[j].SidmScatterFlag     = 1;
}

/*! \brief Evaluates and, if triggered, performs a scatter for particle i.
 *
 *  Walks SidmTree_Nodes locally (no MPI export/import -- v1 local-
 *  partners-only restriction), collecting candidates with
 *  P[cand].ID > P[i].ID (lower-ID-initiates convention) within
 *  P[i].SidmHsml, then follows Vogelsberger et al. (2012)'s
 *  accept/reject and partner-selection procedure.
 */
static void sidm_scatter_evaluate(int i, double dt_i)
{
  if(P[i].SidmHsml <= 0 || P[i].SidmDensity <= 0)
    return; /* density never computed for this particle yet */

  double pos_x = P[i].Pos[0];
  double pos_y = P[i].Pos[1];
  double pos_z = P[i].Pos[2];
  double h     = P[i].SidmHsml;
  double h2    = h * h;
  double hinv  = 1.0 / h;
  double hinv3 = hinv * hinv * hinv;

  /* Fixed-size stack scratch buffer -- Hsml is tuned for ~SidmDesNumNgb
   * neighbours (order 10s), so 512 is a generous cap; using a plain
   * local array here rather than mymalloc_movable avoids any question
   * of this small, per-call scratch buffer interacting with Arepo's
   * custom LIFO allocator's movable-block ordering assumptions. */
  enum
  {
    MAX_SCATTER_CANDIDATES = 512
  };
  sidm_candidate candidates[MAX_SCATTER_CANDIDATES];
  int            ncand = 0;

  int no = SidmTree_MaxPart; /* root */
  while(no >= 0)
    {
      if(no < SidmTree_MaxPart) /* single particle */
        {
          int p = no;
          no    = SidmTree_Nextnode[no];

          if(p == i || P[p].Type != 1)
            continue;

          if(P[p].ID <= P[i].ID)
            continue; /* lower-ID-initiates: not our candidate to evaluate */

          if(P[p].Ti_Current != All.Ti_Current)
            drift_particle(p, All.Ti_Current);

          double dx = P[p].Pos[0] - pos_x;
          double dy = P[p].Pos[1] - pos_y;
          double dz = P[p].Pos[2] - pos_z;
          double r2 = dx * dx + dy * dy + dz * dz;

          if(r2 >= h2)
            continue;

          double r = sqrt(r2);
          double u = r * hinv;
          double wk;

          if(u < 0.5)
            wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
          else
            wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

          double vx = P[i].Vel[0] - P[p].Vel[0];
          double vy = P[i].Vel[1] - P[p].Vel[1];
          double vz = P[i].Vel[2] - P[p].Vel[2];
          double v_ij = sqrt(vx * vx + vy * vy + vz * vz);

          double p_ij = P[i].Mass * wk * All.SidmCrossSection * v_ij * dt_i;

          if(ncand < MAX_SCATTER_CANDIDATES)
            {
              candidates[ncand].index = p;
              candidates[ncand].r     = r;
              candidates[ncand].p_ij  = p_ij;
              ncand++;
            }
          /* else: silently drop excess candidates beyond MAX_SCATTER_CANDIDATES --
           * acceptable for v1 given Hsml targets ~SidmDesNumNgb << 512;
           * revisit if this cap is ever actually hit. */
        }
      else if(no < SidmTree_MaxPart + SidmTree_MaxNodes) /* internal node */
        {
          struct SidmNODE *current = &SidmTree_Nodes[no];

          double dist_x = current->u.d.range_max[0] < pos_x   ? pos_x - current->u.d.range_max[0]
                           : current->u.d.range_min[0] > pos_x ? current->u.d.range_min[0] - pos_x
                                                                 : 0.0;
          double dist_y = current->u.d.range_max[1] < pos_y   ? pos_y - current->u.d.range_max[1]
                           : current->u.d.range_min[1] > pos_y ? current->u.d.range_min[1] - pos_y
                                                                 : 0.0;
          double dist_z = current->u.d.range_max[2] < pos_z   ? pos_z - current->u.d.range_max[2]
                           : current->u.d.range_min[2] > pos_z ? current->u.d.range_min[2] - pos_z
                                                                 : 0.0;

          if(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z > h2)
            {
              no = current->u.d.sibling;
              continue;
            }

          no = current->u.d.nextnode;
        }
      else /* pseudo-particle: v1 restriction -- local partners only, skip */
        {
          no = SidmTree_Nextnode[no - SidmTree_MaxNodes];
        }
    }

  if(ncand == 0)
    return;

  double P_i = 0.0;
  for(int c = 0; c < ncand; c++)
    P_i += candidates[c].p_ij;

  double x = get_random_number_aux();

  if(x >= P_i)
    return; /* no scatter this step */

  qsort(candidates, ncand, sizeof(sidm_candidate), sidm_candidate_compare);

  double cumsum = 0.0;
  int    chosen = -1;
  for(int c = 0; c < ncand; c++)
    {
      cumsum += candidates[c].p_ij;
      if(x < cumsum)
        {
          chosen = candidates[c].index;
          break;
        }
    }

  if(chosen < 0)
    chosen = candidates[ncand - 1].index; /* floating-point edge case fallback */

  sidm_apply_kick(i, chosen);
}

/*! \brief Driver: evaluates scattering for every DM particle at exactly
 *  its own native timebin.
 *
 *  Called from do_gravity_hydro.c's per-hierarchical-timebin loop --
 *  see the call site comment for why this filters to
 *  P[i].TimeBinGrav == timebin rather than using the cumulative active
 *  list gravity's own kick uses.
 */
void sidm_scatter(int timebin)
{
  double dt_i = sidm_dt_code_units(timebin);

  if(dt_i <= 0)
    return;

  for(int idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      int i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type != 1)
        continue;

      if(P[i].TimeBinGrav != timebin)
        continue; /* not this particle's own native step -- see file header */

      sidm_scatter_evaluate(i, dt_i);
    }

  /* TEMPORARY: MPI-reduced scatter-count and conservation-error summary.
   * Remove once the pipeline is trusted / real diagnostics exist. */
  long long global_count;
  double    global_max_p_err, global_max_ke_err;
  MPI_Allreduce(&sidm_scatter_count_local, &global_count, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&sidm_max_momentum_error, &global_max_p_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&sidm_max_energy_error, &global_max_ke_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  if(global_count > 0)
    mpi_printf("SIDM_SCATTER: timebin=%d cumulative scatters=%lld  max_p_rel_err=%.3e  max_ke_rel_err=%.3e\n", timebin,
               global_count, global_max_p_err, global_max_ke_err);
}

#endif /* #ifdef SIDM */
