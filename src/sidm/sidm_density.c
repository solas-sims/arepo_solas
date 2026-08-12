/*!
 * \file        src/sidm/sidm_density.c
 * \brief       v1 elastic SIDM: local DM density + velocity dispersion via
 *              tree-reuse on the existing gravity tree. Filters Type==1 (DM)
 *              at the leaf level (option b -- no dedicated DM mass moment
 *              yet). Driver/Hsml-iteration structure follows star_density()
 *              in src/stars/star_density_sph.c; per-particle data_in/out and
 *              kernel_local/kernel_imported follow the same file's pattern.
 *              Density uses the standard cubic spline (M4) kernel via the
 *              codebase's own KERNEL_COEFF_* constants (see
 *              src/subfind/subfind_density.c for the matching precedent),
 *              not a top-hat -- matches convention used elsewhere and avoids
 *              the top-hat's known small-N negative density bias.
 *              Neighbour-search overlap test (node center +/- 0.5*len vs
 *              search radius Hsml) is NOT the gravity MAC -- deliberately
 *              different criterion, see Day 1/2 notes.
 */

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"
#include "../domain/domain.h"
#include "sidm.h"
#include "sidm_tree.h"

#ifdef SIDM

int               NumDM;
DM_Particle_Data *DMSP;

static int sidm_density_evaluate(int target, int mode, int threadid);
static int sidm_density_isactive(int n);

static MyFloat *SidmNgbs;

typedef struct
{
  MyDouble Pos[3];
  MyFloat Hsml;
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

static void particle2in(data_in *in, int i, int firstnode)
{
  for(int k = 0; k < 3; k++)
    in->Pos[k] = P[i].Pos[k];
  in->Hsml       = DMPS(i).SidmHsml;
  in->Firstnode  = firstnode;
}

typedef struct
{
  MyFloat Density;
  MyFloat VxSum, VySum, VzSum; /*!< mass*kernel-weighted first moments of velocity -- raw
                                 * sums are additive across local+imported contributions;
                                 * a sqrt'd dispersion is NOT, which was the bug this
                                 * replaces (see file header for the full explanation) */
  MyFloat V2Sum;               /*!< mass*kernel-weighted sum of |v|^2 */
  MyFloat Ngb;
} data_out;

static data_out *DataResult, *DataOut;

/*! Persistent per-particle accumulators for the raw velocity moments,
 *  mirroring SidmNgbs[]'s role: reset by the LOCAL (mode==MODE_LOCAL_PARTICLES)
 *  call each Hsml-test iteration, then accumulated onto by however many
 *  imported/cross-task contributions arrive for that same iteration. */
static MyFloat *SidmVxSum, *SidmVySum, *SidmVzSum, *SidmV2Sum;

static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)
    {
      DMPS(i).SidmDensity = out->Density;
      DMPS(i).SidmNumNgb  = (int)out->Ngb;
      SidmNgbs[i]      = out->Ngb;
      SidmVxSum[i]     = out->VxSum;
      SidmVySum[i]     = out->VySum;
      SidmVzSum[i]     = out->VzSum;
      SidmV2Sum[i]     = out->V2Sum;
    }
  else
    {
      DMPS(i).SidmDensity += out->Density;
      DMPS(i).SidmNumNgb  += (int)out->Ngb;
      SidmNgbs[i]       += out->Ngb;
      SidmVxSum[i]      += out->VxSum;
      SidmVySum[i]      += out->VySum;
      SidmVzSum[i]      += out->VzSum;
      SidmV2Sum[i]      += out->V2Sum;
    }
}

#include "../utils/generic_comm_helpers2.h"

static int sidm_density_isactive(int n)
{
  if(P[n].Type != 1) /* DM only, v1 */
    return 0;
  return 1;
}

static void kernel_local(void)
{
  int idx, i, j;
  int threadid = get_thread_num();

  for(j = 0; j < NTask; j++)
    Thread[threadid].Exportflag[j] = -1;

  while(1)
    {
      if(Thread[threadid].ExportSpace < MinSpace)
        break;

      idx = NextParticle++;

      if(idx >= TimeBinsGravity.NActiveParticles)
        break;

      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(sidm_density_isactive(i))
        sidm_density_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
    }
}

static void kernel_imported(void)
{
  int i, cnt = 0;
  int threadid = get_thread_num();

  while(1)
    {
      i = cnt++;

      if(i >= Nimport)
        break;

      sidm_density_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
    }
}

/*! \brief Neighbour-search tree walk for SIDM density estimate.
 *
 *  Walks the dedicated SidmTree_Nodes[] (src/sidm/sidm_tree.c) rather
 *  than the gravity tree -- see the FLAG_FULL_TREE non-convergence bug
 *  this replaces. AABB overlap test (range_min/range_max vs search
 *  sphere) mirrors ngb_treefind_variable_threads() in
 *  src/ngbtree/ngbtree_walk.c exactly, including periodic wrap via
 *  NGB_PERIODIC_LONG_X/Y/Z (degrades to plain fabs() under
 *  GRAVITY_NOT_PERIODIC, same macro either way).
 *
 *  No separate "imported point" branch: unlike the old gravity-tree
 *  walk (which read foreign particles from a Tree_Points import
 *  buffer), this tree -- like ngbtree -- has no such buffer. A remote
 *  search is computed on the task that actually owns the local data
 *  and shipped back as a small result via DataResult/out2particle,
 *  matching ngbtree's (more memory-efficient) convention.
 *
 *  drift_particle() is called here defensively (matching
 *  ngb_treefind_variable_threads's own pattern) even though, given
 *  this tree is fully rebuilt immediately before every walk (see
 *  sidm_density()'s driver), every DM particle should already be
 *  drifted current from sidm_treebuild_construct() itself. No
 *  drift_node() equivalent is needed or present -- see sidm_tree.h's
 *  header comment for why: no node-level velocity extrapolation is
 *  applicable to a tree that's rebuilt fresh each call rather than
 *  lazily aged between builds.
 */
static int sidm_density_evaluate(int target, int mode, int threadid)
{
  int numnodes, *firstnode;
  data_in local, *target_data;
  data_out out, *target_result;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data   = &local;
      target_result = &out;
      numnodes      = 1;
      firstnode     = NULL;
    }
  else
    {
      target_data   = &DataGet[target];
      target_result = &DataResult[target];
      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  double pos_x = target_data->Pos[0];
  double pos_y = target_data->Pos[1];
  double pos_z = target_data->Pos[2];
  double h     = target_data->Hsml;
  double h2    = h * h;
  double hinv  = 1.0 / h;
  double hinv3 = hinv * hinv * hinv;

  double search_min[3], search_max[3], search_max_Lsub[3], search_min_Ladd[3];
#if !defined(GRAVITY_NOT_PERIODIC)
  double xtmp, ytmp, ztmp;
#endif

  search_min[0] = pos_x - 1.001 * h;
  search_min[1] = pos_y - 1.001 * h;
  search_min[2] = pos_z - 1.001 * h;
  search_max[0] = pos_x + 1.001 * h;
  search_max[1] = pos_y + 1.001 * h;
  search_max[2] = pos_z + 1.001 * h;

#if !defined(GRAVITY_NOT_PERIODIC)
  search_max_Lsub[0] = search_max[0] - boxSize_X;
  search_max_Lsub[1] = search_max[1] - boxSize_Y;
  search_max_Lsub[2] = search_max[2] - boxSize_Z;

  search_min_Ladd[0] = search_min[0] + boxSize_X;
  search_min_Ladd[1] = search_min[1] + boxSize_Y;
  search_min_Ladd[2] = search_min[2] + boxSize_Z;
#else
  /* Non-periodic: neutralize the wrap-around terms entirely rather than
   * relying on boxSize_X/Y/Z being sane for this purpose -- BoxSize=0
   * (the standard, CORRECT convention for a genuinely non-periodic
   * setup, e.g. an isolated halo IC) makes the wrap terms degenerate
   * (search_max_Lsub=search_max, search_min_Ladd=search_min), which
   * causes every pruning test below to become impossible to satisfy --
   * i.e. the tree never rejects ANY node, degrading every query from an
   * O(log N) walk to a full O(N) scan. For ~1e5 particles this is
   * O(N^2) total work: technically still making progress (pinned at
   * ~100% CPU) but never finishing in practical time -- indistinguishable
   * from a hang. Setting these to +/-infinity makes each two-clause
   * pruning condition below collapse to a plain single-sided overlap
   * test (see the arithmetic worked out in the fix commit/PR
   * description), which is what a non-periodic search should be doing
   * anyway. */
  search_max_Lsub[0] = search_max_Lsub[1] = search_max_Lsub[2] = -INFINITY;
  search_min_Ladd[0] = search_min_Ladd[1] = search_min_Ladd[2] = INFINITY;
#endif /* #if !defined(GRAVITY_NOT_PERIODIC) #else */

  double density_sum = 0.0;
  double vx_sum = 0.0, vy_sum = 0.0, vz_sum = 0.0;
  double v2_sum       = 0.0;
  int    numngb       = 0;

  for(int k = 0; k < numnodes; k++)
    {
      int no;

      if(mode == MODE_LOCAL_PARTICLES)
        no = SidmTree_MaxPart;
      else
        {
          no = firstnode[k];
          no = SidmTree_Nodes[no].u.d.nextnode;
        }

      while(no >= 0)
        {
          if(no < SidmTree_MaxPart) /* single particle */
            {
              int p = no;
              no    = SidmTree_Nextnode[no];

              if(P[p].Type != 1)
                continue;

              if(P[p].Ti_Current != All.Ti_Current)
                drift_particle(p, All.Ti_Current);

              double dx = GRAVITY_NEAREST_X(P[p].Pos[0] - pos_x);
              if(fabs(dx) > h)
                continue;
              double dy = GRAVITY_NEAREST_Y(P[p].Pos[1] - pos_y);
              if(fabs(dy) > h)
                continue;
              double dz = GRAVITY_NEAREST_Z(P[p].Pos[2] - pos_z);
              if(fabs(dz) > h)
                continue;

              double r2 = dx * dx + dy * dy + dz * dz;
              if(r2 > h2)
                continue;

              double r = sqrt(r2);
              double u = r * hinv;
              double wk;

              if(u < 0.5)
                wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
              else
                wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

              density_sum += P[p].Mass * wk;
              vx_sum += P[p].Mass * wk * P[p].Vel[0];
              vy_sum += P[p].Mass * wk * P[p].Vel[1];
              vz_sum += P[p].Mass * wk * P[p].Vel[2];
              v2_sum += P[p].Mass * wk *
                        (P[p].Vel[0] * P[p].Vel[0] + P[p].Vel[1] * P[p].Vel[1] + P[p].Vel[2] * P[p].Vel[2]);
              numngb++;
            }
          else if(no < SidmTree_MaxPart + SidmTree_MaxNodes) /* internal node */
            {
              struct SidmNODE *current = &SidmTree_Nodes[no];

              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no < SidmTree_FirstNonTopLevelNode)
                    break;
                }

              no = current->u.d.sibling; /* assume we can discard, restore below if not */

              if(search_min[0] > current->u.d.range_max[0] && search_max_Lsub[0] < current->u.d.range_min[0])
                continue;
              if(search_min_Ladd[0] > current->u.d.range_max[0] && search_max[0] < current->u.d.range_min[0])
                continue;
              if(search_min[1] > current->u.d.range_max[1] && search_max_Lsub[1] < current->u.d.range_min[1])
                continue;
              if(search_min_Ladd[1] > current->u.d.range_max[1] && search_max[1] < current->u.d.range_min[1])
                continue;
              if(search_min[2] > current->u.d.range_max[2] && search_max_Lsub[2] < current->u.d.range_min[2])
                continue;
              if(search_min_Ladd[2] > current->u.d.range_max[2] && search_max[2] < current->u.d.range_min[2])
                continue;

              no = current->u.d.nextnode; /* node might contain neighbours: descend */
            }
          else /* pseudo-particle: export to owning task */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("SIDM_DENSITY: MODE_IMPORTED_PARTICLES should not reach a pseudo-particle here");

              if(target >= 0)
                sidm_treefind_export_node_threads(no, target, threadid);

              no = SidmTree_Nextnode[no - SidmTree_MaxNodes];
            }
        }
    }

  target_result->Density = density_sum;
  target_result->Ngb     = numngb;
  target_result->VxSum   = vx_sum;
  target_result->VySum   = vy_sum;
  target_result->VzSum   = vz_sum;
  target_result->V2Sum   = v2_sum;

  /* MODE_IMPORTED_PARTICLES: do NOT call out2particle here. `target` in
   * this branch is a position in the import scratch buffer on the
   * SERVING task, not a real local particle index -- it has no
   * relation to the particle that actually requested this contribution.
   * Writing to P[target]/SidmNgbs[target] here was the bug: it wrote
   * this remote contribution onto whatever unrelated local particle
   * happened to sit at that arbitrary index. generic_comm_helpers2.h
   * already calls out2particle(&DataOut[off], target, MODE_IMPORTED_PARTICLES)
   * itself, on the ORIGINATING task, after shipping DataResult back over
   * MPI -- at which point `target` correctly means the real local index
   * again. Populating target_result (DataResult[target]) above is all
   * this branch needs to do.
   */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(target_result, target, MODE_LOCAL_PARTICLES);

  return numngb;
}

/*! \brief Driver: Hsml-iteration density pass for all active DM particles.
 *
 *  Structure mirrors star_density() in src/stars/star_density_sph.c
 *  (Left/Right bisection to converge neighbour count).
 */
void sidm_density(void)
{
  MyFloat *Left, *Right;
  int idx, i, npleft, iter = 0;
  long long ntot;

  /* sidm_treebuild() is now called explicitly from run.c, tied to the
   * same domain-decomposition cadence as ngb_treebuild(), immediately
   * before this function -- not from inside here. That fixed the
   * STOP!/domain-mismatch crash: the tree-build assertion (ported from
   * ngb_treebuild_construct) assumes particle-to-task assignment is
   * fresh, which is only true right after domain decomposition, not at
   * the more frequent FLAG_FULL_TREE cadence this used to run at. See
   * the corresponding discussion for why walk-cadence and rebuild-
   * cadence are now tied together rather than decoupled (node-level
   * drift tracking, which this tree deliberately doesn't have).
   */

  SidmNgbs = (MyFloat *)mymalloc("SidmNgbs", NumPart * sizeof(MyFloat));
  Left     = (MyFloat *)mymalloc("Left", NumPart * sizeof(MyFloat));
  Right    = (MyFloat *)mymalloc("Right", NumPart * sizeof(MyFloat));
  SidmVxSum = (MyFloat *)mymalloc("SidmVxSum", NumPart * sizeof(MyFloat));
  SidmVySum = (MyFloat *)mymalloc("SidmVySum", NumPart * sizeof(MyFloat));
  SidmVzSum = (MyFloat *)mymalloc("SidmVzSum", NumPart * sizeof(MyFloat));
  SidmV2Sum = (MyFloat *)mymalloc("SidmV2Sum", NumPart * sizeof(MyFloat));

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0 || !sidm_density_isactive(i))
        continue;

      Left[i] = Right[i] = 0;

      if(DMPS(i).SidmHsml <= 0)
        DMPS(i).SidmHsml = All.ForceSoftening[P[i].SofteningType];
    }

  generic_set_MaxNexport();

  do
    {
      generic_comm_pattern(TimeBinsGravity.NActiveParticles, kernel_local, kernel_imported);

      for(idx = 0, npleft = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
        {
          i = TimeBinsGravity.ActiveParticleList[idx];
          if(i < 0 || !sidm_density_isactive(i))
            continue;

          if(SidmNgbs[i] < (All.SidmDesNumNgb - All.SidmDesNumNgbDev) ||
             SidmNgbs[i] > (All.SidmDesNumNgb + All.SidmDesNumNgbDev))
            {
              npleft++;

              if(Left[i] > 0 && Right[i] > 0)
                if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
                  {
                    npleft--;
                    continue;
                  }

              if(SidmNgbs[i] < (All.SidmDesNumNgb - All.SidmDesNumNgbDev))
                Left[i] = dmax(DMPS(i).SidmHsml, Left[i]);
              else
                {
                  if(Right[i] != 0)
                    {
                      if(DMPS(i).SidmHsml < Right[i])
                        Right[i] = DMPS(i).SidmHsml;
                    }
                  else
                    Right[i] = DMPS(i).SidmHsml;
                }

              if(Right[i] > 0 && Left[i] > 0)
                DMPS(i).SidmHsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
              else
                {
                  if(Right[i] == 0 && Left[i] == 0)
                    terminate("SIDM_DENSITY: should not occur");

                  if(Right[i] == 0 && Left[i] > 0)
                    DMPS(i).SidmHsml *= 1.26;

                  if(Right[i] > 0 && Left[i] == 0)
                    DMPS(i).SidmHsml /= 1.26;
                }
            }
        }

      sumup_large_ints(1, &npleft, &ntot);

      if(ntot > 0)
        {
          iter++;
          mpi_printf("SIDM_DENSITY: ngb iteration %3d: need to repeat for %12lld particles.\n", iter, ntot);

          if(iter > MAXITER)
            terminate("SIDM_DENSITY: failed to converge in neighbour iteration\n");
        }
    }
  while(ntot > 0);

  /* Final VelDisp computation: mean-subtracted, using the LAST
   * (converged) iteration's accumulated raw moments. This is the fix --
   * previously VelDisp was computed as sqrt(<v^2>/3) with no mean
   * subtraction (biased by any bulk/coherent flow, not just genuine
   * local dispersion), AND accumulated across local+imported
   * contributions by summing sqrt() results together, which is not
   * mathematically valid (see file header). Both are fixed by shipping
   * raw first/second moments through out2particle (additive, unlike
   * sqrt'd values) and only taking the sqrt once, here, after all
   * contributions for the final converged pass are in.
   *
   * var_3d = <v^2> - <v>^2 is the trace of the velocity dispersion
   * tensor (sigma_x^2+sigma_y^2+sigma_z^2); dividing by 3 and taking
   * the sqrt gives the 1D dispersion convention already used
   * elsewhere in this code. Clamped to >=0 to guard against tiny
   * negative values from floating-point rounding in near-zero-
   * dispersion regions.
   */
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0 || !sidm_density_isactive(i))
        continue;

      if(DMPS(i).SidmDensity > 0)
        {
          double mean_vx = SidmVxSum[i] / DMPS(i).SidmDensity;
          double mean_vy = SidmVySum[i] / DMPS(i).SidmDensity;
          double mean_vz = SidmVzSum[i] / DMPS(i).SidmDensity;
          double mean_v2 = mean_vx * mean_vx + mean_vy * mean_vy + mean_vz * mean_vz;

          double var_3d = SidmV2Sum[i] / DMPS(i).SidmDensity - mean_v2;
          var_3d         = fmax(0.0, var_3d);

          DMPS(i).SidmVelDisp = sqrt(var_3d / 3.0);
        }
      else
        DMPS(i).SidmVelDisp = 0.0;
    }

  /* --- TEMPORARY diagnostic block: no I/O wiring exists yet for these
   * fields, so this is the only way to see actual values right now.
   * Remove once SidmDensity/SidmHsml/SidmNumNgb are in snapshot output. */
  {
    double loc_min_dens = 1e300, loc_max_dens = -1e300, loc_sum_dens = 0;
    double loc_min_hsml = 1e300, loc_max_hsml = -1e300, loc_sum_hsml = 0;
    double loc_min_ngb  = 1e300, loc_max_ngb  = -1e300, loc_sum_ngb  = 0;
    long long loc_count = 0, loc_nan_count = 0, loc_zero_ngb_count = 0;

    for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
      {
        i = TimeBinsGravity.ActiveParticleList[idx];
        if(i < 0 || !sidm_density_isactive(i))
          continue;

        double d = DMPS(i).SidmDensity;
        double h = DMPS(i).SidmHsml;
        double n = (double)DMPS(i).SidmNumNgb;

        if(isnan(d) || isnan(h) || isinf(d) || isinf(h))
          {
            loc_nan_count++;
            continue;
          }

        if(n == 0)
          loc_zero_ngb_count++;

        loc_min_dens = fmin(loc_min_dens, d);
        loc_max_dens = fmax(loc_max_dens, d);
        loc_sum_dens += d;

        loc_min_hsml = fmin(loc_min_hsml, h);
        loc_max_hsml = fmax(loc_max_hsml, h);
        loc_sum_hsml += h;

        loc_min_ngb = fmin(loc_min_ngb, n);
        loc_max_ngb = fmax(loc_max_ngb, n);
        loc_sum_ngb += n;

        loc_count++;
      }

    double glob_min_dens, glob_max_dens, glob_sum_dens;
    double glob_min_hsml, glob_max_hsml, glob_sum_hsml;
    double glob_min_ngb, glob_max_ngb, glob_sum_ngb;
    long long glob_count, glob_nan_count, glob_zero_ngb_count;

    MPI_Allreduce(&loc_min_dens, &glob_min_dens, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_max_dens, &glob_max_dens, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_sum_dens, &glob_sum_dens, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_min_hsml, &glob_min_hsml, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_max_hsml, &glob_max_hsml, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_sum_hsml, &glob_sum_hsml, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_min_ngb, &glob_min_ngb, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_max_ngb, &glob_max_ngb, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_sum_ngb, &glob_sum_ngb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_count, &glob_count, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_nan_count, &glob_nan_count, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_zero_ngb_count, &glob_zero_ngb_count, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    mpi_printf("SIDM_DENSITY diag: N_active=%lld  NaN/Inf=%lld  zero_ngb=%lld\n", glob_count, glob_nan_count,
               glob_zero_ngb_count);
    if(glob_count > 0)
      {
        mpi_printf("SIDM_DENSITY diag: Density  min=%.4e  max=%.4e  mean=%.4e\n", glob_min_dens, glob_max_dens,
                   glob_sum_dens / glob_count);
        mpi_printf("SIDM_DENSITY diag: Hsml     min=%.4e  max=%.4e  mean=%.4e\n", glob_min_hsml, glob_max_hsml,
                   glob_sum_hsml / glob_count);
        mpi_printf("SIDM_DENSITY diag: NumNgb   min=%.1f  max=%.1f  mean=%.2f  (target=%.1f +/- %.1f)\n", glob_min_ngb,
                   glob_max_ngb, glob_sum_ngb / glob_count, All.SidmDesNumNgb, All.SidmDesNumNgbDev);
      }
  }

  myfree(SidmV2Sum);
  myfree(SidmVzSum);
  myfree(SidmVySum);
  myfree(SidmVxSum);
  myfree(Right);
  myfree(Left);
  myfree(SidmNgbs);

  mpi_printf("SIDM_DENSITY: done.\n");
}

#endif /* #ifdef SIDM */
