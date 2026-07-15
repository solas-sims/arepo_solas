/*!
 * \file        src/sidm/sidm_scatter.c
 * \brief       v1 elastic SIDM: Monte Carlo pairwise scattering/kick,
 *              following Vogelsberger, Zavala & Loeb (2012)'s algorithm
 *              (as described in e.g. Valdarnini 2023, arXiv:2309.10374,
 *              Eq. 18-20 -- not derived independently).
 *
 *              CROSS-TASK PARTNERS: candidates are now gathered from
 *              both the local tree AND, via the pseudo-particle branch
 *              (generic_comm_pattern, mirroring sidm_density.c's own
 *              structure), from whichever remote tasks own a
 *              particle's neighbouring domain. This replaces the
 *              earlier local-partners-only v1 restriction. A remote
 *              task self-reports its own identity in each response
 *              (DataOut.remote_task = ThisTask on the serving side) --
 *              the originating task cannot otherwise determine which
 *              task a given response came from purely from the
 *              generic_comm_helpers2.h framework's own bookkeeping.
 *
 *              KNOWN, ACCEPTED APPROXIMATION (deliberate, discussed):
 *              if the selected partner is remote, the kick is computed
 *              using that partner's velocity AS SAMPLED during
 *              candidate gathering. Between that sampling and the
 *              kick-delivery exchange later in this same call, the
 *              remote particle could in principle have already been
 *              kicked by a *different*, purely-local scatter on its
 *              own task. This is accepted as timestep-bounded and in
 *              the same spirit as the pre-existing (and never flagged
 *              as a bug) possibility that a purely local particle gets
 *              selected by two different lower-ID initiators in one
 *              pass -- both are rare, both are bounded by the SIDM
 *              timestep criterion keeping per-step scattering
 *              probability small, and a fully rigorous fix (atomic
 *              remote claims + retry) was judged not worth the
 *              complexity for v1.
 *
 *              LOWER-ID-INITIATES: for a pair (i,j), only the particle
 *              with the smaller ID evaluates and (if triggered)
 *              applies the scatter. This deliberately REPLACES
 *              Vogelsberger's P_i = sum(P_ij)/2 convention (the /2
 *              there corrects for each pair being evaluated twice,
 *              once from each side) -- since lower-ID-initiates
 *              structurally prevents double-evaluation in the first
 *              place, there is nothing left for a /2 to correct, so it
 *              is dropped. This is a derived consequence of the
 *              parallel-implementation choice, not a verbatim copy of
 *              the source algorithm. The ID comparison for a REMOTE
 *              candidate uses its real particle ID (shipped back
 *              during gathering), same rule as for a local candidate.
 *
 *              Called once per particle at EXACTLY its own native
 *              timebin (P[i].TimeBinGrav == timebin), not from the
 *              cumulative hierarchical active-list used by gravity's
 *              own kick loop -- that list includes a particle on every
 *              coarser timebin's pass too, which would re-roll its
 *              scatter dice multiple times per actual step. See the
 *              call site in do_gravity_hydro.c for where this matters.
 *
 *              P[] ARRAY STABILITY DURING THIS CALL: the kick-delivery
 *              exchange below relies on a remote P[] index sampled
 *              during candidate gathering still being valid by the
 *              time the poke arrives. This holds because nothing
 *              within a single sidm_scatter() call reorders P[]
 *              (kicks only change velocities, not array position) --
 *              it would NOT hold if domain decomposition or particle
 *              migration ran between gathering and delivery, which is
 *              why both happen back-to-back in this same function.
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

enum
{
  MAX_CAND_TOTAL      = 128, /*!< combined local+remote candidates per particle, well above ~SidmDesNumNgb */
  MAX_REMOTE_RESPONSE = 16   /*!< per-query cap on how many candidates a single remote task reports back */
};

/*! \brief One candidate scattering partner, local or remote. */
typedef struct
{
  int    is_remote;
  int    index;       /*!< local P[] index if !is_remote; remote P[] index on remote_task if is_remote */
  int    remote_task; /*!< only meaningful if is_remote */
  double r;
  double p_ij;
  double vel[3]; /*!< only populated/used for remote candidates -- local candidates re-read P[index].Vel live */
} sidm_candidate;

static sidm_candidate *SidmCandBuf;   /*!< flattened [NumPart][MAX_CAND_TOTAL] */
static int            *SidmCandCount; /*!< per-particle count of valid entries in SidmCandBuf */

static double sidm_scatter_current_dt_i; /*!< set once per sidm_scatter() call */
static int    sidm_scatter_current_timebin;

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
 *  timebin, in code time units. See file header discussion in the
 *  original local-only version for why this must go through
 *  get_time_difference_in_Gyr() rather than a raw timebin multiply.
 */
static double sidm_dt_code_units(int timebin)
{
  integertime ti_step = timebin ? (((integertime)1) << timebin) : 0;
  integertime tend    = All.Ti_begstep[timebin];
  integertime tstart  = tend - ti_step;

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

/*! Live conservation diagnostic -- same role as the local-only version,
 *  now shared by both the local and cross-task kick paths. */
static long long sidm_scatter_count_local     = 0;
static long long sidm_scatter_count_cross_task = 0;
static double    sidm_max_momentum_error      = 0.0;
static double    sidm_max_energy_error        = 0.0;

static void sidm_check_conservation(int i, double p_before[3], double ke_before, double p_after[3], double ke_after,
                                     long long other_id)
{
  double p_mag_before = sqrt(p_before[0] * p_before[0] + p_before[1] * p_before[1] + p_before[2] * p_before[2]);
  double dp[3]  = {p_after[0] - p_before[0], p_after[1] - p_before[1], p_after[2] - p_before[2]};
  double dp_mag = sqrt(dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2]);
  double p_rel_err  = (p_mag_before > 0) ? dp_mag / p_mag_before : dp_mag;
  double ke_rel_err = (ke_before > 0) ? fabs(ke_after - ke_before) / ke_before : fabs(ke_after - ke_before);

  if(p_rel_err > sidm_max_momentum_error)
    sidm_max_momentum_error = p_rel_err;
  if(ke_rel_err > sidm_max_energy_error)
    sidm_max_energy_error = ke_rel_err;

  if(p_rel_err > 1e-10 || ke_rel_err > 1e-10)
    printf("SIDM_SCATTER conservation WARNING: task=%d IDs=(%lld,%lld) p_rel_err=%.3e ke_rel_err=%.3e\n", ThisTask,
           (long long)P[i].ID, other_id, p_rel_err, ke_rel_err);

  sidm_scatter_count_local++;
}

/*! \brief Applies the isotropic elastic two-body kick when BOTH particles
 *  are local (Vogelsberger et al. 2012, Eq. 20).
 *
 *  NOTE: assumes m_i == m_j (matches the source paper's own stated
 *  assumption of equal-mass DM particles, and this codebase's
 *  equal-mass DM box) -- unchanged from the original local-only
 *  version.
 */
static void sidm_apply_kick_local(int i, int j)
{
  double p_before[3], ke_before;
  for(int k = 0; k < 3; k++)
    p_before[k] = P[i].Mass * P[i].Vel[k] + P[j].Mass * P[j].Vel[k];
  ke_before = 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]) +
              0.5 * P[j].Mass * (P[j].Vel[0] * P[j].Vel[0] + P[j].Vel[1] * P[j].Vel[1] + P[j].Vel[2] * P[j].Vel[2]);

  double vij[3], V[3], e[3];
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

  double p_after[3], ke_after;
  for(int k = 0; k < 3; k++)
    p_after[k] = P[i].Mass * P[i].Vel[k] + P[j].Mass * P[j].Vel[k];
  ke_after = 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]) +
             0.5 * P[j].Mass * (P[j].Vel[0] * P[j].Vel[0] + P[j].Vel[1] * P[j].Vel[1] + P[j].Vel[2] * P[j].Vel[2]);

  sidm_check_conservation(i, p_before, ke_before, p_after, ke_after, (long long)P[j].ID);

  DMPS(i).SidmLastScatterTime = All.Ti_Current;
  DMPS(j).SidmLastScatterTime = All.Ti_Current;
  DMPS(i).SidmScatterFlag     = 1;
  DMPS(j).SidmScatterFlag     = 1;
}

/*! \brief Pending instruction to a remote task: "set this local particle's
 *  velocity to this value, and mark it as scattered" -- the one-way
 *  delivery half of a cross-task kick. Batched and sent once per
 *  destination task, after all local particles this call have been
 *  processed, rather than one message per kick. */
typedef struct
{
  int         remote_index;
  MyFloat     new_vel[3];
  integertime scatter_time;
} sidm_kick_poke;

/*! \brief Computes both new velocities for a local-remote pair, applies
 *  the local half immediately, and queues the remote half into
 *  per-destination-task buffers for the batched delivery exchange.
 *
 *  poke_buf/poke_count are indexed by destination task; each poke_buf[t]
 *  is a fixed-capacity buffer sized generously up front by the caller
 *  (see sidm_scatter()) -- a dynamically-grown-per-append buffer would
 *  violate the LIFO mymalloc discipline this codebase requires.
 */
static void sidm_apply_kick_cross_task(int i, const sidm_candidate *cand, sidm_kick_poke **poke_buf, int *poke_count,
                                        int poke_capacity)
{
  double p_before[3], ke_before;
  for(int k = 0; k < 3; k++)
    p_before[k] = P[i].Mass * P[i].Vel[k] + P[i].Mass * cand->vel[k];
  ke_before = 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]) +
              0.5 * P[i].Mass * (cand->vel[0] * cand->vel[0] + cand->vel[1] * cand->vel[1] + cand->vel[2] * cand->vel[2]);

  double vij[3], V[3], e[3];
  for(int k = 0; k < 3; k++)
    {
      vij[k] = P[i].Vel[k] - cand->vel[k];
      V[k]   = 0.5 * (P[i].Vel[k] + cand->vel[k]);
    }
  double v_rel = sqrt(vij[0] * vij[0] + vij[1] * vij[1] + vij[2] * vij[2]);
  sidm_random_unit_vector(e);

  double new_vel_i[3], new_vel_j[3];
  for(int k = 0; k < 3; k++)
    {
      new_vel_i[k] = V[k] + 0.5 * v_rel * e[k];
      new_vel_j[k] = V[k] - 0.5 * v_rel * e[k];
    }

  double p_after[3], ke_after;
  for(int k = 0; k < 3; k++)
    p_after[k] = P[i].Mass * new_vel_i[k] + P[i].Mass * new_vel_j[k];
  ke_after = 0.5 * P[i].Mass * (new_vel_i[0] * new_vel_i[0] + new_vel_i[1] * new_vel_i[1] + new_vel_i[2] * new_vel_i[2]) +
             0.5 * P[i].Mass * (new_vel_j[0] * new_vel_j[0] + new_vel_j[1] * new_vel_j[1] + new_vel_j[2] * new_vel_j[2]);

  /* Remote particle's real ID isn't carried into sidm_candidate to keep
   * that struct lean -- report -1 rather than fabricate one. */
  sidm_check_conservation(i, p_before, ke_before, p_after, ke_after, -1);

  for(int k = 0; k < 3; k++)
    P[i].Vel[k] = new_vel_i[k];
  DMPS(i).SidmLastScatterTime = All.Ti_Current;
  DMPS(i).SidmScatterFlag     = 1;

  int t = cand->remote_task;
  if(t < 0 || t >= NTask)
    terminate("SIDM_SCATTER: cross-task kick with invalid remote_task=%d -- remote_task was not correctly self-reported\n",
              t);

  sidm_scatter_count_cross_task++;

  if(poke_count[t] >= poke_capacity)
    terminate(
        "SIDM_SCATTER: exceeded per-task pending-poke capacity (%d) -- more remote kicks to a single task this step than "
        "provisioned for\n",
        poke_capacity);

  int slot                       = poke_count[t]++;
  poke_buf[t][slot].remote_index = cand->index;
  poke_buf[t][slot].new_vel[0]   = new_vel_j[0];
  poke_buf[t][slot].new_vel[1]   = new_vel_j[1];
  poke_buf[t][slot].new_vel[2]   = new_vel_j[2];
  poke_buf[t][slot].scatter_time = All.Ti_Current;
}

/*! \brief data_in/data_out for the candidate-gathering comm pattern,
 *  mirroring sidm_density.c's structure exactly. */
typedef struct
{
  MyDouble Pos[3];
  MyFloat  Hsml;
  MyIDType SelfID;
  int      Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

static void particle2in(data_in *in, int i, int firstnode)
{
  for(int k = 0; k < 3; k++)
    in->Pos[k] = P[i].Pos[k];
  in->Hsml      = DMPS(i).SidmHsml;
  in->SelfID    = P[i].ID;
  in->Firstnode = firstnode;
}

typedef struct
{
  int      ncand;
  int      remote_task; /*!< self-reported by the serving task: ThisTask at response time */
  MyIDType id[MAX_REMOTE_RESPONSE];
  MyFloat  vel[MAX_REMOTE_RESPONSE][3];
  MyFloat  r[MAX_REMOTE_RESPONSE];
  int      remote_index[MAX_REMOTE_RESPONSE];
} data_out;

static data_out *DataResult, *DataOut;

static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)
    return; /* local candidates are appended directly into SidmCandBuf during the walk itself, not here */

  /* MODE_IMPORTED_PARTICLES: append this remote task's response onto
   * particle i's candidate list, computing p_ij now (we have i's own
   * Mass/Vel/dt available here on the originating task). */
  for(int c = 0; c < out->ncand; c++)
    {
      if(out->id[c] <= P[i].ID)
        continue; /* defensive: the serving side should already enforce this */

      if(SidmCandCount[i] >= MAX_CAND_TOTAL)
        break; /* silently drop excess -- same acceptable-cap convention as the local walk */

      double vx   = P[i].Vel[0] - out->vel[c][0];
      double vy   = P[i].Vel[1] - out->vel[c][1];
      double vz   = P[i].Vel[2] - out->vel[c][2];
      double v_ij = sqrt(vx * vx + vy * vy + vz * vz);

      double hinv  = 1.0 / DMPS(i).SidmHsml;
      double hinv3 = hinv * hinv * hinv;
      double u     = out->r[c] * hinv;
      double wk;
      if(u < 0.5)
        wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
      else
        wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

      double p_ij = P[i].Mass * wk * All.SidmCrossSection * v_ij * sidm_scatter_current_dt_i;

      int    idx  = SidmCandCount[i];
      size_t base = (size_t)i * MAX_CAND_TOTAL + idx;
      SidmCandBuf[base].is_remote   = 1;
      SidmCandBuf[base].index       = out->remote_index[c];
      SidmCandBuf[base].remote_task = out->remote_task;
      SidmCandBuf[base].r           = out->r[c];
      SidmCandBuf[base].p_ij        = p_ij;
      SidmCandBuf[base].vel[0]      = out->vel[c][0];
      SidmCandBuf[base].vel[1]      = out->vel[c][1];
      SidmCandBuf[base].vel[2]      = out->vel[c][2];
      SidmCandCount[i]++;
    }
}

#include "../utils/generic_comm_helpers2.h"

/*! \brief The shared walk, used both for a particle's own local search
 *  (MODE_LOCAL_PARTICLES) and for serving a remote task's query
 *  (MODE_IMPORTED_PARTICLES).
 *
 *  Local mode: appends local candidates directly into SidmCandBuf, and
 *  now genuinely follows pseudo-particle branches (exporting a query)
 *  instead of skipping them, as the original local-only version did.
 *
 *  Imported mode: walks starting from the given firstnode(s) -- i.e.
 *  this task's own locally-owned portion of the tree -- and packs
 *  matches into the data_out response, self-reporting remote_task =
 *  ThisTask so the originating task knows which task the response is
 *  from.
 */
static int sidm_scatter_walk(int target, int mode, int threadid)
{
  int      numnodes, *firstnode;
  data_in  local, *target_data;
  data_out out, *target_result;

  MyIDType self_id;
  double   pos_x, pos_y, pos_z, h;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data   = &local;
      target_result = &out;
      numnodes      = 1;
      firstnode     = NULL;
      self_id       = P[target].ID;
    }
  else
    {
      target_data   = &DataGet[target];
      target_result = &DataResult[target];
      generic_get_numnodes(target, &numnodes, &firstnode);
      self_id = target_data->SelfID;
    }

  pos_x     = target_data->Pos[0];
  pos_y     = target_data->Pos[1];
  pos_z     = target_data->Pos[2];
  h         = target_data->Hsml;
  double h2 = h * h;

  if(mode == MODE_IMPORTED_PARTICLES)
    {
      target_result->ncand       = 0;
      target_result->remote_task = ThisTask;
    }

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

              if(p == target && mode == MODE_LOCAL_PARTICLES)
                continue;
              if(P[p].Type != 1)
                continue;
              if(P[p].ID <= self_id)
                continue; /* lower-ID-initiates */

              if(P[p].Ti_Current != All.Ti_Current)
                drift_particle(p, All.Ti_Current);

              double dx = P[p].Pos[0] - pos_x;
              double dy = P[p].Pos[1] - pos_y;
              double dz = P[p].Pos[2] - pos_z;
              double r2 = dx * dx + dy * dy + dz * dz;

              if(r2 >= h2)
                continue;

              double r = sqrt(r2);

              if(mode == MODE_LOCAL_PARTICLES)
                {
                  if(SidmCandCount[target] >= MAX_CAND_TOTAL)
                    continue;

                  double hinv  = 1.0 / h;
                  double hinv3 = hinv * hinv * hinv;
                  double u     = r * hinv;
                  double wk;
                  if(u < 0.5)
                    wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                  else
                    wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

                  double vx   = P[target].Vel[0] - P[p].Vel[0];
                  double vy   = P[target].Vel[1] - P[p].Vel[1];
                  double vz   = P[target].Vel[2] - P[p].Vel[2];
                  double v_ij = sqrt(vx * vx + vy * vy + vz * vz);

                  double p_ij = P[target].Mass * wk * All.SidmCrossSection * v_ij * sidm_scatter_current_dt_i;

                  int    idx  = SidmCandCount[target];
                  size_t base = (size_t)target * MAX_CAND_TOTAL + idx;
                  SidmCandBuf[base].is_remote   = 0;
                  SidmCandBuf[base].index       = p;
                  SidmCandBuf[base].remote_task = ThisTask;
                  SidmCandBuf[base].r           = r;
                  SidmCandBuf[base].p_ij        = p_ij;
                  SidmCandCount[target]++;
                }
              else
                {
                  if(target_result->ncand >= MAX_REMOTE_RESPONSE)
                    continue;

                  int c = target_result->ncand;
                  target_result->id[c]           = P[p].ID;
                  target_result->vel[c][0]        = P[p].Vel[0];
                  target_result->vel[c][1]        = P[p].Vel[1];
                  target_result->vel[c][2]        = P[p].Vel[2];
                  target_result->r[c]             = r;
                  target_result->remote_index[c]  = p;
                  target_result->ncand++;
                }
            }
          else if(no < SidmTree_MaxPart + SidmTree_MaxNodes) /* internal node */
            {
              struct SidmNODE *current = &SidmTree_Nodes[no];

              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no < SidmTree_FirstNonTopLevelNode)
                    break;
                }

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
          else /* pseudo-particle */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("SIDM_SCATTER: MODE_IMPORTED_PARTICLES should not reach a pseudo-particle here");

              if(target >= 0)
                sidm_treefind_export_node_threads(no, target, threadid);

              no = SidmTree_Nextnode[no - SidmTree_MaxNodes];
            }
        }
    }

  if(mode == MODE_IMPORTED_PARTICLES)
    out2particle(target_result, target, MODE_IMPORTED_PARTICLES);

  return 0;
}

static void kernel_local(void)
{
  int idx, i;
  int threadid = get_thread_num();

  for(int j = 0; j < NTask; j++)
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
      if(P[i].Type != 1)
        continue;
      if(P[i].TimeBinGrav != sidm_scatter_current_timebin)
        continue;
      if(DMPS(i).SidmHsml <= 0 || DMPS(i).SidmDensity <= 0)
        continue;

      SidmCandCount[i] = 0;
      sidm_scatter_walk(i, MODE_LOCAL_PARTICLES, threadid);
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

      sidm_scatter_walk(i, MODE_IMPORTED_PARTICLES, threadid);
    }
}

/*! \brief Driver: gathers local+remote candidates for every DM particle
 *  at exactly its own native timebin, performs accept/reject + partner
 *  selection, then delivers any remote kicks in one batched exchange
 *  per destination task.
 */
void sidm_scatter(int timebin)
{
  double dt_i = sidm_dt_code_units(timebin);
  if(dt_i <= 0)
    return;

  sidm_scatter_current_dt_i    = dt_i;
  sidm_scatter_current_timebin = timebin;

  SidmCandBuf   = (sidm_candidate *)mymalloc("SidmCandBuf", (size_t)NumPart * MAX_CAND_TOTAL * sizeof(sidm_candidate));
  SidmCandCount = (int *)mymalloc("SidmCandCount", NumPart * sizeof(int));
  memset(SidmCandCount, 0, NumPart * sizeof(int));

  generic_set_MaxNexport();
  generic_comm_pattern(TimeBinsGravity.NActiveParticles, kernel_local, kernel_imported);

  /* Per-destination-task pending-poke buffers. Fixed capacity per task
   * (rather than a dynamically-grown list) to stay within the LIFO
   * mymalloc discipline this codebase requires -- sized to
   * NActiveParticles per task, a safe (if generous) upper bound since
   * at most one cross-task kick can originate per active particle. */
  int              poke_capacity = TimeBinsGravity.NActiveParticles > 0 ? TimeBinsGravity.NActiveParticles : 1;
  sidm_kick_poke **poke_buf      = (sidm_kick_poke **)mymalloc("sidm_poke_buf_ptrs", NTask * sizeof(sidm_kick_poke *));
  int             *poke_count    = (int *)mymalloc("sidm_poke_count", NTask * sizeof(int));
  memset(poke_count, 0, NTask * sizeof(int));
  for(int t = 0; t < NTask; t++)
    poke_buf[t] = (sidm_kick_poke *)mymalloc("sidm_poke_buf_task", poke_capacity * sizeof(sidm_kick_poke));

  for(int idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      int i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0 || P[i].Type != 1 || P[i].TimeBinGrav != timebin)
        continue;
      if(DMPS(i).SidmHsml <= 0 || DMPS(i).SidmDensity <= 0)
        continue;

      int ncand = SidmCandCount[i];
      if(ncand == 0)
        continue;

      sidm_candidate *cand = &SidmCandBuf[(size_t)i * MAX_CAND_TOTAL];

      double P_i = 0.0;
      for(int c = 0; c < ncand; c++)
        P_i += cand[c].p_ij;

      double x = get_random_number_aux();
      if(x >= P_i)
        continue; /* no scatter this step */

      qsort(cand, ncand, sizeof(sidm_candidate), sidm_candidate_compare);

      double cumsum = 0.0;
      int    chosen = -1;
      for(int c = 0; c < ncand; c++)
        {
          cumsum += cand[c].p_ij;
          if(x < cumsum)
            {
              chosen = c;
              break;
            }
        }
      if(chosen < 0)
        chosen = ncand - 1;

      if(!cand[chosen].is_remote)
        sidm_apply_kick_local(i, cand[chosen].index);
      else
        sidm_apply_kick_cross_task(i, &cand[chosen], poke_buf, poke_count, poke_capacity);
    }

  /* Deliver the batched pokes: a plain Sendrecv loop over tasks is
   * sufficient here (not an Alltoallv) since payload sizes are small
   * and this isn't a hot path relative to the tree walk above. */
  for(int t = 0; t < NTask; t++)
    {
      if(t == ThisTask)
        continue;

      int nsend = poke_count[t];
      int nrecv;
      MPI_Sendrecv(&nsend, 1, MPI_INT, t, TAG_N, &nrecv, 1, MPI_INT, t, TAG_N, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      sidm_kick_poke *recv_buf = NULL;
      if(nrecv > 0)
        recv_buf = (sidm_kick_poke *)mymalloc("sidm_poke_recv", nrecv * sizeof(sidm_kick_poke));

      MPI_Sendrecv(poke_buf[t], nsend * sizeof(sidm_kick_poke), MPI_BYTE, t, TAG_DMDATA, recv_buf,
                   nrecv * sizeof(sidm_kick_poke), MPI_BYTE, t, TAG_DMDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      for(int p = 0; p < nrecv; p++)
        {
          int ri = recv_buf[p].remote_index;
          for(int k = 0; k < 3; k++)
            P[ri].Vel[k] = recv_buf[p].new_vel[k];
          DMPS(ri).SidmLastScatterTime = recv_buf[p].scatter_time;
          DMPS(ri).SidmScatterFlag     = 1;
        }

      if(nrecv > 0)
        myfree(recv_buf);
    }

  /* Any pokes destined for ThisTask itself (poke_buf[ThisTask]) --
   * shouldn't occur in practice, since a cross-task kick by definition
   * targets a remote_task != ThisTask (a same-task partner is always
   * found and applied via the local-candidate path instead), but guard
   * against it explicitly rather than silently drop it. */
  if(poke_count[ThisTask] > 0)
    terminate("SIDM_SCATTER: %d pending poke(s) targeted ThisTask itself -- should be unreachable\n",
              poke_count[ThisTask]);

  for(int t = NTask - 1; t >= 0; t--) /* reverse order: matches allocation order for LIFO myfree */
    myfree(poke_buf[t]);
  myfree(poke_count);
  myfree(poke_buf);

  myfree(SidmCandCount);
  myfree(SidmCandBuf);

  /* MPI-reduced scatter-count and conservation-error summary. */
  long long global_count, global_cross_task_count;
  double    global_max_p_err, global_max_ke_err;
  MPI_Allreduce(&sidm_scatter_count_local, &global_count, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&sidm_scatter_count_cross_task, &global_cross_task_count, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&sidm_max_momentum_error, &global_max_p_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&sidm_max_energy_error, &global_max_ke_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  if(global_count > 0)
    mpi_printf("SIDM_SCATTER: timebin=%d cumulative scatters=%lld (cross_task=%lld)  max_p_rel_err=%.3e  max_ke_rel_err=%.3e\n",
               timebin, global_count, global_cross_task_count, global_max_p_err, global_max_ke_err);
}

#endif /* #ifdef SIDM */
