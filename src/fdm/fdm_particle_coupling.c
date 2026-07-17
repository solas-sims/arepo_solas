/*!
 * \file        src/fdm/fdm_particle_coupling.c
 * \brief       Interpolates FDM_Potential and FDM_ForceX/Y/Z (the mesh
 *              quantities computed in fdm_poisson.c/fdm_gradient.c) to
 *              star particle positions.
 *
 *              This is the mesh->star half of Phase 2a's coupling; the
 *              star->mesh half (feeding stellar mass density into the
 *              Poisson solve's source term) and hooking the resulting
 *              force into Arepo's own gravity accumulation are
 *              separate, later pieces.
 *
 *              COMMUNICATION PATTERN: deliberately NOT the tree-based
 *              generic_comm_pattern (src/utils/generic_comm_helpers2.h)
 *              -- that machinery is built around NodeList/Firstnode
 *              tree-walk bookkeeping for particles that need to SEARCH
 *              which tree nodes (local or remote) to visit. We have no
 *              search to do: a star's position directly determines
 *              which task(s) own the relevant FDM mesh slab(s), via
 *              FDM_plan.slab_to_task[] -- exactly the same direct-
 *              lookup situation gravity's own PM code is in. Mirrors
 *              pmforce_nonperiodic_uniform_optimized_prepare_density's
 *              /readout_forces_or_potential's exact three-phase
 *              pattern instead (pm_nonperiodic.c): count which task(s)
 *              each particle's mesh footprint touches, exchange
 *              positions via a plain myMPI_Alltoallv, compute locally
 *              on whichever task owns each corner, exchange results
 *              back using the same send/recv role swap.
 *
 *              Since FDM_plan is purely slab-decomposed (only X is
 *              distributed; Y and Z are fully local to every task,
 *              never split across tasks), a star's trilinear stencil
 *              touches at most 2 tasks (whichever own slab_x and
 *              slab_x+1) -- simpler than gravity's own FFT_COLUMN_BASED
 *              case (up to 4 tasks, since that decomposition splits
 *              both X and Y). We only need gravity's simpler,
 *              non-column-based branch.
 *
 *              Potential and all three force components are exchanged
 *              together in one round, not as four separate passes --
 *              a star needs all four simultaneously, and communication
 *              overhead is the same whether the payload is 1 double or
 *              4.
 *
 *              SCOPE: Type==4 (star) particles only, matching Phase
 *              2a's explicit scoping (stars first, gas/BH deferred).
 */

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"
#include "fdm.h"

#ifdef FDM

/*! \brief Per-star results, filled in by fdm_interpolate_to_stars().
 *  Indexed by the same particle index as P[]/SphP[] -- allocated in
 *  allocate.c, mirroring DMSP[]'s exact pattern (grows in lockstep
 *  with P[] itself via reallocate_memory_maxpart(), no independent
 *  reallocate_memory_* function needed).
 */
fdm_star_result *FDM_StarResult;

/*! \brief Stellar mass density deposited onto the mesh by
 *  fdm_deposit_star_mass() -- same real-space distributed layout as
 *  FDM_Potential. UNLIKE FDM_StarResult, this IS mesh data (not per-
 *  particle), so it follows this module's usual "allocate once, early,
 *  before particle arrays exist" pattern (fdm_particle_coupling_allocate(),
 *  called from fdm_allocate()) rather than DMSP[]'s per-particle
 *  pattern. Zeroed and refilled by fdm_deposit_star_mass() every time
 *  it's called -- an accumulator, not a running total across calls.
 *  NOT static -- fdm_poisson.c's fdm_update_potential() reads this
 *  directly to add the stellar source term, a genuine cross-file
 *  dependency (unlike the test-only debug accessors below, which
 *  exist specifically because FDM_StarMassDensity was originally
 *  private and only needed exposing for validation). */
fft_real *FDM_StarMassDensity;

struct fdm_star_partbuf
{
  MyDouble Pos[3];
};

struct fdm_star_massbuf
{
  MyDouble Pos[3];
  MyDouble Mass;
};

static struct fdm_star_partbuf *partout, *partin;
static fdm_star_result *flistout, *flistin;
static size_t *Sndpm_count, *Sndpm_offset, *Rcvpm_count, *Rcvpm_offset;
static size_t nexport, nimport;

/*! \brief Allocates FDM_StarMassDensity. Called from fdm_allocate() --
 *  before particle arrays exist, same movable-block-safety reasoning
 *  as every other mesh array in this module (see fdm.h's header
 *  comment on FDM_RhoLocal etc). */
void fdm_particle_coupling_allocate(void)
{
  int N = All.FDMGrid;
  size_t local_size = (size_t)FDM_plan.nslab_x * N * N;
  FDM_StarMassDensity = (fft_real *)mymalloc("FDM_StarMassDensity", local_size * sizeof(fft_real));
  memset(FDM_StarMassDensity, 0, local_size * sizeof(fft_real));
}

void fdm_particle_coupling_free(void)
{
  myfree(FDM_StarMassDensity);
}

/*! \brief Trilinear interpolation of Potential/ForceX/Y/Z from this
 *  task's own local mesh arrays, at the given LOCAL slab coordinates
 *  (slab_x already converted to this task's local indexing) and
 *  fractional offsets. Mirrors gravity's own corner-weighted formula
 *  exactly (pm_nonperiodic.c's readout function).
 */
static void fdm_trilinear_accumulate(int slab_x_local, int slab_y, int slab_z, double wx, double dy, double dz, int N,
                                      fdm_star_result *out)
{
  int yy = slab_y + 1, zz = slab_z + 1;
  /* Boundary note: a star within one cell of the box edge in Y or Z
   * would need yy/zz clamped or a neighbor fallback -- deferred here
   * since Phase 2a's own science target keeps stars well within the
   * box, matching this project's existing "central galaxy" scoping;
   * worth revisiting if stars are ever allowed to roam to the edge. */

  size_t base = (size_t)slab_x_local * N * N;

  double w000 = wx * (1.0 - dy) * (1.0 - dz);
  double w001 = wx * (1.0 - dy) * (dz);
  double w010 = wx * (dy) * (1.0 - dz);
  double w011 = wx * (dy) * (dz);

  size_t i000 = base + (size_t)slab_y * N + slab_z;
  size_t i001 = base + (size_t)slab_y * N + zz;
  size_t i010 = base + (size_t)yy * N + slab_z;
  size_t i011 = base + (size_t)yy * N + zz;

  out->Potential += FDM_Potential[i000] * w000 + FDM_Potential[i001] * w001 + FDM_Potential[i010] * w010 + FDM_Potential[i011] * w011;
  out->ForceX += FDM_ForceX[i000] * w000 + FDM_ForceX[i001] * w001 + FDM_ForceX[i010] * w010 + FDM_ForceX[i011] * w011;
  out->ForceY += FDM_ForceY[i000] * w000 + FDM_ForceY[i001] * w001 + FDM_ForceY[i010] * w010 + FDM_ForceY[i011] * w011;
  out->ForceZ += FDM_ForceZ[i000] * w000 + FDM_ForceZ[i001] * w001 + FDM_ForceZ[i010] * w010 + FDM_ForceZ[i011] * w011;
}

/*! \brief Is this particle in scope for FDM coupling? Phase 2a: stars
 *  only (Type==4, Arepo's standard convention). */
static inline int fdm_particle_in_scope(int i) { return P[i].Type == 4; }

void fdm_interpolate_to_stars(void)
{
  int N = All.FDMGrid;
  double dx_cell = All.FDMBoxSize / N;
  double to_slab_fac = 1.0 / dx_cell;

  int multiNtask = NTask; /* single-threaded version -- no per-thread cache-line padding needed here, unlike gravity's OpenMP-parallel prepare_density */

  Sndpm_count  = (size_t *)mymalloc("Sndpm_count", multiNtask * sizeof(size_t));
  Sndpm_offset = (size_t *)mymalloc("Sndpm_offset", multiNtask * sizeof(size_t));
  Rcvpm_count  = (size_t *)mymalloc("Rcvpm_count", NTask * sizeof(size_t));
  Rcvpm_offset = (size_t *)mymalloc("Rcvpm_offset", NTask * sizeof(size_t));

  for(int j = 0; j < NTask; j++)
    Sndpm_count[j] = 0;

  /* Pass 1: count. */
  for(int i = 0; i < NumPart; i++)
    {
      if(!fdm_particle_in_scope(i))
        continue;

      double px = P[i].Pos[0], py = P[i].Pos[1], pz = P[i].Pos[2];
      if(px < 0 || px >= All.FDMBoxSize || py < 0 || py >= All.FDMBoxSize || pz < 0 || pz >= All.FDMBoxSize)
        continue; /* outside the FDM box entirely -- no force contribution from this module */

      int slab_x  = (int)(to_slab_fac * px);
      int slab_xx = slab_x + 1;
      if(slab_xx >= N)
        slab_xx = slab_x; /* clamp at the box edge rather than wrap -- isolated domain, matches fdm_gradient.c's own boundary treatment */

      int task0 = FDM_plan.slab_to_task[slab_x];
      int task1 = FDM_plan.slab_to_task[slab_xx];

      Sndpm_count[task0]++;
      if(task1 != task0)
        Sndpm_count[task1]++;
    }

  MPI_Alltoall(Sndpm_count, sizeof(size_t), MPI_BYTE, Rcvpm_count, sizeof(size_t), MPI_BYTE, MPI_COMM_WORLD);

  nexport = 0;
  nimport = 0;
  Sndpm_offset[0] = 0;
  Rcvpm_offset[0] = 0;
  for(int j = 0; j < NTask; j++)
    {
      nexport += Sndpm_count[j];
      nimport += Rcvpm_count[j];
      if(j > 0)
        {
          Sndpm_offset[j] = Sndpm_offset[j - 1] + Sndpm_count[j - 1];
          Rcvpm_offset[j] = Rcvpm_offset[j - 1] + Rcvpm_count[j - 1];
        }
    }

  partin  = (struct fdm_star_partbuf *)mymalloc("fdm_partin", nimport * sizeof(struct fdm_star_partbuf));
  partout = (struct fdm_star_partbuf *)mymalloc("fdm_partout", nexport * sizeof(struct fdm_star_partbuf));

  /* Pass 2: populate export buffer -- SAME iteration order and SAME
   * task0/task1 determination as Pass 1, relied upon again in Pass 4
   * to correctly match returned results back to the right particle
   * (mirrors gravity's own pattern exactly). */
  {
    size_t *send_count = (size_t *)mymalloc("fdm_send_count_tmp", NTask * sizeof(size_t));
    for(int j = 0; j < NTask; j++)
      send_count[j] = 0;

    for(int i = 0; i < NumPart; i++)
      {
        if(!fdm_particle_in_scope(i))
          continue;

        double px = P[i].Pos[0], py = P[i].Pos[1], pz = P[i].Pos[2];
        if(px < 0 || px >= All.FDMBoxSize || py < 0 || py >= All.FDMBoxSize || pz < 0 || pz >= All.FDMBoxSize)
          continue;

        int slab_x  = (int)(to_slab_fac * px);
        int slab_xx = slab_x + 1;
        if(slab_xx >= N)
          slab_xx = slab_x;

        int task0 = FDM_plan.slab_to_task[slab_x];
        int task1 = FDM_plan.slab_to_task[slab_xx];

        size_t ind0 = Sndpm_offset[task0] + send_count[task0]++;
        partout[ind0].Pos[0] = px;
        partout[ind0].Pos[1] = py;
        partout[ind0].Pos[2] = pz;

        if(task1 != task0)
          {
            size_t ind1 = Sndpm_offset[task1] + send_count[task1]++;
            partout[ind1].Pos[0] = px;
            partout[ind1].Pos[1] = py;
            partout[ind1].Pos[2] = pz;
          }
      }

    myfree(send_count);
  }

  int flag_big = 0, flag_big_all;
  for(int j = 0; j < NTask; j++)
    if(Sndpm_count[j] * sizeof(struct fdm_star_partbuf) > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
      flag_big = 1;
  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  myMPI_Alltoallv(partout, Sndpm_count, Sndpm_offset, partin, Rcvpm_count, Rcvpm_offset, sizeof(struct fdm_star_partbuf), flag_big_all,
                  MPI_COMM_WORLD);

  myfree(partout);

  /* Pass 3: compute locally, for each imported position, whichever
   * corners THIS task owns. Since our decomposition only splits X,
   * a given imported position needs data from AT MOST one X-slab on
   * THIS task (the other slab, if any, was independently exported to
   * whichever task owns it -- gravity's own pattern, not something we
   * check for here). */
  flistin  = (fdm_star_result *)mymalloc("fdm_flistin", nimport * sizeof(fdm_star_result));
  flistout = (fdm_star_result *)mymalloc("fdm_flistout", nexport * sizeof(fdm_star_result));

  for(size_t i = 0; i < nimport; i++)
    {
      flistin[i].Potential = flistin[i].ForceX = flistin[i].ForceY = flistin[i].ForceZ = 0.0;

      double px = partin[i].Pos[0], py = partin[i].Pos[1], pz = partin[i].Pos[2];

      int slab_x = (int)(to_slab_fac * px);
      int slab_y = (int)(to_slab_fac * py);
      int slab_z = (int)(to_slab_fac * pz);

      double dxf = to_slab_fac * px - slab_x;
      double dyf = to_slab_fac * py - slab_y;
      double dzf = to_slab_fac * pz - slab_z;

      int slab_xx = slab_x + 1;
      if(slab_xx >= N)
        slab_xx = slab_x;

      if(FDM_plan.slab_to_task[slab_x] == ThisTask)
        {
          int local_x = slab_x - FDM_plan.slabstart_x;
          fdm_trilinear_accumulate(local_x, slab_y, slab_z, (1.0 - dxf), dyf, dzf, N, &flistin[i]);
        }

      if(FDM_plan.slab_to_task[slab_xx] == ThisTask && slab_xx != slab_x)
        {
          int local_xx = slab_xx - FDM_plan.slabstart_x;
          fdm_trilinear_accumulate(local_xx, slab_y, slab_z, dxf, dyf, dzf, N, &flistin[i]);
        }
    }

  myMPI_Alltoallv(flistin, Rcvpm_count, Rcvpm_offset, flistout, Sndpm_count, Sndpm_offset, sizeof(fdm_star_result), flag_big_all,
                  MPI_COMM_WORLD);

  /* flistin/partin are logically done with here, but CANNOT be freed
   * yet -- flistout (allocated after both, used below in Pass 4) is
   * still alive and needed. Deferred to the end of this function, in
   * exact reverse allocation order. Same mistake already made once
   * before in this exact spot with partout/partin (freeing the older
   * of two allocations while the newer was still alive) -- caught
   * immediately by mymalloc's own safety check, not silent
   * corruption, but a repeat of a lesson from earlier work on this
   * module that should have been applied here from the start. */

  /* Pass 4: reassemble -- same iteration order as Passes 1/2 again. */
  {
    size_t *send_count = (size_t *)mymalloc("fdm_send_count_tmp2", NTask * sizeof(size_t));
    for(int j = 0; j < NTask; j++)
      send_count[j] = 0;

    for(int i = 0; i < NumPart; i++)
      {
        FDM_StarResult[i].Potential = FDM_StarResult[i].ForceX = FDM_StarResult[i].ForceY = FDM_StarResult[i].ForceZ = 0.0;

        if(!fdm_particle_in_scope(i))
          continue;

        double px = P[i].Pos[0], py = P[i].Pos[1], pz = P[i].Pos[2];
        if(px < 0 || px >= All.FDMBoxSize || py < 0 || py >= All.FDMBoxSize || pz < 0 || pz >= All.FDMBoxSize)
          continue;

        int slab_x  = (int)(to_slab_fac * px);
        int slab_xx = slab_x + 1;
        if(slab_xx >= N)
          slab_xx = slab_x;

        int task0 = FDM_plan.slab_to_task[slab_x];
        int task1 = FDM_plan.slab_to_task[slab_xx];

        fdm_star_result r = flistout[Sndpm_offset[task0] + send_count[task0]++];
        FDM_StarResult[i].Potential += r.Potential;
        FDM_StarResult[i].ForceX += r.ForceX;
        FDM_StarResult[i].ForceY += r.ForceY;
        FDM_StarResult[i].ForceZ += r.ForceZ;

        if(task1 != task0)
          {
            fdm_star_result r2 = flistout[Sndpm_offset[task1] + send_count[task1]++];
            FDM_StarResult[i].Potential += r2.Potential;
            FDM_StarResult[i].ForceX += r2.ForceX;
            FDM_StarResult[i].ForceY += r2.ForceY;
            FDM_StarResult[i].ForceZ += r2.ForceZ;
          }
      }

    myfree(send_count);
  }

  myfree(flistout);
  myfree(flistin);
  myfree(partin);
  myfree(Rcvpm_offset);
  myfree(Rcvpm_count);
  myfree(Sndpm_offset);
  myfree(Sndpm_count);
}

/*! \brief Deposits star particle masses onto the FDM mesh via CIC,
 *  into FDM_StarMassDensity -- the star->mesh half of the coupling
 *  (fdm_interpolate_to_stars() above is the mesh->star half). One-way:
 *  unlike interpolation, deposition needs no "send results back" step
 *  -- the mesh-owning task simply accumulates the weighted mass
 *  contribution directly into its own local array, mirroring gravity's
 *  own pmforce_nonperiodic_uniform_optimized_prepare_density exactly
 *  (that function's rhogrid[...] += (mass*weight) pattern, verified
 *  directly rather than assumed, given this project's history of
 *  wrong assumptions about "obvious" send/receive symmetry in this
 *  general area).
 *
 *  Does NOT itself feed the result into fdm_update_potential()'s
 *  source term -- that wiring, and calling this function at the right
 *  point in the timestep, are separate, later integration steps.
 */
void fdm_deposit_star_mass(void)
{
  int N = All.FDMGrid;
  double dx_cell = All.FDMBoxSize / N;
  double to_slab_fac = 1.0 / dx_cell;

  memset(FDM_StarMassDensity, 0, (size_t)FDM_plan.nslab_x * N * N * sizeof(fft_real));

  Sndpm_count  = (size_t *)mymalloc("Sndpm_count", NTask * sizeof(size_t));
  Sndpm_offset = (size_t *)mymalloc("Sndpm_offset", NTask * sizeof(size_t));
  Rcvpm_count  = (size_t *)mymalloc("Rcvpm_count", NTask * sizeof(size_t));
  Rcvpm_offset = (size_t *)mymalloc("Rcvpm_offset", NTask * sizeof(size_t));

  for(int j = 0; j < NTask; j++)
    Sndpm_count[j] = 0;

  /* Pass 1: count -- identical logic to fdm_interpolate_to_stars()'s
   * own Pass 1 (same task0/task1 determination from position alone;
   * mass doesn't affect which tasks are touched). */
  for(int i = 0; i < NumPart; i++)
    {
      if(!fdm_particle_in_scope(i))
        continue;

      double px = P[i].Pos[0], py = P[i].Pos[1], pz = P[i].Pos[2];
      if(px < 0 || px >= All.FDMBoxSize || py < 0 || py >= All.FDMBoxSize || pz < 0 || pz >= All.FDMBoxSize)
        continue;

      int slab_x  = (int)(to_slab_fac * px);
      int slab_xx = slab_x + 1;
      if(slab_xx >= N)
        slab_xx = slab_x;

      int task0 = FDM_plan.slab_to_task[slab_x];
      int task1 = FDM_plan.slab_to_task[slab_xx];

      Sndpm_count[task0]++;
      if(task1 != task0)
        Sndpm_count[task1]++;
    }

  MPI_Alltoall(Sndpm_count, sizeof(size_t), MPI_BYTE, Rcvpm_count, sizeof(size_t), MPI_BYTE, MPI_COMM_WORLD);

  nexport = 0;
  nimport = 0;
  Sndpm_offset[0] = 0;
  Rcvpm_offset[0] = 0;
  for(int j = 0; j < NTask; j++)
    {
      nexport += Sndpm_count[j];
      nimport += Rcvpm_count[j];
      if(j > 0)
        {
          Sndpm_offset[j] = Sndpm_offset[j - 1] + Sndpm_count[j - 1];
          Rcvpm_offset[j] = Rcvpm_offset[j - 1] + Rcvpm_count[j - 1];
        }
    }

  /* massin before massout -- same LIFO lesson as
   * fdm_interpolate_to_stars() (massout must be freed immediately
   * after the exchange, so it must be the LATER allocation). */
  struct fdm_star_massbuf *massin, *massout;
  massin  = (struct fdm_star_massbuf *)mymalloc("fdm_massin", nimport * sizeof(struct fdm_star_massbuf));
  massout = (struct fdm_star_massbuf *)mymalloc("fdm_massout", nexport * sizeof(struct fdm_star_massbuf));

  /* Pass 2: populate. */
  {
    size_t *send_count = (size_t *)mymalloc("fdm_send_count_tmp", NTask * sizeof(size_t));
    for(int j = 0; j < NTask; j++)
      send_count[j] = 0;

    for(int i = 0; i < NumPart; i++)
      {
        if(!fdm_particle_in_scope(i))
          continue;

        double px = P[i].Pos[0], py = P[i].Pos[1], pz = P[i].Pos[2];
        if(px < 0 || px >= All.FDMBoxSize || py < 0 || py >= All.FDMBoxSize || pz < 0 || pz >= All.FDMBoxSize)
          continue;

        int slab_x  = (int)(to_slab_fac * px);
        int slab_xx = slab_x + 1;
        if(slab_xx >= N)
          slab_xx = slab_x;

        int task0 = FDM_plan.slab_to_task[slab_x];
        int task1 = FDM_plan.slab_to_task[slab_xx];

        size_t ind0 = Sndpm_offset[task0] + send_count[task0]++;
        massout[ind0].Pos[0] = px;
        massout[ind0].Pos[1] = py;
        massout[ind0].Pos[2] = pz;
        massout[ind0].Mass   = P[i].Mass;

        if(task1 != task0)
          {
            size_t ind1 = Sndpm_offset[task1] + send_count[task1]++;
            massout[ind1].Pos[0] = px;
            massout[ind1].Pos[1] = py;
            massout[ind1].Pos[2] = pz;
            massout[ind1].Mass   = P[i].Mass;
          }
      }

    myfree(send_count);
  }

  int flag_big = 0, flag_big_all;
  for(int j = 0; j < NTask; j++)
    if(Sndpm_count[j] * sizeof(struct fdm_star_massbuf) > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
      flag_big = 1;
  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  myMPI_Alltoallv(massout, Sndpm_count, Sndpm_offset, massin, Rcvpm_count, Rcvpm_offset, sizeof(struct fdm_star_massbuf), flag_big_all,
                  MPI_COMM_WORLD);

  myfree(massout);

  /* Pass 3: deposit -- ONE-WAY, no results sent back. Accumulate
   * (+=, not =) directly into FDM_StarMassDensity, since multiple
   * stars can contribute to the same cell. */
  for(size_t i = 0; i < nimport; i++)
    {
      double px = massin[i].Pos[0], py = massin[i].Pos[1], pz = massin[i].Pos[2];
      double mass = massin[i].Mass;

      int slab_x = (int)(to_slab_fac * px);
      int slab_y = (int)(to_slab_fac * py);
      int slab_z = (int)(to_slab_fac * pz);

      double dxf = to_slab_fac * px - slab_x;
      double dyf = to_slab_fac * py - slab_y;
      double dzf = to_slab_fac * pz - slab_z;

      int slab_xx = slab_x + 1;
      if(slab_xx >= N)
        slab_xx = slab_x;

      int yy = slab_y + 1 < N ? slab_y + 1 : slab_y;
      int zz = slab_z + 1 < N ? slab_z + 1 : slab_z;
      /* Same boundary note as fdm_trilinear_accumulate(): clamped, not
       * wrapped, since stars are expected well within the box for
       * Phase 2a's science target -- worth revisiting if that changes. */

      if(FDM_plan.slab_to_task[slab_x] == ThisTask)
        {
          int lx = slab_x - FDM_plan.slabstart_x;
          size_t base = (size_t)lx * N * N;
          FDM_StarMassDensity[base + (size_t)slab_y * N + slab_z] += mass * (1.0 - dxf) * (1.0 - dyf) * (1.0 - dzf);
          FDM_StarMassDensity[base + (size_t)slab_y * N + zz] += mass * (1.0 - dxf) * (1.0 - dyf) * (dzf);
          FDM_StarMassDensity[base + (size_t)yy * N + slab_z] += mass * (1.0 - dxf) * (dyf) * (1.0 - dzf);
          FDM_StarMassDensity[base + (size_t)yy * N + zz] += mass * (1.0 - dxf) * (dyf) * (dzf);
        }

      if(FDM_plan.slab_to_task[slab_xx] == ThisTask && slab_xx != slab_x)
        {
          int lxx = slab_xx - FDM_plan.slabstart_x;
          size_t base = (size_t)lxx * N * N;
          FDM_StarMassDensity[base + (size_t)slab_y * N + slab_z] += mass * (dxf) * (1.0 - dyf) * (1.0 - dzf);
          FDM_StarMassDensity[base + (size_t)slab_y * N + zz] += mass * (dxf) * (1.0 - dyf) * (dzf);
          FDM_StarMassDensity[base + (size_t)yy * N + slab_z] += mass * (dxf) * (dyf) * (1.0 - dzf);
          FDM_StarMassDensity[base + (size_t)yy * N + zz] += mass * (dxf) * (dyf) * (dzf);
        }
    }

  myfree(massin);
  myfree(Rcvpm_offset);
  myfree(Rcvpm_count);
  myfree(Sndpm_offset);
  myfree(Sndpm_count);
}

/*! TEST-ONLY helpers exposing FDM_StarMassDensity (intentionally
 * static/private otherwise) for validation. */
double fdm_debug_sum_star_mass_density(void)
{
  int N = All.FDMGrid;
  size_t local_size = (size_t)FDM_plan.nslab_x * N * N;
  double sum = 0.0;
  for(size_t idx = 0; idx < local_size; idx++)
    sum += FDM_StarMassDensity[idx];
  return sum;
}

double fdm_debug_get_star_mass_density_at(int gx, int gy, int gz)
{
  int N = All.FDMGrid;
  if(gx < FDM_plan.slabstart_x || gx >= FDM_plan.slabstart_x + FDM_plan.nslab_x)
    return 0.0; /* not this task's slab -- caller sums across tasks via Allreduce */
  int lx = gx - FDM_plan.slabstart_x;
  return FDM_StarMassDensity[(size_t)lx * N * N + (size_t)gy * N + gz];
}

#endif /* #ifdef FDM */
