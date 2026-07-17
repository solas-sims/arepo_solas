/*!
 * \file        src/fdm/fdm_gradient.c
 * \brief       Computes the force field (-grad(Phi)) on the mesh from
 *              FDM_Potential, via 4th-order-accurate finite
 *              differencing -- the first piece of Phase 2a's particle
 *              coupling (mesh->star force interpolation is a separate,
 *              later piece).
 *
 *              STENCIL: mirrors pm_nonperiodic.c's own force-from-
 *              potential finite difference exactly in FORM:
 *                (4/3)*(Phi[i-1]-Phi[i+1]) - (1/6)*(Phi[i-2]-Phi[i+2])
 *              but NOT in normalization -- gravity's own "fac" constant
 *              is entangled with its k-space/erfc-split/TotalMeshSize
 *              convention, unrelated to our own (real-space, zero-
 *              padded) Poisson solve's normalization. Verified
 *              independently via symbolic Taylor expansion that this
 *              stencil equals -2*dx*d/dx(Phi) + O(dx^5); the force is
 *              therefore the stencil divided by 2*dx directly (no
 *              extra sign flip -- an earlier version double-negated
 *              this, exactly cancelling the correct sign, caught by a
 *              quadratic-potential test showing relative error ~2.0,
 *              the signature of a full sign flip).
 *
 *              TRANSPOSE LAYOUT (the hard-won part of this file):
 *              gravity's own TI(x,y,z) macro (pm_nonperiodic.c) is the
 *              authoritative specification of the post-transpose array
 *              layout: TI(x,y,z) = GRID*(x + y*nslab_x) + z. Decomposing
 *              strides -- z fastest (1), x (the LOCAL, nslab_x-ranged
 *              loop variable) has the MIDDLE stride (GRID), and y (the
 *              FULL-RANGE loop variable) has the LARGEST/OUTERMOST
 *              stride (GRID*nslab_x). An earlier version of this file
 *              had these backwards (local dimension outermost, full-
 *              range dimension in the middle) -- confirmed wrong by a
 *              quadratic-potential test (whose derivative genuinely
 *              varies with position, unlike an earlier linear test
 *              that gave a false-positive pass, since a linear
 *              function's slope doesn't depend on where you sample it).
 *              The fix: the full-range (X) dimension must be
 *              OUTERMOST (largest stride), the local Y-slab dimension
 *              is the MIDDLE dimension -- matching gravity's own
 *              layout exactly, which makes sense given
 *              fdm_transpose_x_y_A/B are faithful copies of gravity's
 *              own proven pack/unpack algorithm (just with stride N
 *              instead of Ngrid2, see below) -- the transpose functions
 *              themselves were never the bug; only this file's own
 *              indexing into their output was.
 *
 *              Y and Z differencing is purely local (a task has all
 *              Y/Z data for its own X-slabs) -- no transpose needed for
 *              those two components.
 *
 *              TRANSPOSE FUNCTIONS: gravity's own my_slab_transposeA/B
 *              (pm_mpi_fft.c) compute offsets using plan->Ngrid2 (the
 *              r2c-padding-specific stride) -- NOT plan->NgridX, since
 *              gravity's own rhogrid/forcegrid use that padded layout.
 *              Our FDM_Potential deliberately does NOT use that padding
 *              (genuinely complex/unpadded throughout this module, see
 *              fdm.h) -- reusing those functions directly would
 *              silently compute wrong offsets against our differently-
 *              strided array. fdm_transpose_x_y_A/B below mirror their
 *              exact algorithm and MPI communication pattern, with the
 *              stride fixed to N (this module's own layout) instead of
 *              Ngrid2.
 *
 *              BOUNDARY: cells within 2 grid points of the mesh edge
 *              fall back to a lower-order (2-point, then 1-point at the
 *              absolute edge) centered/one-sided difference -- wrapping
 *              the 4-point stencil around the array boundary would
 *              incorrectly mix data from opposite sides of an isolated
 *              (non-periodic) domain as if they were physically
 *              adjacent.
 *
 *              SCRATCH BUFFERS: dedicated, never reused as temporary
 *              space for a DIFFERENT array's own output -- pre-
 *              allocated once, early (fdm_gradient_allocate, called
 *              from fdm_allocate() before P[]/SphP[] exist), matching
 *              this module's established movable-block-safety pattern
 *              (see fdm.h's header comment on FDM_RhoLocal etc for the
 *              full reasoning -- a genuine crash, not a style choice).
 */

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"
#include "fdm.h"

#ifdef FDM

fft_real *FDM_ForceX;
fft_real *FDM_ForceY;
fft_real *FDM_ForceZ;

static fft_real *FDM_CommScratch;         /* used internally by the transpose functions themselves */
static fft_real *FDM_TransposedPotential; /* Phi, transposed so X is local */
static fft_real *FDM_TransposedForceX;    /* -dPhi/dx, in the transposed layout, before transposing back */

/*! \brief Transposes `field` so that Y becomes the locally-distributed
 *  dimension and X becomes fully local -- mirrors my_slab_transposeA's
 *  exact algorithm and communication pattern (pm_mpi_fft.c), but with
 *  stride N (this module's own unpadded layout) instead of Ngrid2.
 */
static void fdm_transpose_x_y_A(fft_plan *plan, fft_real *field, fft_real *scratch)
{
  int N = All.FDMGrid;
  int prod = NTask * plan->nslab_x;

  for(int n = 0; n < prod; n++)
    {
      int x    = n / NTask;
      int task = n % NTask;

      for(int y = plan->first_slab_y_of_task[task]; y < plan->first_slab_y_of_task[task] + plan->slabs_y_per_task[task]; y++)
        memcpy(scratch + (size_t)N * (plan->first_slab_y_of_task[task] * plan->nslab_x + x * plan->slabs_y_per_task[task] +
                                       (y - plan->first_slab_y_of_task[task])),
               field + (size_t)N * (N * x + y), N * sizeof(fft_real));
    }

  size_t *scount = (size_t *)mymalloc("fdm_transpose_scount", NTask * sizeof(size_t));
  size_t *rcount = (size_t *)mymalloc("fdm_transpose_rcount", NTask * sizeof(size_t));
  size_t *soff   = (size_t *)mymalloc("fdm_transpose_soff", NTask * sizeof(size_t));
  size_t *roff   = (size_t *)mymalloc("fdm_transpose_roff", NTask * sizeof(size_t));

  int flag_big = 0, flag_big_all = 0;
  for(int task = 0; task < NTask; task++)
    {
      scount[task] = (size_t)plan->nslab_x * plan->slabs_y_per_task[task] * (N * sizeof(fft_real));
      rcount[task] = (size_t)plan->nslab_y * plan->slabs_x_per_task[task] * (N * sizeof(fft_real));

      soff[task] = (size_t)plan->first_slab_y_of_task[task] * plan->nslab_x * (N * sizeof(fft_real));
      roff[task] = (size_t)plan->first_slab_x_of_task[task] * plan->nslab_y * (N * sizeof(fft_real));

      if(scount[task] > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
        flag_big = 1;
    }

  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  myMPI_Alltoallv(scratch, scount, soff, field, rcount, roff, 1, flag_big_all, MPI_COMM_WORLD);

  myfree(roff);
  myfree(soff);
  myfree(rcount);
  myfree(scount);
}

/*! \brief Undoes fdm_transpose_x_y_A -- mirrors my_slab_transposeB
 *  exactly (stride N instead of Ngrid2, same reasoning as above). */
static void fdm_transpose_x_y_B(fft_plan *plan, fft_real *field, fft_real *scratch)
{
  int N = All.FDMGrid;

  size_t *scount = (size_t *)mymalloc("fdm_transpose_scount", NTask * sizeof(size_t));
  size_t *rcount = (size_t *)mymalloc("fdm_transpose_rcount", NTask * sizeof(size_t));
  size_t *soff   = (size_t *)mymalloc("fdm_transpose_soff", NTask * sizeof(size_t));
  size_t *roff   = (size_t *)mymalloc("fdm_transpose_roff", NTask * sizeof(size_t));

  int flag_big = 0, flag_big_all = 0;
  for(int task = 0; task < NTask; task++)
    {
      rcount[task] = (size_t)plan->nslab_x * plan->slabs_y_per_task[task] * (N * sizeof(fft_real));
      scount[task] = (size_t)plan->nslab_y * plan->slabs_x_per_task[task] * (N * sizeof(fft_real));

      roff[task] = (size_t)plan->first_slab_y_of_task[task] * plan->nslab_x * (N * sizeof(fft_real));
      soff[task] = (size_t)plan->first_slab_x_of_task[task] * plan->nslab_y * (N * sizeof(fft_real));

      if(scount[task] > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
        flag_big = 1;
    }

  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  myMPI_Alltoallv(field, scount, soff, scratch, rcount, roff, 1, flag_big_all, MPI_COMM_WORLD);

  myfree(roff);
  myfree(soff);
  myfree(rcount);
  myfree(scount);

  int prod = NTask * plan->nslab_x;
  for(int n = 0; n < prod; n++)
    {
      int x    = n / NTask;
      int task = n % NTask;

      for(int y = plan->first_slab_y_of_task[task]; y < plan->first_slab_y_of_task[task] + plan->slabs_y_per_task[task]; y++)
        memcpy(field + (size_t)N * (N * x + y),
               scratch + (size_t)N * (plan->first_slab_y_of_task[task] * plan->nslab_x + x * plan->slabs_y_per_task[task] +
                                       (y - plan->first_slab_y_of_task[task])),
               N * sizeof(fft_real));
    }
}

/*! \brief The verified 4th-order stencil (falling back to lower order
 *  near the boundary), applied to a column with the given stride --
 *  handles Z (stride 1), Y (stride N), and (after transposing) X all
 *  with one implementation. */
static inline double fdm_diff1d_strided(fft_real *field, size_t base, int idx, int N, size_t stride, double dx)
{
  if(idx >= 2 && idx <= N - 3)
    return ((4.0 / 3.0) * (field[base + (size_t)(idx - 1) * stride] - field[base + (size_t)(idx + 1) * stride]) -
            (1.0 / 6.0) * (field[base + (size_t)(idx - 2) * stride] - field[base + (size_t)(idx + 2) * stride])) /
           (2.0 * dx);
  else if(idx >= 1 && idx <= N - 2)
    return (field[base + (size_t)(idx - 1) * stride] - field[base + (size_t)(idx + 1) * stride]) / (2.0 * dx);
  else if(idx == 0)
    return (field[base] - field[base + stride]) / dx;
  else /* idx == N-1 */
    return (field[base + (size_t)(idx - 1) * stride] - field[base + (size_t)idx * stride]) / dx;
}

/*! \brief Allocates FDM_ForceX/Y/Z and this file's dedicated scratch
 *  buffers. Called from fdm_allocate() -- same "allocate once, early,
 *  before particle arrays exist" requirement as everything else in
 *  this module.
 *
 *  FDM_TransposedPotential/FDM_TransposedForceX are sized using
 *  largest_y_slab (the GLOBAL maximum y-slab size across all tasks)
 *  rather than this task's own nslab_y -- every task must size its own
 *  buffer large enough for whatever it might receive, not just its own
 *  share.
 */
void fdm_gradient_allocate(void)
{
  int N = All.FDMGrid;
  size_t local_size_x = (size_t)FDM_plan.nslab_x * N * N;
  size_t local_size_y_max = (size_t)FDM_plan.largest_y_slab * N * N;
  size_t comm_scratch_size = (size_t)(FDM_plan.largest_x_slab > FDM_plan.largest_y_slab ? FDM_plan.largest_x_slab
                                                                                          : FDM_plan.largest_y_slab) *
                             N * N;

  FDM_ForceX = (fft_real *)mymalloc("FDM_ForceX", local_size_x * sizeof(fft_real));
  FDM_ForceY = (fft_real *)mymalloc("FDM_ForceY", local_size_x * sizeof(fft_real));
  FDM_ForceZ = (fft_real *)mymalloc("FDM_ForceZ", local_size_x * sizeof(fft_real));

  FDM_TransposedPotential = (fft_real *)mymalloc("FDM_TransposedPotential", local_size_y_max * sizeof(fft_real));
  FDM_TransposedForceX    = (fft_real *)mymalloc("FDM_TransposedForceX", local_size_y_max * sizeof(fft_real));
  FDM_CommScratch         = (fft_real *)mymalloc("FDM_CommScratch", comm_scratch_size * sizeof(fft_real));
}

/*! \brief Frees this file's allocations, in exact reverse of
 *  fdm_gradient_allocate()'s order. Called from fdm_free(). */
void fdm_gradient_free(void)
{
  myfree(FDM_CommScratch);
  myfree(FDM_TransposedForceX);
  myfree(FDM_TransposedPotential);
  myfree(FDM_ForceZ);
  myfree(FDM_ForceY);
  myfree(FDM_ForceX);
}

/*! \brief Computes FDM_ForceX/Y/Z (= -grad(FDM_Potential)) via finite
 *  differencing. Call after fdm_update_potential() -- consumes
 *  whatever is currently in FDM_Potential.
 */
void fdm_compute_force(void)
{
  int N = All.FDMGrid;
  double dx = All.FDMBoxSize / N;

  /* Z: stride 1, always local. */
  for(int i = 0; i < FDM_plan.nslab_x; i++)
    for(int j = 0; j < N; j++)
      {
        size_t base = (size_t)i * N * N + (size_t)j * N;
        for(int k = 0; k < N; k++)
          FDM_ForceZ[base + k] = fdm_diff1d_strided(FDM_Potential, base, k, N, 1, dx);
      }

  /* Y: stride N, always local. */
  for(int i = 0; i < FDM_plan.nslab_x; i++)
    for(int k = 0; k < N; k++)
      {
        size_t base = (size_t)i * N * N + k;
        for(int j = 0; j < N; j++)
          FDM_ForceY[base + (size_t)j * N] = fdm_diff1d_strided(FDM_Potential, base, j, N, N, dx);
      }

  /* X: needs neighbouring slabs, potentially on a different task --
   * transpose Phi so X is local, difference it, transpose the result
   * back.
   *
   * POST-TRANSPOSE LAYOUT (see file header for the full derivation
   * against gravity's own TI(x,y,z) macro): index = q*N*nslab_y +
   * p*N + k, where q is the GLOBAL X value (0..N-1, OUTERMOST/largest
   * stride N*nslab_y) and p is the LOCAL Y-slab index (0..nslab_y-1,
   * MIDDLE stride N). This is the opposite of the naive "local
   * dimension outermost" assumption an earlier version of this file
   * used -- confirmed via gravity's own macro, not just re-tested
   * empirically, given this project's own history of misleading
   * empirical tests in exactly this spot. */
  size_t local_size = (size_t)FDM_plan.nslab_x * N * N;
  memcpy(FDM_TransposedPotential, FDM_Potential, local_size * sizeof(fft_real));
  fdm_transpose_x_y_A(&FDM_plan, FDM_TransposedPotential, FDM_CommScratch);

  int nslab_y = FDM_plan.nslab_y;
  for(int p = 0; p < nslab_y; p++)
    for(int k = 0; k < N; k++)
      {
        size_t base = (size_t)p * N + k; /* p is the MIDDLE dimension here, stride N */
        for(int q = 0; q < N; q++)
          FDM_TransposedForceX[(size_t)q * N * nslab_y + base] =
              fdm_diff1d_strided(FDM_TransposedPotential, base, q, N, (size_t)N * nslab_y, dx);
      }

  memcpy(FDM_TransposedPotential, FDM_TransposedForceX, local_size * sizeof(fft_real));
  fdm_transpose_x_y_B(&FDM_plan, FDM_TransposedPotential, FDM_CommScratch);
  memcpy(FDM_ForceX, FDM_TransposedPotential, local_size * sizeof(fft_real));
}

#endif /* #ifdef FDM */
