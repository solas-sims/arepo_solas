/*!
 * \file        src/fdm/fdm_poisson.c
 * \brief       Isolated (non-periodic) Poisson solve for the FDM
 *              potential, via the standard Hockney & Eastwood
 *              zero-padding technique: a doubled mesh (2*All.FDMGrid
 *              per dimension), a plain (unsplit) -1/r Green's
 *              function kernel built once and cached, and FFT ->
 *              multiply -> inverse FFT each call.
 *
 *              REDISTRIBUTION: FDM_plan (size N) and FDM_PoissonPlan
 *              (size 2N) have independently-computed slab
 *              decompositions -- a task's local x-rows in one do NOT
 *              generally align with its local x-rows in the other,
 *              since they're different total sizes divided among the
 *              same NTask. Real MPI communication is required to move
 *              data between them; there is no way to avoid this by
 *              clever indexing alone. Both plans' slab_to_task[]
 *              arrays are fully populated on every task (via the
 *              Allgather in my_slab_based_fft_init), so every task can
 *              independently compute the full send/receive plan
 *              without a separate negotiation round.
 */

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"
#include "fdm.h"

#ifdef FDM

fft_plan FDM_PoissonPlan;
fft_complex *FDM_Kernel;
fft_complex *FDM_PoissonWorkspace;

fft_complex *FDM_RhoLocal;
fft_complex *FDM_RhoPadded;
fft_complex *FDM_PhiLocal;
fft_complex *FDM_RedistPackBuf;
void        *FDM_RedistRequests;

/*! \brief Redistributes the [0,N)x[0,N) (y,z) sub-block of each row from
 *  a source plan's decomposition to a destination plan's decomposition,
 *  where N is transfer_width -- the "real" (non-padded) extent common
 *  to both directions this is used in.
 *
 *  CRITICAL: src and dst rows are NOT necessarily the same width.
 *  FDM_plan's rows are N wide (in y and z); FDM_PoissonPlan's rows are
 *  N2=2N wide. An earlier version of this function treated a "row" as
 *  a flat N*N-element contiguous block and copied it directly at the
 *  same stride on both sides -- which is only valid when both sides
 *  have equal row width. Going from N-wide rows into N2-wide rows (or
 *  back), the N*N source block does NOT sit contiguously at the start
 *  of the destination's actual N2*N2-element row; each of the N
 *  "y-slices" needs to land at its own strided offset (j*N2, not
 *  j*N), not just the first N*N elements of the row. This was found
 *  by a diagnostic showing correct data on the send side but zero on
 *  the receive side at the same physical point -- not from reasoning
 *  about the layout directly, which is exactly why it needed a
 *  dedicated single-mode-style check rather than trusting the "obvious"
 *  flat-copy implementation.
 *
 *  Fix: pack the [0,N)x[0,N) sub-block into a contiguous scratch
 *  buffer before sending (when the source row is wider than N), and
 *  unpack a contiguous received buffer into the correct strided
 *  destination positions (when the destination row is wider than N).
 *  When a side's own row width already equals transfer_width, no
 *  packing/unpacking is needed on that side -- used directly.
 */
static void fdm_redistribute_rows(fft_complex *src_rows, fft_plan *src_plan, int src_row_width, fft_complex *dst_rows,
                                   fft_plan *dst_plan, int dst_row_width, int transfer_width)
{
  int nsend       = src_plan->nslab_x;
  int nrecv_local = dst_plan->nslab_x;
  size_t block_elems = (size_t)transfer_width * transfer_width;

  int src_needs_pack   = (src_row_width != transfer_width);
  int dst_needs_unpack = (dst_row_width != transfer_width);

  /* Persistent, pre-allocated buffers (FDM_RedistPackBuf, sized to
   * comfortably cover whichever role -- send-side packing or
   * receive-side unpacking -- is needed; only ever one of the two at
   * once, confirmed by construction: this function is used exactly
   * twice per fdm_update_potential() call, once in each direction, and
   * each direction needs packing on exactly one side, never both) --
   * see fdm.h's header comment on why these are pre-allocated once,
   * early, rather than allocated-and-freed per call as an earlier
   * version did (a genuine crash, not a style preference). */
  fft_complex *send_buf = src_needs_pack ? FDM_RedistPackBuf : NULL;
  fft_complex *recv_buf = dst_needs_unpack ? FDM_RedistPackBuf : NULL;

  if(src_needs_pack)
    {
      for(int i = 0; i < nsend; i++)
        for(int j = 0; j < transfer_width; j++)
          memcpy(send_buf + (size_t)i * block_elems + (size_t)j * transfer_width,
                 src_rows + (size_t)i * src_row_width * src_row_width + (size_t)j * src_row_width,
                 (size_t)transfer_width * sizeof(fft_complex));
    }

  MPI_Request *requests = (MPI_Request *)FDM_RedistRequests;
  int nreq = 0;

  for(int i = 0; i < nsend; i++)
    {
      int gx = src_plan->slabstart_x + i;
      if(gx >= dst_plan->NgridX)
        continue;

      int dest = dst_plan->slab_to_task[gx];
      fft_complex *sendptr = src_needs_pack ? send_buf + (size_t)i * block_elems : src_rows + (size_t)i * block_elems;
      MPI_Isend(sendptr, block_elems * sizeof(fft_complex), MPI_BYTE, dest, gx, MPI_COMM_WORLD, &requests[nreq++]);
    }

  for(int i = 0; i < nrecv_local; i++)
    {
      int gx = dst_plan->slabstart_x + i;
      if(gx >= src_plan->NgridX)
        continue; /* outside the source mesh's range -- zero-padded region, no data to receive */

      int source = src_plan->slab_to_task[gx];
      fft_complex *recvptr = dst_needs_unpack ? recv_buf + (size_t)i * block_elems : dst_rows + (size_t)i * block_elems;
      MPI_Irecv(recvptr, block_elems * sizeof(fft_complex), MPI_BYTE, source, gx, MPI_COMM_WORLD, &requests[nreq++]);
    }

  MPI_Waitall(nreq, requests, MPI_STATUSES_IGNORE);

  if(dst_needs_unpack)
    {
      for(int i = 0; i < nrecv_local; i++)
        {
          int gx = dst_plan->slabstart_x + i;
          if(gx >= src_plan->NgridX)
            continue; /* zero-padded region: destination already zeroed by caller, nothing to unpack */

          for(int j = 0; j < transfer_width; j++)
            memcpy(dst_rows + (size_t)i * dst_row_width * dst_row_width + (size_t)j * dst_row_width,
                   recv_buf + (size_t)i * block_elems + (size_t)j * transfer_width,
                   (size_t)transfer_width * sizeof(fft_complex));
        }
    }
}

/*! \brief Builds the cached Green's function kernel on the doubled
 *  mesh. Called once, from fdm_allocate().
 *
 *  Plain, unsplit -1/r Green's function -- deliberately NOT gravity's
 *  own kernel (see the header comment on FDM_Kernel in fdm.h for why
 *  gravity's erfc-split, CIC-deconvolved kernel does not apply here).
 *
 *  r=0 self-term: uses a simple half-cell softening (r_soft=0.5*dx),
 *  NOT a precisely-derived analytic self-energy constant -- flagged
 *  explicitly as an approximate, "basic version first" choice. Only
 *  affects a single cell's self-contribution; a secondary refinement
 *  if it ever matters, not a first-order correctness concern for the
 *  overall solve.
 */
void fdm_poisson_kernel_init(void)
{
  int N2 = 2 * All.FDMGrid;
  double dx = All.FDMBoxSize / All.FDMGrid; /* original (non-padded) cell spacing, preserved in the padded mesh */

  my_slab_based_fft_init(&FDM_PoissonPlan, N2, N2, N2);

  size_t local_size = (size_t)FDM_PoissonPlan.nslab_x * N2 * N2;
  size_t max_local_size;
  MPI_Allreduce(&local_size, &max_local_size, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);

  FDM_Kernel           = (fft_complex *)mymalloc_movable(&FDM_Kernel, "FDM_Kernel", max_local_size * sizeof(fft_complex));
  FDM_PoissonWorkspace = (fft_complex *)mymalloc_movable(&FDM_PoissonWorkspace, "FDM_PoissonWorkspace",
                                                         max_local_size * sizeof(fft_complex));

  memset(FDM_Kernel, 0, max_local_size * sizeof(fft_complex));
  memset(FDM_PoissonWorkspace, 0, max_local_size * sizeof(fft_complex));

  /* --- FFTW plan creation for FDM_PoissonPlan --- this was missing
   * entirely in an earlier version of this function: my_slab_based_fft_init
   * only sets up layout bookkeeping (confirmed earlier, documented in
   * fdm.h), NOT the actual FFTW plan objects -- calling
   * my_slab_based_fft_c2c without this caused a segfault inside FFTW's
   * own execute call (uninitialized plan fields). Mirrors fdm_field.c's
   * fdm_allocate() plan creation exactly, just sized to N2 instead of N. */
  int ndim[1] = {N2};

#ifdef DOUBLEPRECISION_FFTW
  int alignflag = 0;
#else  /* #ifdef DOUBLEPRECISION_FFTW */
  int alignflag = FFTW_UNALIGNED;
#endif /* #ifdef DOUBLEPRECISION_FFTW #else */

#ifndef FFT_COLUMN_BASED
  int stride = N2;
#else  /* #ifndef FFT_COLUMN_BASED */
#error "FFT_COLUMN_BASED not yet handled for FDM -- see the same note in fdm_field.c's fdm_allocate()."
#endif /* #ifndef FFT_COLUMN_BASED #else */

  FDM_PoissonPlan.forward_plan_zdir = FFTW(plan_many_dft)(1, ndim, 1, FDM_Kernel, 0, 1, N2, FDM_PoissonWorkspace, 0, 1, N2,
                                                           FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  FDM_PoissonPlan.forward_plan_xdir =
      FFTW(plan_many_dft)(1, ndim, 1, FDM_Kernel, 0, stride, N2 * N2, FDM_PoissonWorkspace, 0, stride, N2 * N2, FFTW_FORWARD,
                          FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  FDM_PoissonPlan.forward_plan_ydir =
      FFTW(plan_many_dft)(1, ndim, 1, FDM_Kernel, 0, stride, N2 * N2, FDM_PoissonWorkspace, 0, stride, N2 * N2, FFTW_FORWARD,
                          FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  FDM_PoissonPlan.backward_plan_zdir = FFTW(plan_many_dft)(1, ndim, 1, FDM_Kernel, 0, 1, N2, FDM_PoissonWorkspace, 0, 1, N2,
                                                            FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  FDM_PoissonPlan.backward_plan_xdir =
      FFTW(plan_many_dft)(1, ndim, 1, FDM_Kernel, 0, stride, N2 * N2, FDM_PoissonWorkspace, 0, stride, N2 * N2, FFTW_BACKWARD,
                          FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  FDM_PoissonPlan.backward_plan_ydir =
      FFTW(plan_many_dft)(1, ndim, 1, FDM_Kernel, 0, stride, N2 * N2, FDM_PoissonWorkspace, 0, stride, N2 * N2, FFTW_BACKWARD,
                          FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  /* Build the real-space kernel, using the same "wraparound" indexing
   * convention as gravity's own kernel (index >= N2/2 represents a
   * negative offset) -- but with plain -1/r, physical units via dx. */
  for(int i = 0; i < FDM_PoissonPlan.nslab_x; i++)
    {
      int gx = FDM_PoissonPlan.slabstart_x + i;
      int wx = (gx < N2 / 2) ? gx : gx - N2;
      for(int j = 0; j < N2; j++)
        {
          int wy = (j < N2 / 2) ? j : j - N2;
          for(int k = 0; k < N2; k++)
            {
              int wz = (k < N2 / 2) ? k : k - N2;

              double x = wx * dx, y = wy * dx, z = wz * dx;
              double r = sqrt(x * x + y * y + z * z);

              size_t idx = (size_t)i * N2 * N2 + j * N2 + k;

              if(r > 0)
                FDM_Kernel[idx][0] = -1.0 / r;
              else
                FDM_Kernel[idx][0] = -1.0 / (0.5 * dx); /* simple half-cell softening, see header comment */

              FDM_Kernel[idx][1] = 0.0;
            }
        }
    }

  /* FFT the kernel once; result overwrites FDM_Kernel in place (matches
   * my_slab_based_fft_c2c's documented in-place-result behaviour). */
  my_slab_based_fft_c2c(&FDM_PoissonPlan, FDM_Kernel, FDM_PoissonWorkspace, 1);

  /* Persistent scratch buffers for fdm_update_potential()/
   * fdm_redistribute_rows() -- allocated HERE (called from
   * fdm_allocate(), which runs in begrun(), before P[]/SphP[] exist)
   * rather than per-call, fixing the movable-block corruption
   * described in fdm.h's header comment on these globals. */
  int N = All.FDMGrid;
  size_t local_size_N  = (size_t)FDM_plan.nslab_x * N * N;
  size_t local_size_N2 = (size_t)FDM_PoissonPlan.nslab_x * N2 * N2;

  FDM_RhoLocal  = (fft_complex *)mymalloc("FDM_RhoLocal", local_size_N * sizeof(fft_complex));
  FDM_RhoPadded = (fft_complex *)mymalloc("FDM_RhoPadded", local_size_N2 * sizeof(fft_complex));
  FDM_PhiLocal  = (fft_complex *)mymalloc("FDM_PhiLocal", local_size_N * sizeof(fft_complex));

  /* Sized to FDM_PoissonPlan.nslab_x*N*N -- covers both directions'
   * packing/unpacking need (see fdm_redistribute_rows itself: the only
   * side that ever needs packing/unpacking is whichever one touches
   * FDM_PoissonPlan's N2-wide rows, and FDM_PoissonPlan.nslab_x is the
   * relevant row count in both the forward and backward call). */
  size_t pack_buf_size = (size_t)FDM_PoissonPlan.nslab_x * N * N;
  FDM_RedistPackBuf = (fft_complex *)mymalloc("FDM_RedistPackBuf", pack_buf_size * sizeof(fft_complex));

  size_t requests_size = (size_t)(FDM_plan.nslab_x + FDM_PoissonPlan.nslab_x);
  FDM_RedistRequests = (void *)mymalloc("FDM_RedistRequests", requests_size * sizeof(MPI_Request));

  mpi_printf("FDM: Poisson kernel initialized (doubled mesh N=%d)\n", N2);
}

/*! \brief The potential-update step (Eq. 20c): rho = 4*pi*G*m*|psi|^2,
 *  zero-padded onto the doubled mesh, FFT'd, multiplied by the cached
 *  kernel, inverse FFT'd, then the [0,N) region extracted back into
 *  FDM_Potential.
 *
 *  No mean-density subtraction -- that correction is specifically for
 *  the periodic/cosmological case (solving for the peculiar potential
 *  relative to a uniform background); we solve the true, non-periodic
 *  Poisson equation directly, standard for an isolated system (same
 *  choice already made for gravity's own GRAVITY_NOT_PERIODIC path).
 */
void fdm_update_potential(void)
{
  int N  = All.FDMGrid;
  int N2 = 2 * All.FDMGrid;

  /* All.FDMMass is mc^2 in eV -- the source term needs the actual
   * mass, not the rest-mass energy. Converting once here rather than
   * threading eV-to-mass conversion through the Poisson-solve loop. */
  double m_grams  = (All.FDMMass * ELECTRONVOLT_IN_ERGS) / (CLIGHT * CLIGHT);
  double m_code   = m_grams / All.UnitMass_in_g;

  size_t local_size_N = (size_t)FDM_plan.nslab_x * N * N;

  /* Build rho = 4*pi*G*m*|psi|^2 on our LOCAL rows of the ORIGINAL
   * (non-padded) mesh, into the persistent FDM_RhoLocal buffer (see
   * fdm.h for why this is pre-allocated rather than allocated here). */
  double dx_cell = All.FDMBoxSize / N;
  double cell_volume = dx_cell * dx_cell * dx_cell;
  for(size_t idx = 0; idx < local_size_N; idx++)
    {
      double re = FDM_Psi[idx][0], im = FDM_Psi[idx][1];
      double dens = re * re + im * im;
      /* Mass density rho = m*|psi|^2 (Eq. 5) -- NOT multiplied by 4*pi*G
       * here. An earlier version of this function incorrectly included
       * a 4*pi*G factor at this stage: that's wrong because the -1/r
       * Green's function used below already has an implicit 4*pi
       * normalization built in (it solves grad^2(-1/r) = 4*pi*delta^3,
       * not delta^3 directly) -- multiplying density by 4*pi*G here AND
       * convolving with -1/r double-counts the 4*pi. The correct
       * G_Newton*dx^3 scaling is applied once, after the inverse FFT,
       * in fdm_update_potential() below -- see that comment for the
       * missing dx^3 (cell-volume) factor this also required, found
       * and derived numerically (not from memory) against a known
       * analytic uniform-sphere solution before trusting it. */
      double rho_wavefunction = m_code * dens;

      /* Phase 2a: stellar mass density, added into the SAME source
       * term the wavefunction's own density feeds. FDM_StarMassDensity
       * (fdm_particle_coupling.c, populated by fdm_deposit_star_mass(),
       * called once per fdm_advance_to_time() before this function's
       * own sub-cycling begins -- star positions don't change during
       * that sub-cycling, so redepositing on every inner fdm_step()
       * call would be wasted communication) holds CIC-deposited MASS
       * per cell, not density -- dividing by cell_volume converts it
       * to the same density units FDM_RhoLocal already uses, matching
       * how psi=sqrt(rho/m_code) itself defines that convention. */
      double rho_stars = FDM_StarMassDensity[idx] / cell_volume;

      FDM_RhoLocal[idx][0] = rho_wavefunction + rho_stars;
      FDM_RhoLocal[idx][1] = 0.0;
    }

  /* Redistribute into the padded mesh's decomposition. Zero the
   * destination buffer FIRST -- this is what actually implements the
   * zero-padding (the [N,2N) region of the padded mesh never receives
   * any data from the redistribution, since FDM_plan's grid only spans
   * [0,N)). */
  size_t local_size_N2 = (size_t)FDM_PoissonPlan.nslab_x * N2 * N2;
  memset(FDM_RhoPadded, 0, local_size_N2 * sizeof(fft_complex));

  /* NOTE: fdm_redistribute_rows moves whole rows of N*N (not N2*N2)
   * elements -- rho_local's rows are N*N (y,z both un-padded), and
   * they land at the START of each destination row in rho_padded
   * (indices 0..N*N-1 out of that row's N2*N2 total), which correctly
   * places them at the [0,N)x[0,N) corner of the y-z plane too, same
   * zero-padding convention applied consistently across all three
   * axes. */
  fdm_redistribute_rows(FDM_RhoLocal, &FDM_plan, N, FDM_RhoPadded, &FDM_PoissonPlan, N2, N);
  my_slab_based_fft_c2c(&FDM_PoissonPlan, FDM_RhoPadded, FDM_PoissonWorkspace, 1);

  for(size_t idx = 0; idx < local_size_N2; idx++)
    {
      double re = FDM_RhoPadded[idx][0] * FDM_Kernel[idx][0] - FDM_RhoPadded[idx][1] * FDM_Kernel[idx][1];
      double im = FDM_RhoPadded[idx][0] * FDM_Kernel[idx][1] + FDM_RhoPadded[idx][1] * FDM_Kernel[idx][0];
      FDM_RhoPadded[idx][0] = re;
      FDM_RhoPadded[idx][1] = im;
    }

  my_slab_based_fft_c2c(&FDM_PoissonPlan, FDM_RhoPadded, FDM_PoissonWorkspace, -1);

  double norm = 1.0 / ((double)N2 * N2 * N2);
  for(size_t idx = 0; idx < local_size_N2; idx++)
    FDM_RhoPadded[idx][0] *= norm; /* only the real part is physically meaningful from here on */

  /* Extract back: redistribute in reverse (src=padded plan, dst=original
   * plan) -- reusing the same helper, since it's symmetric in src/dst
   * roles by construction (each side independently computes its own
   * send/recv plan from the two slab_to_task[] arrays). */
  fdm_redistribute_rows(FDM_RhoPadded, &FDM_PoissonPlan, N2, FDM_PhiLocal, &FDM_plan, N, N);

  /* G_Newton * dx^3: the discrete FFT convolution (rho * kernel, via
   * FFT-multiply-inverseFFT) is an UNWEIGHTED discrete sum
   * (Sum_j rho_j * kernel_{i-j}), not the continuous convolution
   * integral Phi(x) = G_Newton * Integral[rho(x')/|x-x'| d^3x'] we
   * actually want -- approximating that integral as a Riemann sum
   * over grid cells needs an explicit dx^3 (cell volume) factor,
   * which is easy to miss since the FFT convolution theorem itself
   * doesn't include it. Both this and the (separate) incorrect 4*pi
   * factor removed above were found by comparing against a known
   * analytic solution (uniform-density sphere) numerically, not
   * derived from memory and trusted blindly -- see project validation
   * notes. */
  double dx = All.FDMBoxSize / All.FDMGrid;
  double scale = All.G * dx * dx * dx;

  for(size_t idx = 0; idx < local_size_N; idx++)
    FDM_Potential[idx] = FDM_PhiLocal[idx][0] * scale;
}

#endif /* #ifdef FDM */
