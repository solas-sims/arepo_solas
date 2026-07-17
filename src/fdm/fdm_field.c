/*!
 * \file        src/fdm/fdm_field.c
 * \brief       Allocation and FFT-plan setup for the FDM wavefunction
 *              field. Step 1 of the FDM module -- field representation
 *              and the distributed complex-to-complex FFT machinery
 *              the split-step integrator (kick-drift-kick, May &
 *              Springel 2021 Eqs. 17-20) will be built on top of.
 *              Kick/drift/potential-update integration itself is a
 *              later, separate piece -- this file only gets the field
 *              correctly allocated and the FFT plans correctly created.
 *
 *              KEY FINDING this is built on (verified by reading the
 *              actual implementations, not assumed): a genuine,
 *              distributed complex-to-complex FFT engine already
 *              exists in this codebase (my_slab_based_fft_c2c,
 *              src/gravity/pm/pm_mpi_fft.c) but nothing currently
 *              creates the complex-specific FFTW plan objects it
 *              expects -- my_slab_based_fft_init() only sets up
 *              distributed-layout bookkeeping (slab assignment across
 *              tasks), genuinely agnostic to real vs complex, and does
 *              NOT create any FFTW plan objects at all. The actual
 *              plan creation happens separately, in each gravity PM
 *              file's own init function (e.g. pm_init_nonperiodic() in
 *              pm_nonperiodic.c), which creates REAL-to-complex plans
 *              for the z-direction specifically (since gravity's
 *              density field starts real) and genuine complex-to-
 *              complex plans for y/x (operating on the already-complex
 *              intermediate result). Our own plan creation below
 *              mirrors that same structure, but uses complex-to-complex
 *              (fftw_plan_many_dft) for ALL three axes including z,
 *              since psi is complex from the start -- never real at
 *              any stage. This is confirmed necessary (not a
 *              simplification) by reading my_slab_based_fft_c2c's own
 *              execute code directly: it calls FFTW(execute_dft) (the
 *              generic complex-to-complex execute call) for z, y, AND
 *              x uniformly, using plan->NgridZ (the FULL grid
 *              dimension) throughout -- never plan->Ngridz (gravity's
 *              r2c-specific NgridZ/2+1 Hermitian-symmetry-halved
 *              field). Our own array sizing below matches this: no
 *              Ngrid2-style real/complex in-place padding anywhere,
 *              since there is no real-valued stage to pad for.
 */

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"
#include "fdm.h"

#ifdef FDM

fft_plan FDM_plan;
fft_complex *FDM_Psi;
fft_complex *FDM_PsiWorkspace;
fft_real *FDM_Potential;

/*! \brief Allocates the FDM wavefunction field and creates the
 *  distributed complex-to-complex FFT plans.
 *
 *  Mirrors pm_init_nonperiodic()'s structure (src/gravity/pm/
 *  pm_nonperiodic.c) closely, but:
 *    - sized off All.FDMGrid, not PMGRID/GRID (independent mesh
 *      resolution requirement, see project sizing notes),
 *    - genuinely complex-to-complex for all three axes (plan_many_dft
 *      throughout, not plan_many_dft_r2c/c2r for z),
 *    - no Ngrid2-style padding (no real-valued stage exists to pad
 *      for).
 */
void fdm_allocate(void)
{
  int N = All.FDMGrid;

  my_slab_based_fft_init(&FDM_plan, N, N, N);

  /* NOTE: my_slab_based_fft_init() sets FDM_plan.Ngridz = N/2+1 and
   * FDM_plan.Ngrid2 = 2*Ngridz internally (see its own implementation --
   * these fields exist unconditionally in the fft_plan struct, computed
   * assuming a real-to-complex convention). We deliberately do NOT use
   * FDM_plan.Ngridz/Ngrid2 anywhere in this module -- only
   * FDM_plan.NgridZ (the full, un-halved dimension), matching what
   * my_slab_based_fft_c2c's own execute code uses. Leaving those two
   * fields populated but unused is harmless; using them by mistake
   * elsewhere in this module would be a real bug, so this is flagged
   * here explicitly as something to double check in any future code
   * added to this file. */

  size_t local_size = (size_t)FDM_plan.nslab_x * N * N; /* complex elements, no padding */

  size_t max_local_size;
  MPI_Allreduce(&local_size, &max_local_size, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);

  FDM_Psi          = (fft_complex *)mymalloc_movable(&FDM_Psi, "FDM_Psi", max_local_size * sizeof(fft_complex));
  FDM_PsiWorkspace = (fft_complex *)mymalloc_movable(&FDM_PsiWorkspace, "FDM_PsiWorkspace", max_local_size * sizeof(fft_complex));
  FDM_Potential    = (fft_real *)mymalloc_movable(&FDM_Potential, "FDM_Potential", max_local_size * sizeof(fft_real));

  memset(FDM_Psi, 0, max_local_size * sizeof(fft_complex));
  memset(FDM_PsiWorkspace, 0, max_local_size * sizeof(fft_complex));
  memset(FDM_Potential, 0, max_local_size * sizeof(fft_real));

  /* --- FFTW plan creation --- */
  int ndim[1] = {N};

#ifdef DOUBLEPRECISION_FFTW
  int alignflag = 0;
#else  /* #ifdef DOUBLEPRECISION_FFTW */
  int alignflag = FFTW_UNALIGNED;
#endif /* #ifdef DOUBLEPRECISION_FFTW #else */

#ifndef FFT_COLUMN_BASED
  int stride = N; /* NgridZ, matching my_slab_based_fft_c2c's own ngridz=plan->NgridZ convention */
#else  /* #ifndef FFT_COLUMN_BASED */
#error "FFT_COLUMN_BASED not yet handled for FDM -- fdm_allocate() only implements the slab-based path so far. Needed if task count ever grows large enough that column-based decomposition becomes necessary (see pm_mpi_fft.c's own two-strategy design)."
#endif /* #ifndef FFT_COLUMN_BASED #else */

  /* z-direction: complex-to-complex, NOT plan_many_dft_r2c/c2r -- see
   * file header. This is the one genuine departure from
   * pm_nonperiodic.c's own template. */
  FDM_plan.forward_plan_zdir = FFTW(plan_many_dft)(1, ndim, 1, FDM_Psi, 0, 1, N, FDM_PsiWorkspace, 0, 1, N, FFTW_FORWARD,
                                                    FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  FDM_plan.forward_plan_xdir = FFTW(plan_many_dft)(1, ndim, 1, FDM_Psi, 0, stride, N * N, FDM_PsiWorkspace, 0, stride, N * N,
                                                    FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  FDM_plan.forward_plan_ydir = FFTW(plan_many_dft)(1, ndim, 1, FDM_Psi, 0, stride, N * N, FDM_PsiWorkspace, 0, stride, N * N,
                                                    FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  FDM_plan.backward_plan_zdir = FFTW(plan_many_dft)(1, ndim, 1, FDM_Psi, 0, 1, N, FDM_PsiWorkspace, 0, 1, N, FFTW_BACKWARD,
                                                     FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  FDM_plan.backward_plan_xdir = FFTW(plan_many_dft)(1, ndim, 1, FDM_Psi, 0, stride, N * N, FDM_PsiWorkspace, 0, stride, N * N,
                                                     FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  FDM_plan.backward_plan_ydir = FFTW(plan_many_dft)(1, ndim, 1, FDM_Psi, 0, stride, N * N, FDM_PsiWorkspace, 0, stride, N * N,
                                                     FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  mpi_printf("FDM: allocated wavefunction field, N=%d, local complex elements (max over tasks)=%zu\n", N, max_local_size);

  fdm_poisson_kernel_init(); /* allocated LAST here, so freed FIRST in fdm_free() below -- LIFO */
  fdm_gradient_allocate();   /* now allocated last -- freed first in fdm_free() */
  fdm_particle_coupling_allocate(); /* now allocated last -- freed first in fdm_free() */
}

void fdm_free(void)
{
  /* fdm_particle_coupling_allocate() ran LAST in fdm_allocate() --
   * freed FIRST here, per LIFO. */
  fdm_particle_coupling_free();

  /* fdm_gradient_allocate() ran LAST in fdm_allocate() -- freed FIRST
   * here, per LIFO. */
  fdm_gradient_free();

  /* FDM_RedistRequests/FDM_RedistPackBuf/FDM_PhiLocal/FDM_RhoPadded/
   * FDM_RhoLocal were allocated LAST (at the end of
   * fdm_poisson_kernel_init) -- freed FIRST here, per LIFO. */
  myfree(FDM_RedistRequests);
  myfree(FDM_RedistPackBuf);
  myfree(FDM_PhiLocal);
  myfree(FDM_RhoPadded);
  myfree(FDM_RhoLocal);

  /* Poisson-solve allocations were made LAST (fdm_poisson_kernel_init,
   * called at the end of fdm_allocate) -- freed FIRST here, per LIFO. */
  FFTW(destroy_plan)(FDM_PoissonPlan.backward_plan_ydir);
  FFTW(destroy_plan)(FDM_PoissonPlan.backward_plan_xdir);
  FFTW(destroy_plan)(FDM_PoissonPlan.backward_plan_zdir);
  FFTW(destroy_plan)(FDM_PoissonPlan.forward_plan_ydir);
  FFTW(destroy_plan)(FDM_PoissonPlan.forward_plan_xdir);
  FFTW(destroy_plan)(FDM_PoissonPlan.forward_plan_zdir);

  myfree_movable(FDM_PoissonWorkspace);
  myfree_movable(FDM_Kernel);

  myfree(FDM_PoissonPlan.first_slab_y_of_task);
  myfree(FDM_PoissonPlan.slabs_y_per_task);
  myfree(FDM_PoissonPlan.first_slab_x_of_task);
  myfree(FDM_PoissonPlan.slabs_x_per_task);
  myfree(FDM_PoissonPlan.slab_to_task);

  FFTW(destroy_plan)(FDM_plan.backward_plan_ydir);
  FFTW(destroy_plan)(FDM_plan.backward_plan_xdir);
  FFTW(destroy_plan)(FDM_plan.backward_plan_zdir);
  FFTW(destroy_plan)(FDM_plan.forward_plan_ydir);
  FFTW(destroy_plan)(FDM_plan.forward_plan_xdir);
  FFTW(destroy_plan)(FDM_plan.forward_plan_zdir);

  myfree_movable(FDM_Potential);
  myfree_movable(FDM_PsiWorkspace);
  myfree_movable(FDM_Psi);

  myfree(FDM_plan.first_slab_y_of_task);
  myfree(FDM_plan.slabs_y_per_task);
  myfree(FDM_plan.first_slab_x_of_task);
  myfree(FDM_plan.slabs_x_per_task);
  myfree(FDM_plan.slab_to_task);
}

#endif /* #ifdef FDM */
