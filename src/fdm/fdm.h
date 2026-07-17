#ifndef FDM_H
#define FDM_H

#include "../main/allvars.h"

#ifdef FDM

#ifndef PMGRID
#error "FDM requires PMGRID to be set in Config.sh -- fft_real/fft_complex/fft_plan (dtypes.h) and the FFTW3 linkage are gated behind it, and we reuse them rather than duplicating. This doesn't force FDM to use gravity's own PM solver for anything -- it just needs the same type definitions and library linkage."
#endif

/*! \brief Distributed FFT plan for the FDM wavefunction field.
 *
 *  A dedicated fft_plan instance, separate from gravity's own PM solver
 *  (myplan, file-static in pm_nonperiodic.c/pm_periodic.c) -- FDM's mesh
 *  resolution requirement is set independently (by the de Broglie
 *  wavelength / structural resolution needs, see the sizing discussion
 *  in project notes), not tied to PMGRID's own long-range-gravity
 *  resolution.
 *
 *  IMPORTANT: this is a genuine complex-to-complex transform throughout
 *  -- psi is complex at every stage, unlike gravity's PM solver which
 *  starts from a real density field and uses a real-to-complex (r2c)
 *  transform with the associated NgridZ/2+1 Hermitian-symmetry halving
 *  and in-place real/complex reinterpretation padding (Ngrid2). None of
 *  that applies here. FDM_plan.NgridZ (the FULL grid size, not the r2c
 *  Ngridz=NgridZ/2+1 field) is what all our own array sizing must use,
 *  matching what my_slab_based_fft_c2c's execute code itself uses
 *  internally (confirmed by reading its implementation directly, not
 *  assumed).
 *
 *  The distributed layout/transpose machinery (my_slab_based_fft_init,
 *  my_slab_transpose) is genuinely shared/reused as-is -- it is layout
 *  bookkeeping only, agnostic to real vs complex. What we DO need to
 *  build ourselves is the actual FFTW plan creation (populating
 *  FDM_plan.forward_plan_zdir etc. with genuine complex-to-complex
 *  fftw_plan_many_dft plans) -- mirroring pm_nonperiodic.c's own
 *  plan-creation code structurally, but using plan_many_dft for ALL
 *  three axes (including z), not plan_many_dft_r2c/c2r for z as
 *  gravity's real-input case does.
 */
extern fft_plan FDM_plan;

/*! \brief The wavefunction field itself, psi_c (comoving-scaled complex
 *  field per Eq. 6-7 of May & Springel 2021), distributed across tasks
 *  per FDM_plan's slab decomposition. Sized NgridZ (full grid size,
 *  see note above), NOT the r2c-halved Ngridz -- no padding. */
extern fft_complex *FDM_Psi;

/*! \brief Scratch workspace of the same size as FDM_Psi, needed by the
 *  transpose/FFT execution functions (mirrors gravity's own forcegrid/
 *  workspace role). */
extern fft_complex *FDM_PsiWorkspace;

/*! \brief The gravitational potential, real-valued, on the SAME
 *  real-space grid layout as FDM_Psi (nslab_x local rows of x, full
 *  y, full z -- i.e. the layout FDM_Psi is in before any forward
 *  transform, or after the drift step's inverse transform). Computed
 *  by the (separate, not-yet-implemented) potential-update step from
 *  |psi|^2 via a Poisson solve; consumed by fdm_kick() below. */
extern fft_real *FDM_Potential;

void fdm_allocate(void);
void fdm_free(void);
double fdm_hbar_over_m_code(void);
void fdm_drift(double dt);

/*! \brief Applies the kick step: psi(x) *= exp(-i*(m/hbar)*dt*Phi(x)),
 *  pointwise in real space, using the CURRENT contents of
 *  FDM_Potential. dt is used directly (no implicit halving) -- the
 *  caller passes dt/2 for a boundary half-kick or the full merged dt
 *  for an interior kick, matching the paper's own kick-merging
 *  optimization (Eq. 19's note on combining consecutive half-kicks).
 */
void fdm_kick(double dt);

/*! \brief Timestep criterion (Eq. 21, non-cosmological a=1):
 *  dt < min[(4/3pi)*(hbar/m)*dx^2, 2*pi*(hbar/m)/|Phi_max|]. Requires
 *  FDM_Potential to already be populated. */
double fdm_get_timestep(void);

/*! \brief One full kick-drift-kick timestep (Eqs. 17-20): half-kick,
 *  drift, potential update, half-kick. See fdm_integrator.c for the
 *  precondition (FDM_Potential must be valid before the first call)
 *  and the note on the deliberately-not-implemented consecutive-
 *  half-kick merging optimization.
 */
void fdm_step(double dt);

/*! \brief Advances FDM (via internal fdm_step sub-cycling) to catch up
 *  with target_ti. Called once per main-loop iteration, right after
 *  find_next_sync_point() -- see run.c and fdm_integrator.c for the
 *  full reasoning (non-cosmological only for now; independent-clock
 *  design specific to Phase 1's lack of baryonic coupling).
 */
void fdm_advance_to_time(integertime target_ti);

/*! \brief Gathers FDM_Psi to task 0 and writes it to an HDF5 file
 *  (see fdm_io.c for the full format specification -- deliberately the
 *  same format an external IC generator must produce for
 *  fdm_read_field() to consume). */
void fdm_write_field(const char *fname);

/*! \brief Reads an HDF5 file in fdm_write_field()'s format and scatters
 *  it into FDM_Psi. Checks the file's N/FDMMass against the current
 *  run's own values, terminating on mismatch rather than silently
 *  proceeding. */
void fdm_read_field(const char *fname);

/*! \brief Separate, doubled-mesh (2*All.FDMGrid per dimension) FFT plan
 *  used ONLY for the Poisson solve -- the standard Hockney & Eastwood
 *  zero-padding technique for isolated (non-periodic) boundary
 *  conditions via FFT. Distinct from FDM_plan (which stays at
 *  All.FDMGrid, used for psi's own drift/kick).
 *
 *  Deliberately reuses our existing, already-validated c2c machinery
 *  (treating the real-valued padded density/kernel as complex with
 *  zero imaginary part) rather than adding a second, r2c-specific
 *  plan-creation path -- costs 2x memory/compute versus a "proper"
 *  r2c implementation, but is dramatically simpler and lower-risk
 *  given everything already validated for c2c specifically. "Basic
 *  working version first" -- r2c is a documented future optimization,
 *  not needed now.
 *
 *  NOTE ON COST: the doubled mesh means 8x the memory/compute of the
 *  base FDM_plan mesh for this step specifically -- a real, separate
 *  cost consideration from the base-mesh sizing done earlier, worth
 *  remembering once this is run for real.
 */
extern fft_plan FDM_PoissonPlan;

/*! \brief The FFT'd Green's function kernel for the isolated Poisson
 *  solve, on the doubled mesh. Built ONCE (fdm_poisson_kernel_init),
 *  cached, reused every potential-update call.
 *
 *  DELIBERATELY NOT a copy of gravity's own kernel (pm_nonperiodic.c's
 *  pm_setup_nonperiodic_kernel()) -- that kernel is long-range-only
 *  (an Ewald/TreePM short-range/long-range split, fac=1-erfc(u)/r,
 *  since gravity's tree already handles short-range forces) and
 *  includes a CIC deconvolution correction (compensating for
 *  Cloud-In-Cell particle-mass-deposition error). Neither applies to
 *  us: we need the FULL, unsplit potential from |psi|^2 (no separate
 *  short-range handling of the wavefunction's self-gravity exists),
 *  and our density is already exactly on the grid, no particle
 *  deposition/CIC involved. Same general TECHNIQUE as gravity's kernel
 *  (build in real space, FFT once, cache), genuinely different
 *  content.
 */
extern fft_complex *FDM_Kernel;

/*! \brief Scratch workspace for the Poisson-solve FFTs, doubled-mesh
 *  sized, mirroring FDM_PsiWorkspace's role for FDM_plan. */
extern fft_complex *FDM_PoissonWorkspace;

/*! \brief Builds FDM_Kernel -- the one-time, cached Green's function
 *  setup. Call once during fdm_allocate(), not per-timestep.
 */
void fdm_poisson_kernel_init(void);

/*! \brief The potential-update step: computes FDM_Potential from the
 *  current FDM_Psi via 4*pi*G*m*|psi|^2 -> zero-pad -> FFT -> multiply
 *  kernel -> inverse FFT -> extract, per Eq. 20c (no mean-density
 *  subtraction -- that's specifically for the periodic/cosmological
 *  case; we solve the true non-periodic Poisson equation directly, as
 *  is standard for an isolated system).
 */
void fdm_update_potential(void);

/*! \brief Force field (-grad(FDM_Potential)), same real-space
 *  distributed layout as FDM_Potential itself. Computed by
 *  fdm_compute_force() (src/fdm/fdm_gradient.c) via 4th-order finite
 *  differencing -- see that file's header for the full stencil
 *  derivation, the transpose-based handling of X (which needs
 *  neighbouring slabs potentially on a different task), and the
 *  boundary fallback for isolated (non-periodic) domains.
 */
extern fft_real *FDM_ForceX;
extern fft_real *FDM_ForceY;
extern fft_real *FDM_ForceZ;

void fdm_gradient_allocate(void);
void fdm_gradient_free(void);
void fdm_compute_force(void);

/*! \brief Per-particle interpolated FDM potential/force, filled by
 *  fdm_interpolate_to_stars() (src/fdm/fdm_particle_coupling.c).
 *  Unlike every other array in this module (fixed-size mesh data,
 *  allocated once), this is PER-PARTICLE -- sized to All.MaxPart and
 *  grown in lockstep with P[] itself, mirroring DMSP[]'s own
 *  established pattern exactly (src/utils/allocate.c): no independent
 *  reallocate_memory_* function needed, since it's always tied 1:1 to
 *  All.MaxPart the same way DMSP[] is. */
typedef struct
{
  double Potential;
  double ForceX;
  double ForceY;
  double ForceZ;
} fdm_star_result;

extern fdm_star_result *FDM_StarResult;

void fdm_interpolate_to_stars(void);

/*! \brief Stellar mass density deposited by fdm_deposit_star_mass()
 *  (src/fdm/fdm_particle_coupling.c) -- read directly by
 *  fdm_update_potential() (fdm_poisson.c) to add the stellar
 *  contribution to the Poisson solve's source term. Same real-space
 *  distributed layout as FDM_Potential. */
extern fft_real *FDM_StarMassDensity;

void fdm_particle_coupling_allocate(void);
void fdm_particle_coupling_free(void);
void fdm_deposit_star_mass(void);

/*! \brief Persistent scratch buffers for fdm_update_potential() and
 *  fdm_redistribute_rows(), allocated ONCE in fdm_poisson_kernel_init()
 *  (called from fdm_allocate(), which runs in begrun() -- BEFORE
 *  particle arrays P[]/SphP[] get allocated during IC reading) rather
 *  than allocated-and-freed on every call as an earlier version did.
 *
 *  This isn't a performance optimization -- it fixes a genuine, real
 *  crash: mymalloc.c's own documentation states "For a movable block
 *  to be successfully shifted it is required that all the subsequent
 *  allocated blocks are movable." FDM_Psi/FDM_Potential/FDM_Kernel/
 *  FDM_PoissonWorkspace are all mymalloc_movable-allocated once at
 *  startup; P[]/SphP[] (Arepo's own particle arrays) are ALSO movable
 *  and get resized during domain decomposition throughout the run. An
 *  earlier version of fdm_update_potential()/fdm_redistribute_rows()
 *  allocated their scratch buffers with plain (non-movable) mymalloc()
 *  repeatedly, DURING the run -- meaning these non-movable buffers
 *  existed AFTER P[]/SphP[] in the allocation stack whenever
 *  fdm_update_potential() happened to be executing, which would
 *  violate the above requirement if P[]/SphP[] needed to shift at that
 *  exact moment, corrupting the movable-block bookkeeping (observed as
 *  FDM_PoissonWorkspace's own global pointer becoming NULL mid-run,
 *  crashing inside my_slab_based_fft_c2c on a later call -- confirmed
 *  via an lldb backtrace on Apple Silicon; a real, live bug, not a
 *  theoretical one). Allocating these ONCE, early, before P[]/SphP[]
 *  exist, means their presence can never violate this constraint for
 *  any later P[]/SphP[] resize, regardless of whether these particular
 *  buffers are themselves movable or not (they aren't -- plain
 *  mymalloc(), since they're never resized or moved individually,
 *  simplest choice given they only need to exist once, early, and be
 *  freed once, at shutdown).
 *
 *  Sized to comfortably cover BOTH directions fdm_redistribute_rows()
 *  is used in (forward: FDM_plan -> FDM_PoissonPlan; backward: the
 *  reverse) -- see fdm_poisson.c for the exact sizing reasoning.
 */
extern fft_complex *FDM_RhoLocal;
extern fft_complex *FDM_RhoPadded;
extern fft_complex *FDM_PhiLocal;
extern fft_complex *FDM_RedistPackBuf;
extern void        *FDM_RedistRequests; /* MPI_Request*, void* here to avoid requiring mpi.h in every fdm.h includer */

#endif /* #ifdef FDM */

#endif /* #ifndef FDM_H */
