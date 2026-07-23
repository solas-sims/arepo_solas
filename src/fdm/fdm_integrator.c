/*!
 * \file        src/fdm/fdm_integrator.c
 * \brief       Split-step kick-drift-kick integrator for the FDM
 *              wavefunction (May & Springel 2021, Eqs. 17-20).
 *              This file: the drift step only (Eq. 19's middle term)
 *              -- kick and potential-update are separate, later pieces.
 *              Non-cosmological (a=1 throughout) for the current
 *              central-galaxy target; the cosmological a(t)-dependent
 *              terms are a documented simplification, not yet needed.
 */

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"
#include "fdm.h"

#ifdef FDM

/*! \brief Computes hbar/m in code units.
 *
 *  All.FDMMass is mc^2 in eV (the boson rest-mass energy, matching the
 *  literature's own convention for quoting FDM particle masses -- NOT
 *  the mass itself). Converts through CGS (erg*s / g = cm^2/s) to code
 *  units (code_length^2/code_time), cross-checked against the source
 *  paper's own quoted lambda_dB value (1.21 kpc at mc^2=1e-22 eV,
 *  v=100 km/s) before trusting this -- see project validation notes.
 *  PLANCK (allvars.h) is h, not hbar=h/(2 pi) -- divided out here.
 */
double fdm_hbar_over_m_code(void)
{
  double hbar_erg_s = PLANCK / (2.0 * M_PI);
  double m_grams    = (All.FDMMass * ELECTRONVOLT_IN_ERGS) / (CLIGHT * CLIGHT);
  double hbar_over_m_cgs = hbar_erg_s / m_grams; /* cm^2/s */

  double UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;

  return hbar_over_m_cgs * UnitTime_in_s / (All.UnitLength_in_cm * All.UnitLength_in_cm);
}

/*! \brief Converts a grid index (0..N-1) to its physical wavenumber,
 *  using the standard FFT frequency convention (negative frequencies
 *  wrap around above the Nyquist index). L is the code-units box size.
 */
static inline double fdm_k_of_index(int idx, int N, double L)
{
  int wrapped = (idx <= N / 2) ? idx : idx - N;
  return 2.0 * M_PI * wrapped / L;
}

/*! \brief The drift step: free-particle (kinetic-term-only) evolution
 *  of the wavefunction over a half (or full) timestep dt, applied
 *  entirely in Fourier space -- exp(-i*(hbar/m)*dt*k^2/2) per Eq. 19's
 *  central drift factor (a=1, non-cosmological).
 *
 *  Index-to-k mapping used below (index = i_local_y*N*N + kx*N + kz,
 *  ky_global = FDM_plan.slabstart_y + i_local_y) was verified directly
 *  against the actual post-transform array layout using a single-mode
 *  test (see fdm_test_single_mode.c) BEFORE being trusted here -- not
 *  assumed from reading my_slab_transpose alone, since getting this
 *  wrong would silently miscompute every mode's phase without any
 *  obvious symptom.
 */
void fdm_drift(double dt)
{
  int    N   = All.FDMGrid;
  double L   = All.FDMBoxSize;
  double hbar_over_m = fdm_hbar_over_m_code();

  my_slab_based_fft_c2c(&FDM_plan, FDM_Psi, FDM_PsiWorkspace, 1);

  for(int i = 0; i < FDM_plan.nslab_y; i++)
    {
      int    gy = FDM_plan.slabstart_y + i;
      double ky = fdm_k_of_index(gy, N, L);

      for(int kx_idx = 0; kx_idx < N; kx_idx++)
        {
          double kx = fdm_k_of_index(kx_idx, N, L);

          for(int kz_idx = 0; kz_idx < N; kz_idx++)
            {
              double kz = fdm_k_of_index(kz_idx, N, L);
              double k2 = kx * kx + ky * ky + kz * kz;

              double phase = -hbar_over_m * dt * k2 / 2.0;
              double c = cos(phase), s = sin(phase);

              size_t idx = (size_t)i * N * N + (size_t)kx_idx * N + kz_idx;
              double re  = FDM_Psi[idx][0];
              double im  = FDM_Psi[idx][1];

              FDM_Psi[idx][0] = re * c - im * s;
              FDM_Psi[idx][1] = re * s + im * c;
            }
        }
    }

  my_slab_based_fft_c2c(&FDM_plan, FDM_Psi, FDM_PsiWorkspace, -1);

  double norm        = 1.0 / ((double)N * N * N);
  size_t local_size = (size_t)FDM_plan.nslab_x * N * N;
  for(size_t idx = 0; idx < local_size; idx++)
    {
      FDM_Psi[idx][0] *= norm;
      FDM_Psi[idx][1] *= norm;
    }
}

/*! \brief The kick step: psi(x) *= exp(-i*(m/hbar)*dt*Phi(x)), pointwise
 *  in real space, using the current contents of FDM_Potential.
 *
 *  Unlike the drift step, this has no FFT and no k-space indexing at
 *  all -- every element is transformed independently using only its
 *  own local value, so there is no risk of the kind of index-mapping
 *  error the drift step needed a dedicated single-mode test to catch.
 *  A pointwise |psi| conservation check (stronger than the drift
 *  step's total-norm check, since it holds at every single grid point
 *  individually, not just in aggregate) is the natural validation
 *  here.
 *
 *  dt is used directly, no implicit halving -- see the header comment
 *  in fdm.h for the half-kick/merged-kick calling convention.
 */
void fdm_kick(double dt)
{
  double m_over_hbar = 1.0 / fdm_hbar_over_m_code();

  int    N          = All.FDMGrid;
  size_t local_size = (size_t)FDM_plan.nslab_x * N * N;

  for(size_t idx = 0; idx < local_size; idx++)
    {
      double phase = -m_over_hbar * dt * FDM_Potential[idx];
      double c = cos(phase), s = sin(phase);

      double re = FDM_Psi[idx][0];
      double im = FDM_Psi[idx][1];

      FDM_Psi[idx][0] = re * c - im * s;
      FDM_Psi[idx][1] = re * s + im * c;
    }
}

/*! \brief The timestep criterion (Eq. 21), non-cosmological (a=1):
 *    dt < min[(4/3pi)*(hbar/m)*dx^2, 2*pi*(hbar/m)/|Phi_max|]
 *
 *  The first (drift) term bounds the phase accumulated by the largest
 *  k-mode the mesh can represent, so it doesn't alias; the second
 *  (kick) term bounds the phase accumulated from the potential itself.
 *  Requires FDM_Potential to already be populated (via
 *  fdm_update_potential()) -- callers should ensure that before calling
 *  this, same as fdm_kick().
 *
 *  |Phi_max| is a global (MPI-reduced) maximum over the whole
 *  distributed FDM_Potential array, not just this task's local rows.
 */
double fdm_get_timestep(void)
{
  int    N  = All.FDMGrid;
  double dx = All.FDMBoxSize / N;
  double hbar_over_m = fdm_hbar_over_m_code();

  size_t local_size = (size_t)FDM_plan.nslab_x * N * N;
  double phi_max_local = 0.0;
  for(size_t idx = 0; idx < local_size; idx++)
    {
      double a = fabs(FDM_Potential[idx]);
      if(a > phi_max_local)
        phi_max_local = a;
    }

  double phi_max_global;
  MPI_Allreduce(&phi_max_local, &phi_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  double dt_drift = (4.0 / (3.0 * M_PI)) * hbar_over_m * dx * dx;

  if(phi_max_global <= 0.0)
    return dt_drift; /* no meaningful kick constraint yet (e.g. before any potential has been computed) */

  double dt_kick = 2.0 * M_PI * hbar_over_m / phi_max_global;

  return (dt_drift < dt_kick) ? dt_drift : dt_kick;
}

/*! \brief One full split-step kick-drift-kick timestep (Eqs. 17-20),
 *  non-cosmological (a=1): half-kick, full drift, potential update
 *  (from the newly-drifted psi), half-kick.
 *
 *  DELIBERATELY does not implement the paper's consecutive-half-kick
 *  merging optimization (combining the trailing half-kick of one step
 *  with the leading half-kick of the next into a single full-dt kick
 *  when calling this repeatedly in a loop) -- that's a performance
 *  optimization (saves one kick's worth of redundant computation per
 *  step boundary), not a correctness requirement, and keeping each
 *  call to this function fully self-contained (valid on its own,
 *  simple to test in isolation) matches this project's "basic working
 *  version first" approach throughout. A documented, deferred
 *  optimization, not an oversight.
 *
 *  PRECONDITION: FDM_Potential must already be valid (from a prior
 *  fdm_update_potential() call) before the FIRST call to this function
 *  in a simulation -- the opening half-kick uses whatever potential is
 *  currently there. Every call after the first is self-sufficient
 *  (the potential gets refreshed, from the newly-drifted psi, before
 *  the closing half-kick, and that refreshed potential is exactly
 *  what the NEXT call's opening half-kick correctly wants to use).
 */
void fdm_step(double dt)
{
  fdm_kick(dt / 2.0);
  fdm_drift(dt);
  fdm_update_potential();
  fdm_kick(dt / 2.0);
}

/*! \brief FDM's own continuous physical-time clock, in the same code
 *  time units fdm_step()/fdm_get_timestep() operate in -- distinct
 *  from All.Ti_Current, which is Arepo's integer timebin bookkeeping,
 *  not a physical time directly. Initialized to All.TimeBegin (see
 *  fdm_advance_to_time()'s first-call handling), matching the same
 *  zero-point convention All.Ti_Current=0 uses.
 */
static double FDM_CurrentTime;
static int    FDM_TimeInitialized = 0;

/*! \brief Converts an integer timebin value to physical code time.
 *
 *  Non-cosmological only (All.ComovingIntegrationOn == 0) -- matches
 *  the same a=1 simplification documented throughout this module
 *  (fdm_drift, fdm_kick, fdm_get_timestep). Verified directly against
 *  how All.Timebase_interval itself gets set up (init.c): for a
 *  non-cosmological run, Timebase_interval = (TimeMax-TimeBegin)/
 *  TIMEBASE (linear), with Ti_Current=0 corresponding to time=TimeBegin
 *  -- so time = TimeBegin + Ti_Current*Timebase_interval directly, no
 *  log-scale-factor complication (that convention is cosmological-only
 *  and would need separate, not-yet-implemented handling here).
 */
static double fdm_time_of_ti(integertime ti)
{
  if(All.ComovingIntegrationOn)
    terminate("FDM: fdm_time_of_ti() only implements the non-cosmological conversion -- All.ComovingIntegrationOn=1 needs "
              "separate, not-yet-implemented handling (same a=1 simplification documented throughout this module).\n");

  return All.TimeBegin + ti * All.Timebase_interval;
}

/*! \brief Advances the FDM wavefunction, via internal fdm_step() sub-
 *  cycling, to catch up with target_ti (Arepo's own current/next
 *  sync-point time). Called once per outer main-loop iteration, right
 *  after find_next_sync_point() updates All.Ti_Current -- see the call
 *  site in run.c.
 *
 *  Phase 1 (DM-only, no baryonic coupling yet): FDM evolves on its own,
 *  fully independent clock, sub-cycling as many internal fdm_step()
 *  calls as its own timestep criterion demands to reach target_ti,
 *  landing exactly on it (not overshooting) rather than participating
 *  in Arepo's hierarchical timebin system at all -- there is no
 *  coupling to lose synchronization with yet. This will need revisiting
 *  once particles need to feel FDM's potential DURING FDM's own sub-
 *  cycling, not just at sync-point boundaries (Phase 2a+).
 *
 *  First call bootstraps: sets FDM_CurrentTime = All.TimeBegin (the
 *  same zero-point All.Ti_Current=0 uses) and establishes the initial
 *  potential via fdm_update_potential() -- fdm_step()'s own documented
 *  precondition, satisfied here rather than left to the caller.
 */
void fdm_advance_to_time(integertime target_ti)
{
  double target_time = fdm_time_of_ti(target_ti);

  /* Phase 2a: deposit star mass ONCE per call (i.e. once per outer
   * sync point), before ANY potential update this call performs --
   * including the bootstrap one below, so even the very first
   * potential solve already includes the stellar contribution, not
   * just later ones. Star positions don't change during this
   * function's own fdm_step() sub-cycling below (stars only move via
   * Arepo's own gravity/timestep machinery, at the OUTER sync-point
   * cadence) -- so one deposit here is both correct and sufficient;
   * redepositing on every inner fdm_step() call would just be
   * repeated communication for an unchanged mass distribution. This
   * is the same basic-version-first limitation already flagged below
   * (fdm_interpolate_to_stars() call) and in this function's own
   * header comment: stars feel FDM (and vice versa) only at sync-point
   * boundaries, not during FDM's finer internal sub-cycling. */
  fdm_deposit_star_mass();

  if(!FDM_TimeInitialized)
    {
      FDM_CurrentTime     = All.TimeBegin;
      FDM_TimeInitialized = 1;
      fdm_update_potential();
    }

  while(FDM_CurrentTime < target_time)
    {
      double dt = fdm_get_timestep();
      if(FDM_CurrentTime + dt > target_time)
        dt = target_time - FDM_CurrentTime; /* land exactly on target_time, don't overshoot */

      fdm_step(dt);
      FDM_CurrentTime += dt;
    }

  /* fdm_compute_force() computes FDM_ForceX/Y/Z from the potential
   * fdm_update_potential() just left in FDM_Potential -- needed here
   * because fdm_interpolate_to_stars() below reads those force arrays
   * directly. fdm_step()'s own internal fdm_update_potential() calls
   * don't need a matching force computation: fdm_kick() (the FDM
   * field's own evolution) uses the potential directly (the standard
   * Schrodinger-Poisson kick, psi *= exp(-i*V*dt/hbar)), never its
   * gradient -- the force is only ever consumed by the star coupling
   * below, so computing it once here, from the final, fully-updated
   * potential after all sub-cycling above has completed, is both
   * correct and sufficient. A real, confirmed regression: this call
   * was missing entirely from this function -- found from stars
   * showing a genuinely nonzero, reasonable FDM_Potential but an
   * exactly-zero force at their own position, which is not physically
   * consistent for a real, spatially-varying potential (a nonzero,
   * varying potential must have a nonzero gradient almost everywhere)
   * -- confirmed by grep showing fdm_compute_force() was called only
   * from this project's own standalone test files, never from the
   * actual simulation's own call path. */
  fdm_compute_force();

  /* Phase 2a: interpolate the now-current potential/force to stars,
   * once, after sub-cycling has fully caught up to target_ti -- ready
   * for Arepo's own gravity accumulation to pick up FDM_StarResult
   * before the next kick. Same basic-version-first caveat as the
   * deposit call above: this reflects the potential at target_ti, not
   * anything finer-grained during the sub-cycling in between. */
  fdm_interpolate_to_stars();

  /* Output, rate-limited to roughly All.TimeBetSnapshot -- the same
   * cadence particle snapshots use, tracked independently rather than
   * reusing All.Ti_nextoutput directly (that variable has ALREADY been
   * advanced past the current time by create_snapshot_if_desired(),
   * which runs earlier in the same run.c loop iteration -- checking it
   * here would never fire). Genuinely needed, not a nicety: writing a
   * full field snapshot every single sync point (this function's own
   * previous behaviour) accumulates fast -- N=128 real+imag doubles is
   * ~34MB per snapshot, and several hundred sync points, easily
   * reached during a real run, means tens of GB -- confirmed as the
   * actual cause of a real "unable to write dataset" HDF5 failure,
   * traced to exhausted disk quota, not any actual data corruption. */
  {
    static double fdm_next_output_time = -1.0;
    static int fdm_output_counter      = 0;

    if(fdm_next_output_time < 0.0)
      fdm_next_output_time = All.TimeBegin;

    if(FDM_CurrentTime >= fdm_next_output_time)
      {
        char fname[MAXLEN_PATH];
        snprintf(fname, MAXLEN_PATH, "%s/fdm_field_%03d.hdf5", All.OutputDir, fdm_output_counter);
        fdm_write_field(fname);
        fdm_output_counter++;

        fdm_next_output_time += All.TimeBetSnapshot;
      }
  }
}

#endif /* #ifdef FDM */
