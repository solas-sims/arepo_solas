/*!
 * \file        src/fdm/fdm_test_drift.c
 * \brief       Validates fdm_drift(): (1) total norm (sum|psi|^2)
 *              conservation -- a genuine physical invariant of free
 *              Schrodinger evolution, index-order-agnostic; (2) for a
 *              single known Fourier mode, the exact predicted phase
 *              shift exp(-i*(hbar/m)*dt*k^2/2) with unchanged
 *              magnitude -- this DOES rely on the index-to-k mapping,
 *              already independently verified in fdm_test_single_mode.c
 *              before fdm_drift.c was written to use it.
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"
#include "fdm.h"

void endrun(void)
{
  MPI_Finalize();
  exit(1);
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  All.MaxMemSize = 500;
  mymalloc_init();

  All.FDMGrid    = 16;
  All.FDMBoxSize = 10.0; /* code length units, e.g. kpc */
  All.FDMMass    = 1e-22; /* eV */
  All.UnitLength_in_cm        = 3.085678e21; /* 1 kpc */
  All.UnitMass_in_g           = 1.989e43;    /* 1e10 Msun */
  All.UnitVelocity_in_cm_per_s = 1e5;         /* 1 km/s */

  int N = All.FDMGrid;
  int kx0 = 3;
  double dt = 0.01;

  fdm_allocate();

  size_t local_size = (size_t)FDM_plan.nslab_x * N * N;

  /* Single mode along X, same construction as fdm_test_single_mode.c */
  double norm_before_local = 0.0;
  for(int i = 0; i < FDM_plan.nslab_x; i++)
    {
      int gx = FDM_plan.slabstart_x + i;
      double phase = 2 * M_PI * kx0 * gx / N;
      for(int j = 0; j < N; j++)
        for(int k = 0; k < N; k++)
          {
            size_t idx = (size_t)i * N * N + j * N + k;
            FDM_Psi[idx][0] = cos(phase);
            FDM_Psi[idx][1] = sin(phase);
            norm_before_local += FDM_Psi[idx][0] * FDM_Psi[idx][0] + FDM_Psi[idx][1] * FDM_Psi[idx][1];
          }
    }

  double norm_before_global;
  MPI_Allreduce(&norm_before_local, &norm_before_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  fdm_drift(dt);

  double norm_after_local = 0.0;
  for(size_t idx = 0; idx < local_size; idx++)
    norm_after_local += FDM_Psi[idx][0] * FDM_Psi[idx][0] + FDM_Psi[idx][1] * FDM_Psi[idx][1];

  double norm_after_global;
  MPI_Allreduce(&norm_after_local, &norm_after_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  double norm_rel_err = fabs(norm_after_global - norm_before_global) / norm_before_global;

  /* Analytic prediction: drift should leave a pure-X plane wave as a
   * pure-X plane wave (mode is an eigenstate of the kinetic operator),
   * unchanged magnitude, with an additional GLOBAL phase
   * exp(-i*(hbar/m)*dt*kx^2/2) applied uniformly -- so the field
   * should still look like exp(i*(original_phase + extra_phase)) at
   * every point. */
  double hbar_over_m = fdm_hbar_over_m_code();
  double kx_phys = 2 * M_PI * kx0 / All.FDMBoxSize;
  double extra_phase = -hbar_over_m * dt * (kx_phys * kx_phys) / 2.0;

  double max_phase_err = 0.0;
  for(int i = 0; i < FDM_plan.nslab_x; i++)
    {
      int gx = FDM_plan.slabstart_x + i;
      double original_phase = 2 * M_PI * kx0 * gx / N;
      double predicted_phase = original_phase + extra_phase;
      double pre = cos(predicted_phase), pim = sin(predicted_phase);

      for(int j = 0; j < N; j++)
        for(int k = 0; k < N; k++)
          {
            size_t idx = (size_t)i * N * N + j * N + k;
            double dre = FDM_Psi[idx][0] - pre;
            double dim = FDM_Psi[idx][1] - pim;
            double err = sqrt(dre * dre + dim * dim);
            if(err > max_phase_err)
              max_phase_err = err;
          }
    }

  double global_max_phase_err;
  MPI_Allreduce(&max_phase_err, &global_max_phase_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("FDM_TEST_DRIFT: N=%d, NTask=%d, kx0=%d, dt=%.4f\n", N, NTask, kx0, dt);
      printf("FDM_TEST_DRIFT: hbar_over_m_code=%.6e, predicted extra phase=%.6f rad\n", hbar_over_m, extra_phase);
      printf("FDM_TEST_DRIFT: norm before=%.6e, after=%.6e, rel_err=%.3e\n", norm_before_global, norm_after_global,
             norm_rel_err);
      printf("FDM_TEST_DRIFT: max phase-prediction error (abs, on unit circle) = %.3e\n", global_max_phase_err);

      if(norm_rel_err < 1e-10 && global_max_phase_err < 1e-8)
        printf("FDM_TEST_DRIFT: PASS\n");
      else
        printf("FDM_TEST_DRIFT: FAIL\n");
    }

  fdm_free();
  MPI_Finalize();
  return 0;
}
