/*!
 * \file        src/fdm/fdm_test_kick.c
 * \brief       Validates fdm_kick(): (1) |psi| is preserved pointwise,
 *              at every single grid point independently (stronger than
 *              a total-norm check); (2) for a spatially-VARYING known
 *              test potential (not just a uniform one -- chosen
 *              deliberately so different grid points get genuinely
 *              different expected phases, which would catch an
 *              indexing mistake a uniform potential could not), the
 *              resulting phase shift matches the analytic prediction
 *              -(m/hbar)*dt*Phi(x) exactly at every point.
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
  All.FDMBoxSize = 10.0;
  All.FDMMass    = 1e-22;
  All.UnitLength_in_cm         = 3.085678e21;
  All.UnitMass_in_g            = 1.989e43;
  All.UnitVelocity_in_cm_per_s = 1e5;

  int    N  = All.FDMGrid;
  double dt = 0.01;

  fdm_allocate();

  size_t local_size = (size_t)FDM_plan.nslab_x * N * N;

  /* Arbitrary, non-trivial complex test field, and a spatially-varying
   * test potential Phi(x,y,z) = Phi0 * sin(2*pi*gx/N) * cos(2*pi*(gy+gz)/N),
   * genuinely different at (almost) every grid point. */
  double Phi0 = 3.7; /* arbitrary code-units amplitude */
  double *original_mag = (double *)mymalloc("original_mag", local_size * sizeof(double));

  for(int i = 0; i < FDM_plan.nslab_x; i++)
    {
      int gx = FDM_plan.slabstart_x + i;
      for(int j = 0; j < N; j++)
        for(int k = 0; k < N; k++)
          {
            size_t idx = (size_t)i * N * N + j * N + k;

            double re = 0.3 + 0.1 * gx + 0.05 * j - 0.02 * k; /* arbitrary, non-trivial */
            double im = 0.2 - 0.07 * gx + 0.11 * j + 0.03 * k;
            FDM_Psi[idx][0] = re;
            FDM_Psi[idx][1] = im;
            original_mag[idx] = sqrt(re * re + im * im);

            FDM_Potential[idx] = Phi0 * sin(2 * M_PI * gx / N) * cos(2 * M_PI * (j + k) / N);
          }
    }

  fdm_kick(dt);

  double m_over_hbar = 1.0 / fdm_hbar_over_m_code();

  double max_mag_err = 0.0, max_phase_err = 0.0;
  for(int i = 0; i < FDM_plan.nslab_x; i++)
    {
      int gx = FDM_plan.slabstart_x + i;
      for(int j = 0; j < N; j++)
        for(int k = 0; k < N; k++)
          {
            size_t idx = (size_t)i * N * N + j * N + k;

            double new_mag = sqrt(FDM_Psi[idx][0] * FDM_Psi[idx][0] + FDM_Psi[idx][1] * FDM_Psi[idx][1]);
            double mag_err = fabs(new_mag - original_mag[idx]);
            if(mag_err > max_mag_err)
              max_mag_err = mag_err;

            /* predicted new value: original * exp(-i*m_over_hbar*dt*Phi) */
            double phi_here = Phi0 * sin(2 * M_PI * gx / N) * cos(2 * M_PI * (j + k) / N);
            double phase    = -m_over_hbar * dt * phi_here;
            double c = cos(phase), s = sin(phase);

            /* reconstruct original re/im from the same deterministic formula used above */
            double ore = 0.3 + 0.1 * gx + 0.05 * j - 0.02 * k;
            double oim = 0.2 - 0.07 * gx + 0.11 * j + 0.03 * k;

            double pre = ore * c - oim * s;
            double pim = ore * s + oim * c;

            double dre = FDM_Psi[idx][0] - pre;
            double dim = FDM_Psi[idx][1] - pim;
            double err = sqrt(dre * dre + dim * dim);
            if(err > max_phase_err)
              max_phase_err = err;
          }
    }

  double global_max_mag_err, global_max_phase_err;
  MPI_Allreduce(&max_mag_err, &global_max_mag_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&max_phase_err, &global_max_phase_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("FDM_TEST_KICK: N=%d, NTask=%d, dt=%.4f, Phi0=%.2f\n", N, NTask, dt, Phi0);
      printf("FDM_TEST_KICK: max pointwise |psi| error = %.3e\n", global_max_mag_err);
      printf("FDM_TEST_KICK: max pointwise phase-prediction error = %.3e\n", global_max_phase_err);

      if(global_max_mag_err < 1e-10 && global_max_phase_err < 1e-8)
        printf("FDM_TEST_KICK: PASS\n");
      else
        printf("FDM_TEST_KICK: FAIL\n");
    }

  myfree(original_mag);
  fdm_free();
  MPI_Finalize();
  return 0;
}
