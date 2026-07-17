/*!
 * \file        src/fdm/fdm_test_step.c
 * \brief       Validates the full fdm_step() sequence (kick-drift-kick,
 *              Eqs. 17-20) over multiple steps. Total norm
 *              (sum|psi|^2, globally) is a genuine physical invariant
 *              here: both kick (pointwise phase rotation) and drift
 *              (phase rotation in k-space) are individually exactly
 *              norm-preserving, so the composed sequence should be
 *              too, REGARDLESS of the specific potential values -- any
 *              drift in total norm across repeated steps would
 *              specifically indicate a bug in how the pieces are
 *              sequenced/composed, not just ordinary physics
 *              discretization error (which the individual drift/kick/
 *              Poisson tests already characterized separately).
 *
 *              Uses the same uniform-sphere-consistent psi as
 *              fdm_test_poisson.c for a physically-motivated (not
 *              arbitrary) initial condition, reusing validated setup.
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

static double total_norm(int N)
{
  size_t local_size = (size_t)FDM_plan.nslab_x * N * N;
  double local = 0.0;
  for(size_t idx = 0; idx < local_size; idx++)
    local += FDM_Psi[idx][0] * FDM_Psi[idx][0] + FDM_Psi[idx][1] * FDM_Psi[idx][1];

  double global;
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return global;
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  All.MaxMemSize = 1000;
  mymalloc_init();

  All.FDMGrid    = 32;
  All.FDMBoxSize = 20.0;
  All.FDMMass    = 1e-22;
  All.UnitLength_in_cm         = 3.085678e21;
  All.UnitMass_in_g            = 1.989e43;
  All.UnitVelocity_in_cm_per_s = 1e5;

  double UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.G = GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * pow(UnitTime_in_s, 2);

  int    N = All.FDMGrid;
  double L = All.FDMBoxSize;
  double R = L / 4.0;
  double M = 1.0;
  double rho0 = M / ((4.0 / 3.0) * M_PI * R * R * R);

  double m_grams = (All.FDMMass * ELECTRONVOLT_IN_ERGS) / (CLIGHT * CLIGHT);
  double m_code  = m_grams / All.UnitMass_in_g;

  fdm_allocate();

  double dx = L / N;
  double cx = L / 2.0, cy = L / 2.0, cz = L / 2.0;

  for(int i = 0; i < FDM_plan.nslab_x; i++)
    {
      int gx = FDM_plan.slabstart_x + i;
      double x = gx * dx;
      for(int j = 0; j < N; j++)
        {
          double y = j * dx;
          for(int k = 0; k < N; k++)
            {
              double z = k * dx;
              double r = sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy) + (z - cz) * (z - cz));
              double rho_desired = (r <= R) ? rho0 : 0.0;
              double psi_val     = sqrt(rho_desired / m_code);

              size_t idx = (size_t)i * N * N + j * N + k;
              FDM_Psi[idx][0] = psi_val;
              FDM_Psi[idx][1] = 0.0;
            }
        }
    }

  double norm_initial = total_norm(N);

  /* Establish the initial potential -- precondition for fdm_step()'s
   * first call, per its own documented requirement. */
  fdm_update_potential();

  double dt = fdm_get_timestep();

  if(ThisTask == 0)
    printf("FDM_TEST_STEP: N=%d, NTask=%d, dt=%.6e, norm_initial=%.6e\n", N, NTask, dt, norm_initial);

  if(!(dt > 0.0) || !isfinite(dt))
    {
      if(ThisTask == 0)
        printf("FDM_TEST_STEP: FAIL -- timestep criterion produced non-positive or non-finite dt\n");
      fdm_free();
      MPI_Finalize();
      return 1;
    }

  int n_steps = 10;
  double max_norm_rel_drift = 0.0;

  for(int s = 0; s < n_steps; s++)
    {
      fdm_step(dt);
      double norm_now = total_norm(N);
      double rel_drift = fabs(norm_now - norm_initial) / norm_initial;
      if(rel_drift > max_norm_rel_drift)
        max_norm_rel_drift = rel_drift;

      if(ThisTask == 0)
        printf("FDM_TEST_STEP: after step %d, norm=%.6e, rel_drift_from_initial=%.3e\n", s + 1, norm_now, rel_drift);
    }

  if(ThisTask == 0)
    {
      printf("FDM_TEST_STEP: max norm relative drift over %d steps = %.3e\n", n_steps, max_norm_rel_drift);
      if(max_norm_rel_drift < 1e-8)
        printf("FDM_TEST_STEP: PASS\n");
      else
        printf("FDM_TEST_STEP: FAIL\n");
    }

  fdm_free();
  MPI_Finalize();
  return 0;
}
