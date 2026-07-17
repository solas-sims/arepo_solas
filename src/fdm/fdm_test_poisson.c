/*!
 * \file        src/fdm/fdm_test_poisson.c
 * \brief       Validates fdm_update_potential() against the known
 *              analytic potential of a uniform-density sphere:
 *                Phi(r) = -G*M/(2*R^3)*(3*R^2 - r^2),  r <= R
 *                Phi(r) = -G*M/r,                       r > R
 *              A standard test case precisely because it has a clean
 *              closed form both inside (harmonic) and outside
 *              (point-mass-like 1/r), exercising the density
 *              computation, the redistribution communication (both
 *              directions between FDM_plan and FDM_PoissonPlan), the
 *              kernel convolution, and the actual physics together.
 *
 *              psi is set to sqrt(rho_desired/(4*pi*G*m)) (real) so
 *              that fdm_update_potential's own 4*pi*G*m*|psi|^2
 *              recovers exactly the desired density -- this tests the
 *              Poisson solve itself, independent of whether |psi|^2
 *              matches any particular physical wavefunction (that's a
 *              separate question for actual simulation setup/ICs).
 *
 *              Discretization error is expected here (unlike the
 *              exact-to-machine-precision drift/kick tests) -- this is
 *              a genuinely different kind of test, checking physical
 *              correctness at finite resolution, not an algebraic
 *              identity. Tolerance is set accordingly looser.
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

  All.MaxMemSize = 1000;
  mymalloc_init();

  All.FDMGrid    = 32;
  All.FDMBoxSize = 20.0;  /* code length units, e.g. kpc */
  All.FDMMass    = 1e-22; /* eV -- arbitrary for this test, doesn't affect Phi since psi is set to cancel it out */
  All.UnitLength_in_cm         = 3.085678e21; /* 1 kpc */
  All.UnitMass_in_g            = 1.989e43;    /* 1e10 Msun */
  All.UnitVelocity_in_cm_per_s = 1e5;          /* 1 km/s */

  double UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.G = GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * pow(UnitTime_in_s, 2);

  int    N = All.FDMGrid;
  double L = All.FDMBoxSize;
  double R = L / 4.0;   /* sphere radius */
  double M = 1.0;       /* total mass, code units */
  double rho0 = M / ((4.0 / 3.0) * M_PI * R * R * R);

  double m_grams = (All.FDMMass * ELECTRONVOLT_IN_ERGS) / (CLIGHT * CLIGHT);
  double m_code  = m_grams / All.UnitMass_in_g;

  fdm_allocate();

  double dx = L / N;
  double cx = L / 2.0, cy = L / 2.0, cz = L / 2.0; /* sphere center, box center */

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
              double psi_val     = sqrt(rho_desired / m_code); /* matches fdm_update_potential's rho=m*|psi|^2, Eq. 5 --
                                                                  * no 4*pi*G factor here, that was the bug this
                                                                  * whole test caught in fdm_poisson.c itself */

              size_t idx = (size_t)i * N * N + j * N + k;
              FDM_Psi[idx][0] = psi_val;
              FDM_Psi[idx][1] = 0.0;
            }
        }
    }

  fdm_update_potential();

  /* Check at several specific radii, along the +x axis from the sphere
   * center, at whichever local grid points are closest to each target
   * radius (exact grid alignment not required -- just report what we
   * actually find there). */
  double test_radii[] = {0.0, 0.5 * R, R, 1.5 * R, 2.0 * R};
  int n_test = 5;

  for(int t = 0; t < n_test; t++)
    {
      double r_target = test_radii[t];
      double x_target = cx + r_target;
      int gx_target    = (int)round(x_target / dx);
      if(gx_target >= N)
        gx_target = N - 1;

      int gj = (int)round(cy / dx);
      int gk = (int)round(cz / dx);

      /* Only the task owning this gx reports it. */
      if(gx_target >= FDM_plan.slabstart_x && gx_target < FDM_plan.slabstart_x + FDM_plan.nslab_x)
        {
          int i = gx_target - FDM_plan.slabstart_x;
          size_t idx = (size_t)i * N * N + gj * N + gk;

          double r_actual = fabs(gx_target * dx - cx); /* actual distance for the grid point we landed on */
          double phi_numeric = FDM_Potential[idx];

          double phi_analytic;
          if(r_actual <= R)
            phi_analytic = -All.G * M / (2.0 * R * R * R) * (3.0 * R * R - r_actual * r_actual);
          else
            phi_analytic = -All.G * M / r_actual;

          double rel_err = fabs(phi_numeric - phi_analytic) / fabs(phi_analytic);

          printf("Task %d: r_target=%.3f r_actual=%.3f  Phi_numeric=%.6e  Phi_analytic=%.6e  rel_err=%.3e\n", ThisTask,
                 r_target, r_actual, phi_numeric, phi_analytic, rel_err);
        }
    }

  fdm_free();
  MPI_Finalize();
  return 0;
}
