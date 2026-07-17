/*!
 * \file        src/fdm/fdm_test_interpolation.c
 * \brief       Validates fdm_interpolate_to_stars() against a known
 *              analytic potential, Phi(x,y,z) = x^2+y^2+z^2 -- same
 *              test case as fdm_test_force.c, reused here since it's
 *              already proven sensitive to position-dependent indexing
 *              bugs (a linear test function would not be, see that
 *              file's own header for why).
 *
 *              IMPORTANT DIFFERENCE from fdm_test_force.c's tolerance:
 *              trilinear interpolation is NOT exact for a quadratic
 *              (only exact for functions linear in each coordinate
 *              separately) -- unlike the 4th-order finite-difference
 *              stencil, which WAS exact for a quadratic. A genuine,
 *              computable O(dx^2) interpolation error is expected
 *              here, not floating-point-precision agreement. Test
 *              stars are placed at the CENTER of a grid cell (dx=dy=dz
 *              in [0,1) locally) specifically so the expected error
 *              can be bounded analytically and checked against, rather
 *              than just eyeballing "is it small".
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

  MaxThreads = 1;

  All.MaxMemSize = 500;
  mymalloc_init();

  All.FDMGrid    = 32;
  All.FDMBoxSize = 20.0;
  All.FDMMass    = 1e-22;
  All.UnitLength_in_cm         = 3.085678e21;
  All.UnitMass_in_g            = 1.989e43;
  All.UnitVelocity_in_cm_per_s = 1e5;

  int N = All.FDMGrid;
  double dx_cell = All.FDMBoxSize / N;

  fdm_allocate();

  /* Fill FDM_Potential with the known quadratic (NOT centered this
   * time -- fdm_interpolate_to_stars() expects positions in [0, FDMBoxSize),
   * matching the raw grid convention, same as the actual particle
   * coupling will use; centering was specific to fdm_test_force.c's own
   * standalone setup). */
  for(int i = 0; i < FDM_plan.nslab_x; i++)
    {
      int gx = FDM_plan.slabstart_x + i;
      double x = gx * dx_cell;
      for(int j = 0; j < N; j++)
        {
          double y = j * dx_cell;
          for(int k = 0; k < N; k++)
            {
              double z = k * dx_cell;
              size_t idx = (size_t)i * N * N + (size_t)j * N + k;
              FDM_Potential[idx] = x * x + y * y + z * z;
            }
        }
    }

  fdm_compute_force();

  /* Minimal particle setup -- deliberately NOT calling allocate_memory(),
   * which pulls in timebins/other Arepo infrastructure requiring further
   * initialization this standalone test has no need for. Allocate
   * directly just what fdm_interpolate_to_stars() actually touches:
   * P[] and FDM_StarResult, both sized to a small MaxPart. */
  All.MaxPart = 100;
  P = (struct particle_data *)mymalloc("P", All.MaxPart * sizeof(struct particle_data));
  FDM_StarResult = (fdm_star_result *)mymalloc("FDM_StarResult", All.MaxPart * sizeof(fdm_star_result));

  /* Place test stars at cell-center offsets (dx=dy=dz=0.5 locally) at a
   * few different global positions -- deliberately NOT at a grid point
   * itself (which would trivially match exactly, telling us nothing
   * about the interpolation's genuine behavior between points). */
  int n_test = 5;
  double test_gx[5] = {5.5, 10.5, 15.5, 20.5, 25.5}; /* in grid-cell units, i.e. actual position = test_gx*dx_cell */

  NumPart = n_test;
  for(int t = 0; t < n_test; t++)
    {
      P[t].Type   = 4;
      P[t].Pos[0] = test_gx[t] * dx_cell;
      P[t].Pos[1] = test_gx[t] * dx_cell; /* same offset in y, z for simplicity -- still non-grid-aligned in all 3 dims */
      P[t].Pos[2] = test_gx[t] * dx_cell;
    }

  fdm_interpolate_to_stars();

  double max_rel_err_potential = 0.0, max_rel_err_force = 0.0;

  for(int t = 0; t < n_test; t++)
    {
      double x = P[t].Pos[0], y = P[t].Pos[1], z = P[t].Pos[2];
      double analytic_phi = x * x + y * y + z * z;
      double analytic_fx = -2.0 * x, analytic_fy = -2.0 * y, analytic_fz = -2.0 * z;

      double got_phi = FDM_StarResult[t].Potential;
      double got_fx = FDM_StarResult[t].ForceX, got_fy = FDM_StarResult[t].ForceY, got_fz = FDM_StarResult[t].ForceZ;

      /* Expected trilinear interpolation error for a quadratic, at a
       * cell-center offset: the standard error bound for linear
       * interpolation of f(x)=x^2 over an interval of width h, at the
       * midpoint, is exactly h^2/4 per dimension (second derivative is
       * 2, standard interpolation error formula -f''(xi)/8 * h^2 for
       * midpoint, i.e. dx_cell^2/4 here) -- computed as a genuine
       * analytic bound to check against, not just an arbitrary
       * tolerance. */
      double expected_phi_err_bound = 3 * (dx_cell * dx_cell / 4.0); /* one h^2/4 contribution per dimension */

      double err_phi = fabs(got_phi - analytic_phi);
      double err_fx = fabs(got_fx - analytic_fx);
      double err_fy = fabs(got_fy - analytic_fy);
      double err_fz = fabs(got_fz - analytic_fz);

      if(ThisTask == 0)
        printf("Star %d at (%.3f,%.3f,%.3f): phi_got=%.6f phi_analytic=%.6f err=%.6f (bound=%.6f) | "
               "F=(%.4f,%.4f,%.4f) analytic=(%.4f,%.4f,%.4f)\n",
               t, x, y, z, got_phi, analytic_phi, err_phi, expected_phi_err_bound, got_fx, got_fy, got_fz, analytic_fx, analytic_fy,
               analytic_fz);

      if(err_phi / expected_phi_err_bound > max_rel_err_potential)
        max_rel_err_potential = err_phi / expected_phi_err_bound;

      double force_mag = sqrt(analytic_fx * analytic_fx + analytic_fy * analytic_fy + analytic_fz * analytic_fz);
      double force_err_mag = sqrt(err_fx * err_fx + err_fy * err_fy + err_fz * err_fz);
      if(force_err_mag / force_mag > max_rel_err_force)
        max_rel_err_force = force_err_mag / force_mag;
    }

  if(ThisTask == 0)
    {
      printf("FDM_TEST_INTERP: max (err/expected_bound) for potential = %.4f (should be O(1), not >>1)\n", max_rel_err_potential);
      printf("FDM_TEST_INTERP: max relative force error = %.4e\n", max_rel_err_force);
      if(max_rel_err_potential < 3.0 && max_rel_err_force < 0.1)
        printf("FDM_TEST_INTERP: PASS\n");
      else
        printf("FDM_TEST_INTERP: FAIL\n");
    }

  myfree(FDM_StarResult);
  myfree(P);

  fdm_free();
  MPI_Finalize();
  return 0;
}
