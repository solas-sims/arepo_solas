/*!
 * \file        src/fdm/fdm_test_force.c
 * \brief       Validates fdm_compute_force() against a known analytic
 *              potential, Phi(x,y,z) = x^2+y^2+z^2, whose gradient is
 *              exactly (2x,2y,2z) -- Force = -grad(Phi) = (-2x,-2y,-2z).
 *
 *              Deliberately quadratic, not linear -- a linear test
 *              potential's derivative is constant everywhere, so it
 *              would give a "correct-looking" answer even with a
 *              genuine position-dependent indexing bug (this
 *              project's own history: an earlier linear test passed
 *              despite a real bug in the X-direction transpose layout,
 *              precisely because a constant derivative is insensitive
 *              to exactly where you sample it). The quadratic's
 *              derivative genuinely varies with position, so this test
 *              is actually sensitive to indexing correctness.
 *
 *              The stencil's own error term involves the 5th
 *              derivative (verified via Taylor expansion in
 *              fdm_gradient.c's header), which is nonzero for x^3 and
 *              higher but IS exactly zero for a quadratic's own
 *              higher derivatives in the relevant sense -- a correctly
 *              implemented 4-point stencil should reproduce this
 *              gradient to floating-point precision in the interior.
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

  All.FDMGrid    = 32;
  All.FDMBoxSize = 20.0;
  All.FDMMass    = 1e-22;
  All.UnitLength_in_cm         = 3.085678e21;
  All.UnitMass_in_g            = 1.989e43;
  All.UnitVelocity_in_cm_per_s = 1e5;

  int N = All.FDMGrid;
  double dx = All.FDMBoxSize / N;

  fdm_allocate();

  double cx = All.FDMBoxSize / 2.0, cy = cx, cz = cx;
  for(int i = 0; i < FDM_plan.nslab_x; i++)
    {
      int gx = FDM_plan.slabstart_x + i;
      double x = gx * dx - cx;
      for(int j = 0; j < N; j++)
        {
          double y = j * dx - cy;
          for(int k = 0; k < N; k++)
            {
              double z = k * dx - cz;
              size_t idx = (size_t)i * N * N + (size_t)j * N + k;
              FDM_Potential[idx] = x * x + y * y + z * z;
            }
        }
    }

  fdm_compute_force();

  double max_interior_err = 0.0, max_interior_val = 0.0;
  double max_edge_err = 0.0;

  for(int i = 0; i < FDM_plan.nslab_x; i++)
    {
      int gx = FDM_plan.slabstart_x + i;
      double x = gx * dx - cx;
      int x_edge = (gx == 0 || gx == N - 1);
      for(int j = 0; j < N; j++)
        {
          double y = j * dx - cy;
          int y_edge = (j == 0 || j == N - 1);
          for(int k = 0; k < N; k++)
            {
              double z = k * dx - cz;
              int z_edge = (k == 0 || k == N - 1);

              size_t idx = (size_t)i * N * N + (size_t)j * N + k;

              double fx_analytic = -2.0 * x;
              double fy_analytic = -2.0 * y;
              double fz_analytic = -2.0 * z;

              double err_x = fabs(FDM_ForceX[idx] - fx_analytic);
              double err_y = fabs(FDM_ForceY[idx] - fy_analytic);
              double err_z = fabs(FDM_ForceZ[idx] - fz_analytic);

              if(x_edge)
                { if(err_x > max_edge_err) max_edge_err = err_x; }
              else
                {
                  if(err_x > max_interior_err) max_interior_err = err_x;
                  if(fabs(fx_analytic) > max_interior_val) max_interior_val = fabs(fx_analytic);
                }

              if(y_edge)
                { if(err_y > max_edge_err) max_edge_err = err_y; }
              else
                {
                  if(err_y > max_interior_err) max_interior_err = err_y;
                  if(fabs(fy_analytic) > max_interior_val) max_interior_val = fabs(fy_analytic);
                }

              if(z_edge)
                { if(err_z > max_edge_err) max_edge_err = err_z; }
              else
                {
                  if(err_z > max_interior_err) max_interior_err = err_z;
                  if(fabs(fz_analytic) > max_interior_val) max_interior_val = fabs(fz_analytic);
                }
            }
        }
    }

  double global_max_interior_err, global_max_interior_val, global_max_edge_err;
  MPI_Allreduce(&max_interior_err, &global_max_interior_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&max_interior_val, &global_max_interior_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&max_edge_err, &global_max_edge_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("FDM_TEST_FORCE: N=%d, NTask=%d\n", N, NTask);
      printf("FDM_TEST_FORCE: interior (away from absolute edge) max abs error = %.3e (rel to max|F|=%.3e: %.3e)\n",
             global_max_interior_err, global_max_interior_val, global_max_interior_err / global_max_interior_val);
      printf("FDM_TEST_FORCE: absolute-edge (one-sided fallback) max abs error = %.3e "
             "(expected O(dx)=%.3e, NOT a bug if consistent with this)\n",
             global_max_edge_err, dx);

      if(global_max_interior_err / global_max_interior_val < 1e-10)
        printf("FDM_TEST_FORCE: PASS\n");
      else
        printf("FDM_TEST_FORCE: FAIL\n");
    }

  fdm_free();
  MPI_Finalize();
  return 0;
}
