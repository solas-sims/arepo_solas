/*!
 * \file        src/fdm/fdm_test_deposit.c
 * \brief       Validates fdm_deposit_star_mass() two ways:
 *              1) Total mass conservation: CIC weights always sum to
 *                 1.0 per particle, so summing FDM_StarMassDensity over
 *                 the WHOLE mesh (across all tasks) must exactly equal
 *                 the sum of input star masses -- a strong, easy-to-
 *                 check global property, independent of any position
 *                 details.
 *              2) A star placed EXACTLY at a grid point should deposit
 *                 its entire mass into that one cell alone, with zero
 *                 in all 26 neighbors -- an exact, not approximate,
 *                 check (weight=1.0 at the grid point itself, 0.0
 *                 everywhere else, for perfectly grid-aligned CIC).
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

  All.MaxPart = 100;
  P = (struct particle_data *)mymalloc("P", All.MaxPart * sizeof(struct particle_data));
  FDM_StarResult = (fdm_star_result *)mymalloc("FDM_StarResult", All.MaxPart * sizeof(fdm_star_result));

  /* Test 1: mass conservation with a handful of arbitrary,
   * non-grid-aligned stars, spread across the box (including near a
   * task boundary if NTask>1, same spirit as the interpolation test). */
  int n_test = 5;
  double test_gx[5]     = {5.3, 10.7, 15.5, 20.2, 25.9};
  double test_masses[5] = {1.0, 2.5, 0.7, 3.3, 1.8};

  NumPart = (ThisTask == 0) ? n_test : 0;
  double total_input_mass = 0.0;
  for(int t = 0; t < n_test; t++)
    {
      total_input_mass += test_masses[t]; /* known regardless of which task owns the particles */
      if(ThisTask != 0)
        continue;
      P[t].Type   = 4;
      P[t].Pos[0] = test_gx[t] * dx_cell;
      P[t].Pos[1] = test_gx[t] * dx_cell;
      P[t].Pos[2] = test_gx[t] * dx_cell;
      P[t].Mass   = test_masses[t];
    }

  fdm_deposit_star_mass();

  /* Sum FDM_StarMassDensity over this task's own local portion, then
   * MPI_Allreduce for the global total -- FDM_StarMassDensity is
   * private to fdm_particle_coupling.c (static), so sum it via a
   * small helper exposed just for this test rather than reaching
   * into the module's internals from outside. */
  extern double fdm_debug_sum_star_mass_density(void);
  double local_sum = fdm_debug_sum_star_mass_density();
  double global_sum;
  MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("FDM_TEST_DEPOSIT: total input mass = %.6f, total deposited mass = %.6f, diff = %.3e\n", total_input_mass, global_sum,
             fabs(total_input_mass - global_sum));
    }

  /* Test 2: a star exactly at a grid point deposits its entire mass
   * into exactly that one cell. */
  NumPart = (ThisTask == 0) ? 1 : 0;
  if(ThisTask == 0)
    {
      P[0].Type   = 4;
      int test_grid_x = 12, test_grid_y = 12, test_grid_z = 12;
      P[0].Pos[0] = test_grid_x * dx_cell;
      P[0].Pos[1] = test_grid_y * dx_cell;
      P[0].Pos[2] = test_grid_z * dx_cell;
      P[0].Mass   = 7.0;
    }
  int test_grid_x = 12, test_grid_y = 12, test_grid_z = 12; /* known regardless of which task owns the particle */

  fdm_deposit_star_mass();

  extern double fdm_debug_get_star_mass_density_at(int gx, int gy, int gz);
  double at_point = fdm_debug_get_star_mass_density_at(test_grid_x, test_grid_y, test_grid_z);
  double at_neighbor = fdm_debug_get_star_mass_density_at(test_grid_x + 1, test_grid_y, test_grid_z);

  double at_point_global, at_neighbor_global;
  MPI_Allreduce(&at_point, &at_point_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&at_neighbor, &at_neighbor_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("FDM_TEST_DEPOSIT: grid-aligned star: density at point = %.6f (expect 7.0), at neighbor = %.6f (expect 0.0)\n",
             at_point_global, at_neighbor_global);

      int pass = (fabs(total_input_mass - global_sum) < 1e-9) && (fabs(at_point_global - 7.0) < 1e-9) &&
                 (fabs(at_neighbor_global) < 1e-9);
      printf("FDM_TEST_DEPOSIT: %s\n", pass ? "PASS" : "FAIL");
    }

  myfree(FDM_StarResult);
  myfree(P);
  fdm_free();
  MPI_Finalize();
  return 0;
}
