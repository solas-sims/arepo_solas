/*!
 * \file        src/fdm/fdm_test_single_mode.c
 * \brief       Validates the EXACT post-forward-transform array layout
 *              (which array index corresponds to which (kx,ky,kz)
 *              triple) by inputting a single known Fourier mode and
 *              checking the transformed result is a single spike at
 *              the predicted location. This is specifically needed
 *              before implementing the drift operator, which must
 *              multiply each k-space element by a phase depending on
 *              k^2 -- unlike the round-trip/Parseval test, this cannot
 *              be done in an index-order-agnostic way.
 *
 *              Derived layout (from tracing my_slab_transpose and the
 *              final x-direction FFTW stride convention directly, not
 *              assumed): after the full forward c2c transform,
 *              index = i_local_y * N * N + kx_global * N + kz_global
 *              where i_local_y is LOCAL to this task (add
 *              FDM_plan.slabstart_y for the global y-index), and
 *              kx_global, kz_global are GLOBAL grid indices (0..N-1),
 *              since x and z are both fully local after the transpose
 *              (only y remains distributed).
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

  All.FDMGrid = 16; /* small, so we can afford to print/inspect every significant element */
  int N        = All.FDMGrid;

  fdm_allocate();

  size_t local_size = (size_t)FDM_plan.nslab_x * N * N;

  /* Pure plane wave along X only: psi(x,y,z) = exp(i*2*pi*kx0*x/N) */
  int kx0 = 3;
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
          }
    }

  my_slab_based_fft_c2c(&FDM_plan, FDM_Psi, FDM_PsiWorkspace, 1);

  /* Predicted spike location: kx_global=kx0, ky_global=0, kz_global=0.
   * Only the task owning ky_global=0 in its local y-slab should see it. */
  double threshold = 1e-6 * N * N * N; /* spike magnitude should be ~N^3 (single mode, Parseval) */

  printf("Task %d: scanning for spikes (N=%d, kx0=%d, nslab_y_local=%d, slabstart_y=%d)\n", ThisTask, N, kx0,
         FDM_plan.nslab_y, FDM_plan.slabstart_y);

  int nfound = 0;
  for(int i = 0; i < FDM_plan.nslab_y; i++)
    {
      int gy = FDM_plan.slabstart_y + i;
      for(int kx = 0; kx < N; kx++)
        for(int kz = 0; kz < N; kz++)
          {
            size_t idx = (size_t)i * N * N + kx * N + kz;
            double mag = sqrt(FDM_Psi[idx][0] * FDM_Psi[idx][0] + FDM_Psi[idx][1] * FDM_Psi[idx][1]);
            if(mag > threshold)
              {
                printf("Task %d: SPIKE at local_idx=(i=%d,kx=%d,kz=%d) -> global (kx=%d,ky=%d,kz=%d), mag=%.4e, "
                       "phase=%.4f\n",
                       ThisTask, i, kx, kz, kx, gy, kz, mag, atan2(FDM_Psi[idx][1], FDM_Psi[idx][0]));
                nfound++;
              }
          }
    }

  if(nfound == 0)
    printf("Task %d: no spike found locally (expected if this task doesn't own global ky=0)\n", ThisTask);

  int total_found;
  MPI_Allreduce(&nfound, &total_found, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(ThisTask == 0)
    printf("TOTAL spikes found across all tasks: %d (expect exactly 1, at kx=%d, ky=0, kz=0)\n", total_found, kx0);

  fdm_free();
  MPI_Finalize();
  return 0;
}
