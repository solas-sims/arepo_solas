/*!
 * \file        src/fdm/fdm_test_fft.c
 * \brief       Standalone validation test for the FDM wavefunction
 *              field's FFT machinery -- NOT part of the main Arepo
 *              binary. Compiled/linked separately (its own main()),
 *              reusing the already-built object files from the main
 *              build (fdm_field.o, mymalloc.o, pm_mpi_fft.o, etc.)
 *              minus main.o.
 *
 *              Two independent checks, deliberately chosen to avoid
 *              needing to understand or verify the internal slab-
 *              transpose reordering that happens during the
 *              distributed FFT (a real source of indexing risk that
 *              is easy to get subtly wrong when reasoning about it
 *              directly):
 *
 *              1. Round-trip recovery: forward FFT then backward FFT
 *                 should recover the original field (up to the
 *                 standard unnormalized-FFTW factor of N^3) to
 *                 floating-point precision. This alone would not catch
 *                 e.g. a wrong axis being transformed or a scale
 *                 error that happens to cancel between forward and
 *                 backward -- hence check 2.
 *
 *              2. Parseval's theorem: total power in real space
 *                 (sum|psi_x|^2) and total power in k-space
 *                 (sum|psi_k|^2) must be related by exactly a factor
 *                 of N^3 (the unnormalized FFTW convention), summed
 *                 GLOBALLY across all tasks. This is a genuine,
 *                 independent invariant -- it holds regardless of how
 *                 array elements are internally reordered/transposed
 *                 between tasks, since it only depends on every
 *                 element being counted exactly once, not on knowing
 *                 which physical mode ended up at which index.
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"
#include "fdm.h"

/* Minimal stub: endrun() is defined in main.c, which we deliberately
 * exclude from this standalone test's link (main.c has its own
 * main()). Only referenced from an error path in init.c that this
 * test's minimal init sequence doesn't reach. */
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

  All.MaxMemSize = 500; /* MB, plenty for a small test mesh */
  mymalloc_init();

  All.FDMGrid = 32; /* small test mesh -- fast, but non-trivial for multi-task slab decomposition */
  int N        = All.FDMGrid;

  fdm_allocate();

  size_t local_size = (size_t)FDM_plan.nslab_x * N * N;

  /* Save a copy of the original field for later comparison -- the
   * forward transform will overwrite FDM_Psi in place (workspace is
   * scratch, but the final forward result also lands back in the
   * original FDM_Psi pointer, per my_slab_based_fft_c2c's own
   * documented behaviour, matching how gravity's PM solver uses the
   * same pattern). */
  fft_complex *original = (fft_complex *)mymalloc("original", local_size * sizeof(fft_complex));

  /* Fill with a smooth, non-trivial, deterministic complex pattern
   * with genuine frequency content (not just a DC/single-mode term),
   * using the GLOBAL x-index (slabstart_x + local i) so the pattern
   * is well-defined regardless of task count. */
  double real_space_power_local = 0.0;
  for(int i = 0; i < FDM_plan.nslab_x; i++)
    {
      int gx = FDM_plan.slabstart_x + i;
      for(int j = 0; j < N; j++)
        {
          for(int k = 0; k < N; k++)
            {
              size_t idx = (size_t)i * N * N + j * N + k;
              double re  = sin(2 * M_PI * gx / N);
              double im  = cos(2 * M_PI * (j + k) / N);
              FDM_Psi[idx][0] = re;
              FDM_Psi[idx][1] = im;
              original[idx][0] = re;
              original[idx][1] = im;
              real_space_power_local += re * re + im * im;
            }
        }
    }

  double real_space_power_global;
  MPI_Allreduce(&real_space_power_local, &real_space_power_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  /* Forward transform */
  my_slab_based_fft_c2c(&FDM_plan, FDM_Psi, FDM_PsiWorkspace, 1);

  double k_space_power_local = 0.0;
  for(size_t idx = 0; idx < local_size; idx++)
    k_space_power_local += FDM_Psi[idx][0] * FDM_Psi[idx][0] + FDM_Psi[idx][1] * FDM_Psi[idx][1];

  double k_space_power_global;
  MPI_Allreduce(&k_space_power_local, &k_space_power_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  double N3 = (double)N * N * N;
  double parseval_ratio    = k_space_power_global / real_space_power_global;
  double parseval_rel_err  = fabs(parseval_ratio - N3) / N3;

  /* Backward transform */
  my_slab_based_fft_c2c(&FDM_plan, FDM_Psi, FDM_PsiWorkspace, -1);

  double max_abs_err = 0.0, max_val = 0.0;
  for(size_t idx = 0; idx < local_size; idx++)
    {
      double re = FDM_Psi[idx][0] / N3;
      double im = FDM_Psi[idx][1] / N3;
      double dre = re - original[idx][0];
      double dim = im - original[idx][1];
      double err = sqrt(dre * dre + dim * dim);
      double val = sqrt(original[idx][0] * original[idx][0] + original[idx][1] * original[idx][1]);
      if(err > max_abs_err)
        max_abs_err = err;
      if(val > max_val)
        max_val = val;
    }

  double global_max_abs_err, global_max_val;
  MPI_Allreduce(&max_abs_err, &global_max_abs_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&max_val, &global_max_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("FDM_TEST_FFT: N=%d, NTask=%d\n", N, NTask);
      printf("FDM_TEST_FFT: real-space power (global) = %.6e\n", real_space_power_global);
      printf("FDM_TEST_FFT: k-space power (global)     = %.6e\n", k_space_power_global);
      printf("FDM_TEST_FFT: Parseval check: k-space/real-space = %.6e, expected N^3 = %.6e, rel_err = %.3e\n",
             parseval_ratio, N3, parseval_rel_err);
      printf("FDM_TEST_FFT: round-trip max abs error = %.3e (max|original|=%.3e, rel = %.3e)\n", global_max_abs_err,
             global_max_val, global_max_abs_err / global_max_val);

      if(parseval_rel_err < 1e-10 && (global_max_abs_err / global_max_val) < 1e-10)
        printf("FDM_TEST_FFT: PASS\n");
      else
        printf("FDM_TEST_FFT: FAIL\n");
    }

  myfree(original);
  fdm_free();

  MPI_Finalize();
  return 0;
}
