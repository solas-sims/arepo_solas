/*!
 * \file        src/fdm/fdm_io.c
 * \brief       Basic I/O for the FDM wavefunction field: gather-and-
 *              write / read-and-scatter, using plain MPI_Gatherv/
 *              MPI_Scatterv (task 0 assembles or distributes the full
 *              field) rather than parallel HDF5 -- deliberately simple
 *              for "basic version first", not a scalability
 *              optimization deferred by oversight. At the mesh sizes
 *              this project has sized as tractable so far (N up to a
 *              few hundred), gathering the full field to one task is a
 *              bounded, modest cost (N^3 * 16 bytes for double-
 *              precision real+imag; e.g. ~268 MB at N=256), not a
 *              genuine bottleneck yet.
 *
 *              This file defines BOTH the periodic-output format
 *              during a run AND, deliberately the same format, the
 *              interface an external initial-condition generator must
 *              produce to be read back in via fdm_read_field() --
 *              designing these symmetrically avoids needing two
 *              different formats/pieces of code for what is
 *              structurally the same problem (getting a global psi
 *              array into or out of FDM_plan's distributed layout).
 *
 *              File format (HDF5): attributes "N" (int), "L" (double,
 *              box size), "FDMMass" (double, mc^2 in eV -- checked on
 *              read as a sanity check, not silently ignored if it
 *              doesn't match the current run's All.FDMMass); datasets
 *              "PsiReal" and "PsiImag", each N*N*N doubles in plain
 *              global (x,y,z) row-major order (x slowest, z fastest) --
 *              NOT FDM_plan's distributed layout, which is an
 *              implementation detail of this specific run's task
 *              count, not something an external IC generator should
 *              need to know about.
 *
 *              Always double precision in the FILE regardless of
 *              whether this build uses DOUBLEPRECISION_FFTW -- keeps
 *              the file format stable for an external generator
 *              (numpy's natural type) independent of this specific
 *              build's FFTW precision choice. Converted to/from
 *              fft_real (which may be float) during gather/scatter.
 */

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"
#include "fdm.h"

#ifdef FDM

#ifdef DOUBLEPRECISION_FFTW
#define FDM_MPI_REAL_TYPE MPI_DOUBLE
#else /* #ifdef DOUBLEPRECISION_FFTW */
#define FDM_MPI_REAL_TYPE MPI_FLOAT
#endif /* #ifdef DOUBLEPRECISION_FFTW #else */

/*! \brief Gathers FDM_Psi (distributed per FDM_plan) to task 0 and
 *  writes it to an HDF5 file in plain global (x,y,z) order.
 *
 *  Uses FDM_plan's own slabs_x_per_task[]/first_slab_x_of_task[]
 *  bookkeeping (already populated for exactly this kind of purpose) to
 *  build the MPI_Gatherv counts/displacements -- a standard collective
 *  gather, genuinely simpler than fdm_poisson.c's row redistribution,
 *  since there is only one target layout (task 0's fully global
 *  buffer) rather than two independently-decomposed plans to
 *  reconcile.
 */
void fdm_write_field(const char *fname)
{
  int N = All.FDMGrid;

  fft_complex *global_psi = NULL;
  int         *recvcounts = NULL, *displs = NULL;

  if(ThisTask == 0)
    {
      global_psi = (fft_complex *)mymalloc("fdm_io_global_psi", (size_t)N * N * N * sizeof(fft_complex));
      recvcounts = (int *)mymalloc("fdm_io_recvcounts", NTask * sizeof(int));
      displs     = (int *)mymalloc("fdm_io_displs", NTask * sizeof(int));

      for(int t = 0; t < NTask; t++)
        {
          /* *2 for real+imag as separate reals -- fft_complex is a
           * 2-element array (fftw_complex/fftwf_complex), safe to treat
           * as 2*count reals of the matching MPI type for gather/
           * scatter purposes. */
          recvcounts[t] = FDM_plan.slabs_x_per_task[t] * N * N * 2;
          displs[t]     = FDM_plan.first_slab_x_of_task[t] * N * N * 2;
        }
    }

  int sendcount = FDM_plan.nslab_x * N * N * 2;
  MPI_Gatherv(FDM_Psi, sendcount, FDM_MPI_REAL_TYPE, global_psi, recvcounts, displs, FDM_MPI_REAL_TYPE, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      size_t total = (size_t)N * N * N;
      double *psi_real = (double *)mymalloc("fdm_io_psi_real", total * sizeof(double));
      double *psi_imag = (double *)mymalloc("fdm_io_psi_imag", total * sizeof(double));

      for(size_t idx = 0; idx < total; idx++)
        {
          psi_real[idx] = (double)global_psi[idx][0];
          psi_imag[idx] = (double)global_psi[idx][1];
        }

      hid_t file_id = my_H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

      hsize_t dims[3] = {(hsize_t)N, (hsize_t)N, (hsize_t)N};
      hid_t   space_id = my_H5Screate_simple(3, dims, NULL);

      hid_t dset_real = my_H5Dcreate(file_id, "PsiReal", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
      my_H5Dwrite(dset_real, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, psi_real, "PsiReal");
      H5Dclose(dset_real);

      hid_t dset_imag = my_H5Dcreate(file_id, "PsiImag", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
      my_H5Dwrite(dset_imag, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, psi_imag, "PsiImag");
      H5Dclose(dset_imag);

      H5Sclose(space_id);

      hid_t scalar_space = my_H5Screate(H5S_SCALAR);
      hid_t attr_N       = my_H5Acreate(file_id, "N", H5T_NATIVE_INT, scalar_space, H5P_DEFAULT);
      my_H5Awrite(attr_N, H5T_NATIVE_INT, &N, "N");
      H5Aclose(attr_N);

      hid_t attr_L = my_H5Acreate(file_id, "L", H5T_NATIVE_DOUBLE, scalar_space, H5P_DEFAULT);
      my_H5Awrite(attr_L, H5T_NATIVE_DOUBLE, &All.FDMBoxSize, "L");
      H5Aclose(attr_L);

      hid_t attr_mass = my_H5Acreate(file_id, "FDMMass", H5T_NATIVE_DOUBLE, scalar_space, H5P_DEFAULT);
      my_H5Awrite(attr_mass, H5T_NATIVE_DOUBLE, &All.FDMMass, "FDMMass");
      H5Aclose(attr_mass);

      H5Sclose(scalar_space);
      H5Fclose(file_id);

      myfree(psi_imag);
      myfree(psi_real);
      myfree(displs);
      myfree(recvcounts);
      myfree(global_psi);

      mpi_printf("FDM: wrote field to %s (N=%d)\n", fname, N);
    }
}

/*! \brief Reads an HDF5 file in the format fdm_write_field() produces
 *  (or that an external IC generator matching that format produces)
 *  and scatters it into FDM_Psi's distributed layout. Reverse of
 *  fdm_write_field().
 *
 *  Checks the file's "N" attribute against All.FDMGrid and its
 *  "FDMMass" attribute against All.FDMMass -- terminates with a clear
 *  message on mismatch rather than silently proceeding with
 *  inconsistent assumptions (an external generator producing ICs for
 *  the wrong mesh size or mass is a real, easy mistake to make, not a
 *  hypothetical).
 */
void fdm_read_field(const char *fname)
{
  int N = All.FDMGrid;

  fft_complex *global_psi = NULL;
  int         *sendcounts = NULL, *displs = NULL;
  double      *psi_real = NULL, *psi_imag = NULL;

  if(ThisTask == 0)
    {
      hid_t file_id = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

      hid_t scalar_space = my_H5Screate(H5S_SCALAR);

      int file_N;
      hid_t attr_N = my_H5Aopen_name(file_id, "N");
      my_H5Aread(attr_N, H5T_NATIVE_INT, &file_N, "N", 1);
      H5Aclose(attr_N);

      if(file_N != N)
        terminate("FDM: IC file %s has N=%d, but All.FDMGrid=%d -- regenerate the IC for the current mesh size.\n", fname,
                  file_N, N);

      double file_mass;
      hid_t  attr_mass = my_H5Aopen_name(file_id, "FDMMass");
      my_H5Aread(attr_mass, H5T_NATIVE_DOUBLE, &file_mass, "FDMMass", 1);
      H5Aclose(attr_mass);

      if(fabs(file_mass - All.FDMMass) > 1e-10 * All.FDMMass)
        terminate("FDM: IC file %s has FDMMass=%.6e eV, but All.FDMMass=%.6e eV -- these must match.\n", fname, file_mass,
                  All.FDMMass);

      double file_L;
      hid_t  attr_L = my_H5Aopen_name(file_id, "L");
      my_H5Aread(attr_L, H5T_NATIVE_DOUBLE, &file_L, "L", 1);
      H5Aclose(attr_L);

      if(fabs(file_L - All.FDMBoxSize) > 1e-10 * All.FDMBoxSize)
        terminate("FDM: IC file %s has L=%.6e, but All.FDMBoxSize=%.6e -- these must match.\n", fname, file_L,
                  All.FDMBoxSize);

      H5Sclose(scalar_space);

      size_t total = (size_t)N * N * N;
      psi_real = (double *)mymalloc("fdm_io_psi_real", total * sizeof(double));
      psi_imag = (double *)mymalloc("fdm_io_psi_imag", total * sizeof(double));

      hid_t dset_real = my_H5Dopen(file_id, "PsiReal");
      my_H5Dread(dset_real, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, psi_real, "PsiReal");
      H5Dclose(dset_real);

      hid_t dset_imag = my_H5Dopen(file_id, "PsiImag");
      my_H5Dread(dset_imag, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, psi_imag, "PsiImag");
      H5Dclose(dset_imag);

      H5Fclose(file_id);

      global_psi = (fft_complex *)mymalloc("fdm_io_global_psi", total * sizeof(fft_complex));
      for(size_t idx = 0; idx < total; idx++)
        {
          global_psi[idx][0] = (fft_real)psi_real[idx];
          global_psi[idx][1] = (fft_real)psi_imag[idx];
        }
      /* psi_imag/psi_real are logically done with here, but CANNOT be
       * freed yet -- global_psi (and sendcounts/displs, allocated
       * below) were allocated AFTER them, so mymalloc/myfree's strict
       * LIFO discipline requires those be freed FIRST. Deferred to the
       * end of this function, in exact reverse allocation order. Same
       * mistake already made once before in fdm_update_potential()
       * (freeing mid-function while a later-allocated buffer was still
       * alive) -- caught by mymalloc's own safety check again here,
       * not a silent corruption, but a real repeat of a lesson that
       * should have been internalized after the first time. */

      sendcounts = (int *)mymalloc("fdm_io_sendcounts", NTask * sizeof(int));
      displs     = (int *)mymalloc("fdm_io_displs", NTask * sizeof(int));
      for(int t = 0; t < NTask; t++)
        {
          sendcounts[t] = FDM_plan.slabs_x_per_task[t] * N * N * 2;
          displs[t]     = FDM_plan.first_slab_x_of_task[t] * N * N * 2;
        }

      mpi_printf("FDM: read field from %s (N=%d)\n", fname, N);
    }

  int recvcount = FDM_plan.nslab_x * N * N * 2;
  MPI_Scatterv(global_psi, sendcounts, displs, FDM_MPI_REAL_TYPE, FDM_Psi, recvcount, FDM_MPI_REAL_TYPE, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      myfree(displs);
      myfree(sendcounts);
      myfree(global_psi);
      myfree(psi_imag);
      myfree(psi_real);
    }
}

#endif /* #ifdef FDM */
