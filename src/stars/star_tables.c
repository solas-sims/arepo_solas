#include <hdf5.h>

#include "../main/allvars.h"
#include "../main/proto.h"


// Main Sequence
int Z_COUNT = 0;
int M_COUNT = 0;

double *Z_VALUES = NULL;
double *M_VALUES = NULL;

double *logZ_VALUES = NULL;
double *logM_VALUES = NULL;
   
int **N = NULL;

double ***Age = NULL;
double ***Radius = NULL;
double ***Temperature = NULL;

#ifdef WINDS
// Winds
double ***MassLossRate = NULL;

#if GRACKLE_CHEMISTRY >= 1
double ***HLossRate = NULL;
double ***HeLossRate = NULL;
#endif
#ifdef METALS
double ***MetalsLossRate = NULL;
#endif
double ***WindVelocity = NULL;
#endif

#ifdef STAR_RADIATION_ACTIVE
// Radiation
WavebandData ***Flux[WAVEBANDS] = {0};
#endif 

#ifdef SUPERNOVAE
// Supernovae 
double **SN_MassLoss = NULL;
#if GRACKLE_CHEMISTRY >= 1
double **SN_HLoss = NULL;
double **SN_HeLoss = NULL;
#endif  
#ifdef METALS  
double **SN_MetalsLoss = NULL;
#endif
#endif  

#ifdef AGB 
// Asymptotic Giant Branch
double **AGB_MassLoss; 
#ifdef METALS
double **AGB_MetalsLoss; 
#endif 
#endif

void free_stellar_tables(void)
{
  if(N)
    {
      for(int z = 0; z < Z_COUNT; z++)
        for(int m = 0; m < M_COUNT; m++)
          {
            free(Age[z][m]);
            free(Radius[z][m]);
            free(Temperature[z][m]);

#ifdef WINDS
            free(MassLossRate[z][m]);
#if GRACKLE_CHEMISTRY >= 1
            free(HLossRate[z][m]);
            free(HeLossRate[z][m]);
#endif
#ifdef METALS
            free(MetalsLossRate[z][m]);
#endif
            free(WindVelocity[z][m]);
#endif

#ifdef STAR_RADIATION_ACTIVE
            for(int w = 0; w < WAVEBANDS; w++)
              free(Flux[w][z][m]);
#endif
          }

      for(int z = 0; z < Z_COUNT; z++)
        {
          free(N[z]);
          
          free(Age[z]);
          free(Radius[z]);
          free(Temperature[z]);

#ifdef WINDS
          free(MassLossRate[z]);
#if GRACKLE_CHEMISTRY >= 1
          free(HLossRate[z]);
          free(HeLossRate[z]);
#endif
#ifdef METALS
          free(MetalsLossRate[z]);
#endif
          free(WindVelocity[z]);
#endif

#ifdef STAR_RADIATION_ACTIVE
          for(int w = 0; w < WAVEBANDS; w++)
            free(Flux[w][z]);
#endif

#ifdef SUPERNOVAE
          free(SN_MassLoss[z]);
#if GRACKLE_CHEMISTRY >= 1
          free(SN_HLoss[z]);
          free(SN_HeLoss[z]);
#endif
#ifdef METALS
          free(SN_MetalsLoss[z]);
#endif
#endif
        }

      free(Z_VALUES);
      free(M_VALUES);

      free(logZ_VALUES);
      free(logM_VALUES);
      
      free(N);
      
      free(Age);
      free(Radius);
      free(Temperature);

#ifdef WINDS
      free(MassLossRate);
#if GRACKLE_CHEMISTRY >= 1
      free(HLossRate);
      free(HeLossRate);
#endif
#ifdef METALS
      free(MetalsLossRate);
#endif
      free(WindVelocity);
#endif

#ifdef STAR_RADIATION_ACTIVE
      for(int w = 0; w < WAVEBANDS; w++)
        free(Flux[w]);
#endif

#ifdef SUPERNOVAE
      free(SN_MassLoss);
#if GRACKLE_CHEMISTRY >= 1
      free(SN_HLoss);
      free(SN_HeLoss);
#endif
#ifdef METALS
      free(SN_MetalsLoss);
#endif
#endif

      Z_VALUES = NULL;
      M_VALUES = NULL;

      logZ_VALUES = NULL;
      logM_VALUES = NULL;
      
      N = NULL;
      
      Age = NULL;
      Radius = NULL;
      Temperature = NULL;

#ifdef WINDS
      MassLossRate = NULL;
#if GRACKLE_CHEMISTRY >= 1
      HLossRate = NULL;
      HeLossRate = NULL;
#endif
#ifdef METALS
      MetalsLossRate = NULL;
#endif
      WindVelocity = NULL;
#endif

#ifdef STAR_RADIATION_ACTIVE
      for(int w = 0; w < WAVEBANDS; w++)
        Flux[w] = NULL;
#endif    

#ifdef SUPERNOVAE
      SN_MassLoss = NULL;
#if GRACKLE_CHEMISTRY >= 1
      SN_HLoss = NULL;
      SN_HeLoss = NULL;
#endif
#ifdef METALS
      SN_MetalsLoss = NULL;
#endif
#endif
    }
}

void load_star_tables(const char *filename)
{ 
  hid_t file_id = -1;

  if(ThisTask == 0)
    {
      file_id = my_H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);  
      
      hid_t zc = my_H5Aopen_name(file_id, "Z_COUNT");
      hid_t mc = my_H5Aopen_name(file_id, "M_COUNT");
      
      my_H5Aread(zc, H5T_NATIVE_INT, &Z_COUNT, "Z_COUNT", 1);
      my_H5Aread(mc, H5T_NATIVE_INT, &M_COUNT, "M_COUNT", 1);

      my_H5Aclose(zc, "Z_COUNT");
      my_H5Aclose(mc, "M_COUNT");

      Z_VALUES = malloc(Z_COUNT * sizeof(double));
      M_VALUES = malloc(M_COUNT * sizeof(double));

      logZ_VALUES = malloc(Z_COUNT * sizeof(double));
      logM_VALUES = malloc(M_COUNT * sizeof(double));
      
      hid_t zv = my_H5Dopen(file_id, "Z_VALUES");
      hid_t mv = my_H5Dopen(file_id, "M_VALUES");

      my_H5Dread(zv, H5T_NATIVE_DOUBLE, 
              H5S_ALL, H5S_ALL, H5P_DEFAULT, Z_VALUES, "Z_VALUES");
      my_H5Dread(mv, H5T_NATIVE_DOUBLE, 
              H5S_ALL, H5S_ALL, H5P_DEFAULT, M_VALUES, "M_VALUES");

      for(int z = 0; z < Z_COUNT; z++)
        {
          if(z > 0 && Z_VALUES[z] <= Z_VALUES[z - 1])
            terminate("Z_VALUES not strictly increasing at z=%d (%g <= %g)", z, Z_VALUES[z], Z_VALUES[z - 1]);
          if(Z_VALUES[z] <= 0.0)
            terminate("Z_VALUES[%d] = %g; table must not contain Z <= 0", z, Z_VALUES[z]);
          logZ_VALUES[z] = log10(Z_VALUES[z]);
        }

      for(int m = 0; m < M_COUNT; m++)
        {
          if(m > 0 && M_VALUES[m] <= M_VALUES[m - 1])
            terminate("M_VALUES not strictly increasing at m=%d", m);
          if(M_VALUES[m] <= 0.0)
            terminate("M_VALUES[%d] = %g; table must not contain M <= 0", m, M_VALUES[m]);
          logM_VALUES[m] = log10(M_VALUES[m]);
        }

      my_H5Dclose(zv, "Z_VALUES");
      my_H5Dclose(mv, "M_VALUES");
    }

  MPI_Bcast(&Z_COUNT, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&M_COUNT, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(ThisTask != 0)
    {
      Z_VALUES = malloc(Z_COUNT * sizeof(double));
      M_VALUES = malloc(M_COUNT * sizeof(double));

      logZ_VALUES = malloc(Z_COUNT * sizeof(double));
      logM_VALUES = malloc(M_COUNT * sizeof(double));
    }

  MPI_Bcast(Z_VALUES, Z_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(M_VALUES, M_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Bcast(logZ_VALUES, Z_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(logM_VALUES, M_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  N = malloc(Z_COUNT * sizeof(int*));

  Age = malloc(Z_COUNT * sizeof(double**));
  Radius = malloc(Z_COUNT * sizeof(double**));
  Temperature = malloc(Z_COUNT * sizeof(double**));

#ifdef WINDS
  MassLossRate = malloc(Z_COUNT * sizeof(double**));
#if GRACKLE_CHEMISTRY >= 1
  HLossRate = malloc(Z_COUNT * sizeof(double**));
  HeLossRate = malloc(Z_COUNT * sizeof(double**));
#endif
#ifdef METALS
  MetalsLossRate = malloc(Z_COUNT * sizeof(double**));
#endif
  WindVelocity = malloc(Z_COUNT * sizeof(double**));
#endif

#ifdef STAR_RADIATION_ACTIVE
  for(int w = 0; w < WAVEBANDS; w++)
    Flux[w] = malloc(Z_COUNT * sizeof(WavebandData**));
#endif

#ifdef SUPERNOVAE
  SN_MassLoss = malloc(Z_COUNT * sizeof(double *));
#if GRACKLE_CHEMISTRY >= 1
  SN_HLoss = malloc(Z_COUNT * sizeof(double *));
  SN_HeLoss = malloc(Z_COUNT * sizeof(double *));
#endif
#ifdef METALS
  SN_MetalsLoss = malloc(Z_COUNT * sizeof(double *));
#endif
#endif 

  for(int z = 0; z < Z_COUNT; z++)
    {
      N[z] = malloc(M_COUNT * sizeof(int));

      Age[z] = malloc(M_COUNT * sizeof(double*));
      Radius[z] = malloc(M_COUNT * sizeof(double*));
      Temperature[z] = malloc(M_COUNT * sizeof(double*));

#ifdef WINDS
      MassLossRate[z] = malloc(M_COUNT * sizeof(double*));
#if GRACKLE_CHEMISTRY >= 1
      HLossRate[z] = malloc(M_COUNT * sizeof(double*));
      HeLossRate[z] = malloc(M_COUNT * sizeof(double*));
#endif
#ifdef METALS
      MetalsLossRate[z] = malloc(M_COUNT * sizeof(double*));
#endif
      WindVelocity[z] = malloc(M_COUNT * sizeof(double*));
#endif

#ifdef STAR_RADIATION_ACTIVE
      for(int w = 0; w < WAVEBANDS; w++)
        Flux[w][z] = malloc(M_COUNT * sizeof(WavebandData*));
#endif

#ifdef SUPERNOVAE
      SN_MassLoss[z] = malloc(M_COUNT * sizeof(double));
#if GRACKLE_CHEMISTRY >= 1
      SN_HLoss[z] = malloc(M_COUNT * sizeof(double));
      SN_HeLoss[z] = malloc(M_COUNT * sizeof(double));
#endif
#ifdef METALS
      SN_MetalsLoss[z] = malloc(M_COUNT * sizeof(double));
#endif
#endif 
    }

  if(ThisTask == 0)
    {
      for (int z = 0; z < Z_COUNT; z++)
        {
          char zname[64];
          snprintf(zname, sizeof(zname), "Z=%g", Z_VALUES[z]);

          hid_t zgrp = my_H5Gopen(file_id, zname);

          for(int m = 0; m < M_COUNT; m++)
            {
              char mname[64];
              snprintf(mname, sizeof(mname), "M=%03d", (int)round((M_VALUES[m])));

              if (H5Lexists(zgrp, mname, H5P_DEFAULT) <= 0)
                {
                  terminate("Error loading stellar tables!");
                }

              hid_t mgrp = my_H5Gopen(zgrp, mname);
              
              hid_t d_age = my_H5Dopen(mgrp, "Age");
              hid_t d_rad = my_H5Dopen(mgrp, "Radius");
              hid_t d_tem = my_H5Dopen(mgrp, "Temperature");

              hsize_t dims[1];
          
              hid_t space = H5Dget_space(d_age);
              H5Sget_simple_extent_dims(space, dims, NULL);
              H5Sclose(space);
              
              N[z][m] = (int)dims[0];

              Age[z][m] = malloc(N[z][m] * sizeof(double));
              Radius[z][m] = malloc(N[z][m] * sizeof(double));
              Temperature[z][m] = malloc(N[z][m] * sizeof(double));

              my_H5Dread(d_age, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, Age[z][m], "Age");          
              my_H5Dread(d_rad, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, Radius[z][m], "Radius");
              my_H5Dread(d_tem, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, Temperature[z][m], "Temperature");
              
              my_H5Dclose(d_age, "Age");
              my_H5Dclose(d_rad, "Radius");
              my_H5Dclose(d_tem, "Temperature");

#ifdef WINDS
              hid_t d_ml = my_H5Dopen(mgrp, "MassLossRate");
#if GRACKLE_CHEMISTRY >= 1
              hid_t d_Hl  = my_H5Dopen(mgrp, "HLossRate");
              hid_t d_Hel = my_H5Dopen(mgrp, "HeLossRate");
#endif
#ifdef METALS
              hid_t d_mz = my_H5Dopen(mgrp, "MetalsLossRate");
#endif
              hid_t d_wv = my_H5Dopen(mgrp, "WindVelocity");

              MassLossRate[z][m] = malloc(N[z][m] * sizeof(double));
#if GRACKLE_CHEMISTRY >= 1
              HLossRate[z][m] = malloc(N[z][m] * sizeof(double));
              HeLossRate[z][m] = malloc(N[z][m] * sizeof(double));
#endif
#ifdef METALS
              MetalsLossRate[z][m] = malloc(N[z][m] * sizeof(double));
#endif
              WindVelocity[z][m] = malloc(N[z][m] * sizeof(double));

              my_H5Dread(d_ml, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, MassLossRate[z][m], "MassLossRate");
#if GRACKLE_CHEMISTRY >= 1
              my_H5Dread(d_Hl, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, HLossRate[z][m], "HLossRate");
              my_H5Dread(d_Hel, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, HeLossRate[z][m], "HeLossRate");
#endif
#ifdef METALS
              my_H5Dread(d_mz, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, MetalsLossRate[z][m], "MetalsLossRate");
#endif
              my_H5Dread(d_wv, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, WindVelocity[z][m], "WindVelocity");

              my_H5Dclose(d_ml, "MassLossRate");
#if GRACKLE_CHEMISTRY >= 1
              my_H5Dclose(d_Hl, "HLossRate");
              my_H5Dclose(d_Hel, "HeLossRate");
#endif
#ifdef METALS
              my_H5Dclose(d_mz, "MetalsLossRate");
#endif
              my_H5Dclose(d_wv, "WindVelocity");
#endif

#ifdef STAR_RADIATION_ACTIVE
              hid_t d_energy  = my_H5Dopen(mgrp, "Energy");
              hid_t d_photons = my_H5Dopen(mgrp, "Photons");

              double (*energy_buf)[WAVEBANDS]  = malloc(N[z][m] * sizeof(*energy_buf));
              double (*photon_buf)[WAVEBANDS]  = malloc(N[z][m] * sizeof(*photon_buf));

              my_H5Dread(d_energy,  H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, energy_buf,  "Energy");
              my_H5Dread(d_photons, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, photon_buf,  "Photons");

              my_H5Dclose(d_energy,  "Energy");
              my_H5Dclose(d_photons, "Photons");

              for(int w = 0; w < WAVEBANDS; w++)
                {
                  Flux[w][z][m] = malloc(N[z][m] * sizeof(WavebandData));
                  for(int i = 0; i < N[z][m]; i++)
                    {
                      Flux[w][z][m][i].Energy  = energy_buf[i][w];
                      Flux[w][z][m][i].Photons = photon_buf[i][w];
                    }
                }

              free(energy_buf);
              free(photon_buf);
#endif

#ifdef SUPERNOVAE
              hid_t d_snml = my_H5Dopen(mgrp, "SN_MassLoss");
#if GRACKLE_CHEMISTRY >= 1
              hid_t d_snHl = my_H5Dopen(mgrp, "SN_HLoss");
              hid_t d_snHel = my_H5Dopen(mgrp, "SN_HeLoss");
#endif
#ifdef METALS
              hid_t d_snmz = my_H5Dopen(mgrp, "SN_MetalsLoss");
#endif
              my_H5Dread(d_snml, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, &SN_MassLoss[z][m], "SN_MassLoss");
#if GRACKLE_CHEMISTRY >= 1
              my_H5Dread(d_snHl, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, &SN_HLoss[z][m], "SN_HLoss");
              my_H5Dread(d_snHel, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, &SN_HeLoss[z][m], "SN_HeLoss");
#endif
#ifdef METALS
              my_H5Dread(d_snmz, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, &SN_MetalsLoss[z][m], "SN_MetalsLoss");
#endif
              my_H5Dclose(d_snml, "SN_MassLoss");
#if GRACKLE_CHEMISTRY >= 1
              my_H5Dclose(d_snHl, "SN_HLoss");
              my_H5Dclose(d_snHel, "SN_HeLoss");
#endif
#ifdef METALS
              my_H5Dclose(d_snmz, "SN_MetalsLoss");
#endif
#endif
              my_H5Gclose(mgrp, mname);
            }
          my_H5Gclose(zgrp, zname);
        }
      my_H5Fclose(file_id, filename);
    }

  for(int z = 0; z < Z_COUNT; z++)
    MPI_Bcast(N[z], M_COUNT, MPI_INT, 0, MPI_COMM_WORLD);

  for(int z = 0; z < Z_COUNT; z++)
    {
#ifdef SUPERNOVAE
      MPI_Bcast(SN_MassLoss[z], M_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#if GRACKLE_CHEMISTRY >= 1
      MPI_Bcast(SN_HLoss[z], M_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(SN_HeLoss[z], M_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
#ifdef METALS
      MPI_Bcast(SN_MetalsLoss[z], M_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
#endif
      for(int m = 0; m < M_COUNT; m++)
        if(N[z][m] > 0)
          {
            if(ThisTask != 0)
              {
                Age[z][m] = malloc(N[z][m] * sizeof(double));
                Radius[z][m] = malloc(N[z][m] * sizeof(double));
                Temperature[z][m] = malloc(N[z][m] * sizeof(double));

#ifdef WINDS
                MassLossRate[z][m] = malloc(N[z][m] * sizeof(double));
#if GRACKLE_CHEMISTRY >= 1
                HLossRate[z][m] = malloc(N[z][m] * sizeof(double));
                HeLossRate[z][m] = malloc(N[z][m] * sizeof(double));
#endif
#ifdef METALS
                MetalsLossRate[z][m] = malloc(N[z][m] * sizeof(double));
#endif
                WindVelocity[z][m] = malloc(N[z][m] * sizeof(double));
#endif

#ifdef STAR_RADIATION_ACTIVE
                for(int w = 0; w < WAVEBANDS; w++)
                  Flux[w][z][m] = malloc(N[z][m] * sizeof(WavebandData));
#endif
              }

            MPI_Bcast(Age[z][m], N[z][m], MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(Radius[z][m], N[z][m], MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(Temperature[z][m], N[z][m], MPI_DOUBLE, 0, MPI_COMM_WORLD);         

#ifdef WINDS
            MPI_Bcast(MassLossRate[z][m], N[z][m], MPI_DOUBLE, 0, MPI_COMM_WORLD);
#if GRACKLE_CHEMISTRY >= 1
            MPI_Bcast(HLossRate[z][m], N[z][m], MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(HeLossRate[z][m], N[z][m], MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
#ifdef METALS
            MPI_Bcast(MetalsLossRate[z][m], N[z][m], MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
            MPI_Bcast(WindVelocity[z][m], N[z][m], MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

#ifdef STAR_RADIATION_ACTIVE
            for(int w = 0; w < WAVEBANDS; w++)
              MPI_Bcast(Flux[w][z][m], N[z][m] * sizeof(WavebandData), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
          }
    }
}