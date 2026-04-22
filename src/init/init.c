/*!
 * \copyright   This file is part of the public version of the AREPO code.
 * \copyright   Copyright (C) 2009-2019, Max-Planck Institute for Astrophysics
 * \copyright   Developed by Volker Springel (vspringel@MPA-Garching.MPG.DE) and
 *              contributing authors.
 * \copyright   Arepo is free software: you can redistribute it and/or modify
 *              it under the terms of the GNU General Public License as published by
 *              the Free Software Foundation, either version 3 of the License, or
 *              (at your option) any later version.
 *
 *              Arepo is distributed in the hope that it will be useful,
 *              but WITHOUT ANY WARRANTY; without even the implied warranty of
 *              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *              GNU General Public License for more details.
 *
 *              A copy of the GNU General Public License is available under
 *              LICENSE as part of this program.  See also
 *              <https://www.gnu.org/licenses/>.
 *
 * \file        src/init/init.c
 * \date        05/2018
 * \brief       Initialization of a simulation from initial conditions.
 * \details     contains functions:
 *                int init(void)
 *                void check_omega(void)
 *                void setup_smoothinglengths(void)
 *                void setup_smoothinglengths_particles(void) *
 *                void test_id_uniqueness(void)
 *                void calculate_maxid(void)
 *                int compare_IDs(const void *a, const void *b)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 04.05.2018 Prepared file for public release -- Rainer Weinberger
 *  - 10/10/2025 Added functionality to compute initial smoothing lengths for stars and black holes - Chris Power
 */

#include <gsl/gsl_sf_gamma.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"
#include "../mesh/voronoi/voronoi.h"

#include "../fof/fof_seeding.h"

#include "../../celib/src/config.h"

/*! \brief Prepares the loaded initial conditions for the run.
 *
 *  It is only called if RestartFlag !=1. Various counters and variables are
 *  initialized. Entries of the particle data structures not read from initial
 *  conditions are initialized or converted and a initial domain decomposition
 *  is performed. If gas cells are present, the initial SPH smoothing lengths
 *  are determined.
 *
 *  \return status code: <0 if finished without errors and run can start,
 *          0 code ends after calling init()   0 an error occurred, terminate.
 */
int init(void)
{
  int i, j;
  double mass;

  assert(RestartFlag != 1);

  if(All.ComovingIntegrationOn)
    if(All.PeriodicBoundariesOn == 1)
      {
        if(RestartFlag < 3)
          /* can't do this check when not all particles are loaded */
          check_omega();
        else
          mpi_printf("INIT: Skipping Omega check since we are not doing a dynamical evolution (not all particles may be loaded)\n");
      }

/* IonizeParams initalises the Treecool file and such, 
* which we don't need when using grackle
*/
#if defined(COOLING) && !defined(USE_GRACKLE)
    IonizeParams();
#endif /* defined(COOLING) && !defined(USE_GRACKLE) */

  if(All.ComovingIntegrationOn)
    {
      All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
      All.Ti_Current        = 0;
    }
  else
    {
      All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
      All.Ti_Current        = 0;
    }

  set_cosmo_factors_for_current_time();

  for(j = 0; j < 3; j++)
    All.GlobalDisplacementVector[j] = 0;

  All.NumCurrentTiStep  = 0; /* setup some counters */
  All.SnapshotFileCount = 0;

  if(RestartFlag == 2)
    {
      if(RestartSnapNum < 0)
        All.SnapshotFileCount = atoi(All.InitCondFile + strlen(All.InitCondFile) - 3) + 1;
      else
        All.SnapshotFileCount = RestartSnapNum + 1;
    }

  All.TotNumOfForces     = 0;
  All.TopNodeAllocFactor = 0.08;
  All.TreeAllocFactor    = 0.7;
  All.NgbTreeAllocFactor = 0.7;

  if(NumPart < 1000)
    All.TreeAllocFactor = 10.0;

  DeRefMesh.Indi.AllocFacNdp = MIN_ALLOC_NUMBER;
  DeRefMesh.Indi.AllocFacNdt = MIN_ALLOC_NUMBER;

  Mesh.Indi.AllocFacNdp = 1.2 * NumGas + MIN_ALLOC_NUMBER;
  Mesh.Indi.AllocFacNdt = 8.0 * NumGas + MIN_ALLOC_NUMBER;
  Mesh.Indi.AllocFacNvf = 8.0 * NumGas + MIN_ALLOC_NUMBER;

  Mesh.Indi.AllocFacNvc = 16.0 * NumGas + MIN_ALLOC_NUMBER;
  Nvc                   = 0;

  Mesh.Indi.AllocFacNinlist     = 1.2 * NumGas + MIN_ALLOC_NUMBER;
  Mesh.Indi.AllocFacN_DP_Buffer = 0.2 * NumGas + MIN_ALLOC_NUMBER;
  Mesh.Indi.AllocFacNflux       = 0.01 * NumGas + MIN_ALLOC_NUMBER;
  Mesh.Indi.AllocFacNradinflux  = 0.01 * NumGas + MIN_ALLOC_NUMBER;

#ifdef MHD_POWELL
  for(j = 0; j < 3; j++)
    {
      All.Powell_Momentum[j]         = 0;
      All.Powell_Angular_Momentum[j] = 0;
    }
  All.Powell_Energy = 0;
#endif /* #ifdef MHD_POWELL */

  All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;

  set_softenings();

#ifdef ADAPTIVE_HYDRO_SOFTENING
  mpi_printf("INIT: Adaptive hydro softening, minimum gravitational softening for cells: %g\n", All.MinimumComovingHydroSoftening);
  mpi_printf("INIT: Adaptive hydro softening, maximum gravitational softening for cells: %g\n",
             All.MinimumComovingHydroSoftening * pow(All.AdaptiveHydroSofteningSpacing, NSOFTTYPES_HYDRO - 1));
  mpi_printf("INIT: Adaptive hydro softening, number of softening values: %d\n", NSOFTTYPES_HYDRO);
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
  init_individual_softenings();
#endif /* #ifdef INDIVIDUAL_GRAVITY_SOFTENING */

#ifdef SHIFT_BY_HALF_BOX
  for(i = 0; i < NumPart; i++)
    for(j = 0; j < 3; j++)
      P[i].Pos[j] += 0.5 * All.BoxSize;
#endif /* #ifdef SHIFT_BY_HALF_BOX */

  for(i = 0; i < GRAVCOSTLEVELS; i++)
    All.LevelToTimeBin[i] = -1;

  for(i = 0; i < NumPart; i++)
    for(j = 0; j < GRAVCOSTLEVELS; j++)
      P[i].GravCost[j] = 0;

      /* set unused coordinate values in 1d and 2d simulations to zero; this is needed for correct interfaces */
  int nonzero_vel = 0;
#ifdef ONEDIMS
  for(i = 0; i < NumPart; i++)
    {
      P[i].Pos[1] = 0.0;
      P[i].Pos[2] = 0.0;

      if(P[i].Vel[1] != 0.0 || P[i].Vel[2] != 0.0)
      {
   	    nonzero_vel = 1;
      }
    }
  if(nonzero_vel > 0)
  {
    warn("Initial y or z velocity nonzero in 1d simulation! Make sure you really want this!");
  }
#endif /* #ifdef ONEDIMS */

#ifdef TWODIMS
  for(i = 0; i < NumPart; i++)
    {
      P[i].Pos[2] = 0;

      if(P[i].Vel[2] != 0.0)
      {
        nonzero_vel = 1;
      }
    }
  if(nonzero_vel > 0)
  {
	warn("Initial z velocity nonzero in 2d simulation! Make sure you really want this!");
  }
#endif /* #ifdef TWODIMS */

  if(All.ComovingIntegrationOn) /*  change to new velocity variable */
    {
      for(i = 0; i < NumPart; i++)
        {
          for(j = 0; j < 3; j++)
            P[i].Vel[j] *= sqrt(All.Time) * All.Time; /* for dm/gas particles, p = a^2 xdot */
        }
    }

  /* measure mean cell mass */
  int num = 0;
  long long glob_num;
  double glob_mass;
  mass = 0;

  for(i = 0; i < NumGas; i++)
#ifdef REFINEMENT_HIGH_RES_GAS
    if(SphP[i].AllowRefinement != 0)
#endif /* #ifdef REFINEMENT_HIGH_RES_GAS */
      {
        num += 1;
        mass += P[i].Mass;
      }

  sumup_large_ints(1, &num, &glob_num);
  MPI_Allreduce(&mass, &glob_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#ifndef REFINEMENT_HIGH_RES_GAS
  if(glob_num != All.TotNumGas)
    terminate("glob_num(=%lld) != All.TotNumGas(=%lld)", glob_num, All.TotNumGas);
#endif /* #ifndef REFINEMENT_HIGH_RES_GAS */

  if(All.TotNumGas > 0 && (glob_num == 0 || glob_mass == 0))
    terminate("All.TotNumGas(=%lld) > 0 && (glob_num(=%lld) == 0 || glob_mass(=%g) == 0)", All.TotNumGas, glob_num, glob_mass);

    /* assign global variables that depend on the mean cell mass */
#if defined(REFINEMENT)
  if(All.ReferenceGasPartMass == 0)
    {
      if(!All.ComovingIntegrationOn)
        terminate("In non-comoving runs, ReferenceGasPartMass must be set to a non-zero value");

      All.ReferenceGasPartMass = glob_mass / glob_num;

      mpi_printf("REFINEMENT: The mean cell mass, which is used as a reference, is %g\n", All.ReferenceGasPartMass);
    }
  else
    mpi_printf("REFINEMENT: The given reference cell mass is %g\n", All.ReferenceGasPartMass);
  All.TargetGasMass = All.TargetGasMassFactor * All.ReferenceGasPartMass;
  mpi_printf("REFINEMENT: setting All.TargetGasMass=%g\n", All.TargetGasMass);
#endif /* #if defined(REFINEMENT) */

  for(i = 0; i < TIMEBINS; i++)
    All.Ti_begstep[i] = 0;

  for(i = 0; i < NumPart; i++) /*  start-up initialization */
    {
      for(j = 0; j < 3; j++)
        P[i].GravAccel[j] = 0;

#ifdef PMGRID
      for(j = 0; j < 3; j++)
        P[i].GravPM[j] = 0;
#endif /* #ifdef PMGRID */
      P[i].TimeBinHydro = 0;
      P[i].TimeBinGrav  = 0;
      P[i].OldAcc       = 0; /* Do not zero as masses are stored here */

#ifdef SELFGRAVITY
#ifdef EVALPOTENTIAL
      if(RestartFlag == 0)
        P[i].Potential = 0;
#endif /* #ifdef EVALPOTENTIAL */
#endif /* #ifdef SELFGRAVITY */

#ifdef USE_SFR
      if(RestartFlag == 0 && P[i].Type == 0)
        SphP[i].Sfr = 0;
#endif /* #ifdef USE_SFR */

#ifdef REFINEMENT_AROUND_BH
      if(P[i].Type == 0)
        {
          SphP[i].RefBHFlag = 0;
        #ifdef OUTPUT_REFBHCOUNTER
          SphP[i].RefBHCounter = 0;
        #endif //* #ifdef OUTPUT_REFBHCOUNTER */
        }
#endif /* #ifdef REFINEMENT_AROUND_BH*/
  }

  for(i = 0; i < TIMEBINS; i++)
    TimeBinSynchronized[i] = 1;

  reconstruct_timebins();

#ifdef PMGRID
  All.PM_Ti_endstep = All.PM_Ti_begstep = 0;
#endif /* #ifdef PMGRID */

  for(i = 0; i < NumGas; i++) /* initialize sph_properties */
    {
      if(RestartFlag == 2 || RestartFlag == 3)
        for(j = 0; j < 3; j++)
          SphP[i].Center[j] = P[i].Pos[j];

#if defined(CELL_CENTER_GRAVITY) && !defined(OUTPUT_CENTER_OF_MASS)
      if(RestartFlag == 17 || RestartFlag == 18)
        for(j = 0; j < 3; j++)
          SphP[i].Center[j] = P[i].Pos[j];
#endif /* #if defined(CELL_CENTER_GRAVITY) && !defined(OUTPUT_CENTER_OF_MASS) */

      if(RestartFlag == 0)
        {
          for(j = 0; j < 3; j++)
            SphP[i].Center[j] = P[i].Pos[j];

          SphP[i].Hsml = 0;
#if defined(COOLING)
          SphP[i].Ne = 1.0;
#endif /* #if defined(COOLING)  */
        }
    }

#ifndef NODEREFINE_BACKGROUND_GRID
  double mvol = 0;
  if(All.TotNumGas)
    {
#ifdef TWODIMS
      mvol = boxSize_X * boxSize_Y / All.TotNumGas;
#else /* #ifdef TWODIMS */
#ifdef ONEDIMS
      mvol                  = boxSize_X / All.TotNumGas;
#else  /* #ifdef ONEDIMS */
      mvol = boxSize_X * boxSize_Y * boxSize_Z / All.TotNumGas;
#endif /* #ifdef ONEDIMS #else */
#endif /* #ifdef TWODIMS #else */
    }

  All.MeanVolume = mvol;
#endif /* #ifndef NODEREFINE_BACKGROUND_GRID */

  mpi_printf("INIT: MeanVolume=%g\n", All.MeanVolume);

#ifndef NO_ID_UNIQUE_CHECK
  test_id_uniqueness();
#endif /* #ifndef NO_ID_UNIQUE_CHECK */

#ifdef REFINEMENT_MERGE_CELLS
  for(i = 0; i < NumPart; i++)
    if(P[i].Type == 0 && P[i].ID == 0)
      terminate("INIT: Cannot use ID==0 for gas in ICs with derefinement enabled.");
#endif /* #ifdef REFINEMENT_MERGE_CELLS */

  voronoi_init_connectivity(&Mesh);

#ifdef ADDBACKGROUNDGRID
  prepare_domain_backgroundgrid();
#endif /* #ifdef ADDBACKGROUNDGRID */

  domain_Decomposition(); /* do initial domain decomposition (gives equal numbers of particles) */

  if(RestartFlag == 18) /* recalculation of potential */
    {
      mark_active_timebins();
      open_logfiles();
#if defined(USE_SFR) && defined(EEOS_SF)
      sfr_init();
#endif /* #if defined(USE_SFR) */
      set_non_standard_physics_for_current_time();

#ifdef PMGRID
      long_range_init_regionsize();
#endif /* #ifdef PMGRID */

      compute_grav_accelerations(All.HighestActiveTimeBin, FLAG_FULL_TREE);

#if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT) && defined(FOF)
      PS = (struct subfind_data *)mymalloc_movable(&PS, "PS", All.MaxPart * sizeof(struct subfind_data));
      fof_prepare_output_order(); /* sort by type and Fileorder */
      fof_subfind_exchange(MPI_COMM_WORLD);
#endif /* #if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT) && defined(FOF) */

      sprintf(All.SnapshotFileBase, "%s_potupdated", All.SnapshotFileBase);
      mpi_printf("Start writing file %s\nRestartSnapNum %d\n", All.SnapshotFileBase, RestartSnapNum);
      savepositions(RestartSnapNum, 0);

      endrun();
    }

  /* will build tree */
  ngb_treeallocate();
  ngb_treebuild(NumGas);

#ifdef HALO_SEEDING
#ifndef FOF
#error "HALO_SEEDING is only implemented for FOF."
#endif /* #ifndef FOF */
  fof_seeding_init(RestartFlag);
#endif
 
  if(RestartFlag == 3)
    {
#ifdef FOF
      fof_fof(RestartSnapNum);
      DumpFlag = 1;
      savepositions(RestartSnapNum, 0);
#endif /* #ifdef FOF */
      return (0);
    }

  All.Ti_Current = 0;

  if(RestartFlag == 0 || RestartFlag == 2 || RestartFlag == 14 || RestartFlag == 17)
  {
      mpi_printf("Finding smoothing lengths on local processor %d\n", ThisTask);
      setup_smoothinglengths();
#if defined(STARS) || defined(BLACKHOLES)
      mpi_printf("Finding smoothing lengths for stars and black holes on local processor %d\n", ThisTask);
      setup_smoothinglengths_particles();
#endif /* #if defined(STARS) || defined(BLACKHOLES) */
  }

#ifdef ADDBACKGROUNDGRID
  // This return more clearly shows that this function terminates the run
  return add_backgroundgrid();
#endif /* #ifdef ADDBACKGROUNDGRID */

  create_mesh();
  mesh_setup_exchange();

  if(RestartFlag == 14)
    {
      char tess_name[1024];
      sprintf(tess_name, "%s/tess_%03d", All.OutputDir, RestartSnapNum);
      write_voronoi_mesh(&Mesh, tess_name, 0, NTask - 1);
      return 0;
    }

  for(i = 0, mass = 0; i < NumGas; i++)
    {
      if(RestartFlag == 0)
        {
#ifdef READ_MASS_AS_DENSITY_IN_INPUT
          P[i].Mass *= SphP[i].Volume;
#endif /* #ifdef READ_MASS_AS_DENSITY_IN_INPUT */
        }

      SphP[i].Density = P[i].Mass / SphP[i].Volume;

      if(SphP[i].Density < All.MinimumDensityOnStartUp)
        {
          SphP[i].Density = All.MinimumDensityOnStartUp;

          P[i].Mass = SphP[i].Volume * SphP[i].Density;
        }

      SphP[i].Momentum[0] = P[i].Mass * P[i].Vel[0];
      SphP[i].Momentum[1] = P[i].Mass * P[i].Vel[1];
      SphP[i].Momentum[2] = P[i].Mass * P[i].Vel[2];

        /* utherm has been loaded from IC file */
#ifdef MESHRELAX
      SphP[i].Energy = P[i].Mass * SphP[i].Utherm;
#else  /* #ifdef MESHRELAX */
      SphP[i].Energy = P[i].Mass * All.cf_atime * All.cf_atime * SphP[i].Utherm +
                       0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
#endif /* #ifdef MESHRELAX #else */

#ifdef MHD
      SphP[i].Energy += 0.5 * (SphP[i].B[0] * SphP[i].B[0] + SphP[i].B[1] * SphP[i].B[1] + SphP[i].B[2] * SphP[i].B[2]) *
                        SphP[i].Volume * All.cf_atime;
#endif /* #ifdef MHD */

      for(j = 0; j < 3; j++)
        SphP[i].VelVertex[j] = P[i].Vel[j];

      mass += P[i].Mass;
    }

/* NOTE: The metals have to be initialised before the PASSIVE_SCALARS.
 * The value in the PScalars are set to zero during reading ICs. If PConservedScalars are set to zero,
 * this the same as no advection!
 * */
#ifdef METALS
  for(i=0; i<NumGas; i++)
      SphP[i].Metals = All.InitMetallicityinSolar * SOLAR_ABUNDANCE;
#endif /* ifdef METALS */

#ifdef PASSIVE_SCALARS
  for(i = 0; i < NumGas; i++)
      for(j = 0; j < PASSIVE_SCALARS; j++)
        SphP[i].PConservedScalars[j] = SphP[i].PScalars[j] * P[i].Mass;
#endif /* #ifdef PASSIVE_SCALARS */

  if(RestartFlag == 17)
    {
      update_primitive_variables();
      exchange_primitive_variables();
      calculate_gradients();
      exchange_primitive_variables_and_gradients();
      DumpFlag = 1;
      savepositions(RestartSnapNum + 1, 0);
      return (0);
    }

  update_primitive_variables();

#ifdef TREE_BASED_TIMESTEPS
  tree_based_timesteps_setsoundspeeds();
#endif /* #ifdef TREE_BASED_TIMESTEPS */

  /* initialize star formation rate */
#if defined(USE_SFR) && defined(EEOS_SF)
  sfr_init();
#endif /* #if defined(USE_SFR) */

#if defined(USE_SFR)
  for(i = 0; i < NumGas; i++)
    SphP[i].Sfr = get_starformation_rate(i);
#endif /* #if defined(USE_SFR) */

  update_primitive_variables();

  exchange_primitive_variables();

  calculate_gradients();

  exchange_primitive_variables_and_gradients();

#if !defined(ONEDIMS) && !defined(TWODIMS)
  int xaxis, yaxis, zaxis, weight_flag = 0;
  double xmin, xmax, ymin, ymax, zmin, zmax;
#endif /* #if !defined(ONEDIMS) && !defined(TWODIMS) */

  free_mesh();

#if defined(STARS) || defined(BH_WITH_FEEDBACK)
  /* initialize feedback variables */
  All.FeedbackFlag = 1;
  All.EnergyExchange[0] = All.EnergyExchange[1] = 0;
  All.EnergyExchange[2] = All.EnergyExchange[3] = 0;
  All.EnergyExchange[4] = All.EnergyExchange[5] = 0;

  double *exch = All.EnergyExchangeTot;
  exch = malloc(6 * sizeof(double));
#endif 
/*
#ifdef STARS
   
  for(i=0; i<NumStars; i++)
  {
    SP[i].Metals = 0.0004;

    struct CELibStructNextEventTimeInput Input = 
      {
        .R = (double)rand()/(double)RAND_MAX,
        .InitialMass_in_Msun = (PPS(i).Mass * All.UnitMass_in_g / SOLAR_MASS),
        .Metallicity = SP[i].Metals
      };

    SP[i].SNIITime = CELibGetNextEventTime(Input, CELibFeedbackType_SNII) 
      / (1.e6) / All.UnitTime_in_Megayears;

    //timebin_add_particle(&TimeBinsStar, NumStars, -1, 0, 1);  
  }

#endif
*/
  return -1;  // return -1 means we ran to completion, i.e. not an endrun code
}

/*! \brief This routine computes the mass content of the box and compares it
 *         to the specified value of Omega-matter.
 *
 *  If discrepant, the run is terminated.
 *
 *  \return void
 */
void check_omega(void)
{
  double mass   = 0, masstot, omega;
  double mass_b = 0, masstot_b, omega_b;
  int i, n_b = 0;

  for(i = 0; i < NumPart; i++)
    {
      mass += P[i].Mass;
      if(P[i].Type == 0)
        {
          mass_b += P[i].Mass;
          n_b += 1;
        }
#ifdef USE_SFR
      if(P[i].Type == 4)
        {
          mass_b += P[i].Mass;
          n_b += 1;
        }
#endif /* #ifdef USE_SFR */
    }
  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&mass_b, &masstot_b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  omega   = masstot / (All.BoxSize * All.BoxSize * All.BoxSize) / (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G));
  omega_b = masstot_b / (All.BoxSize * All.BoxSize * All.BoxSize) / (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G));

  if(n_b > 0)
    {
      if(fabs((omega - All.Omega0) / omega) > 1.0e-1 || fabs((omega_b - All.OmegaBaryon) / omega_b) > 1.0e-1)
        {
#ifndef TWODIMS
          mpi_terminate(
              "\n\nI've found something odd!\nThe mass content accounts for Omega=%g and OmegaBaryon=%g,\nbut you specified Omega=%g "
              "and OmegaBaryon=%g in the parameterfile.\n\nI better stop.\n",
              omega, omega_b, All.Omega0, All.OmegaBaryon);
#endif /* #ifndef TWODIMS */
        }

      if(fabs((omega - All.Omega0) / omega) > 1.0e-3 || fabs((omega_b - All.OmegaBaryon) / omega_b) > 1.0e-3)
        if(ThisTask == 0)
          warn(
              "I've found something odd! The mass content accounts for Omega=%g and OmegaBaryon=%g, but you specified Omega=%g and "
              "OmegaBaryon=%g in the parameterfile.",
              omega, omega_b, All.Omega0, All.OmegaBaryon);
    }
  else
    {
      if(All.OmegaBaryon != 0)
        if(ThisTask == 0)
          warn(
              "We are running with no baryons, even though you have specified OmegaBaryon=%g in the parameterfile. Please make sure "
              "you really want this.\n\n",
              All.OmegaBaryon);

      if(fabs((omega - All.Omega0) / omega) > 1.0e-1)
        {
#ifndef TWODIMS
          mpi_terminate(
              "\n\nI've found something odd!\nThe mass content accounts for Omega=%g and OmegaBaryon=%g,\nbut you specified Omega=%g "
              "and OmegaBaryon=%g in the parameterfile.\n\nI better stop.\n",
              omega, omega_b, All.Omega0, All.OmegaBaryon);
#endif /* #ifndef TWODIMS */
        }

      if(fabs((omega - All.Omega0) / omega) > 1.0e-3)
        if(ThisTask == 0)
          warn(
              "I've found something odd! The mass content accounts for Omega=%g and OmegaBaryon=%g, but you specified Omega=%g and "
              "OmegaBaryon=%g in the parameterfile.",
              omega, omega_b, All.Omega0, All.OmegaBaryon);
    }
}

/*! \brief This function is used to find an initial SPH smoothing length for
 *         each cell.
 *
 *  It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the smoothing length is provided to the function density(), which will
 *  then iterate if needed to find the right smoothing length.
 *
 *  \return void
 */
void setup_smoothinglengths(void)
{
  int i, no, p;
  double *save_masses = mymalloc("save_masses", NumGas * sizeof(double));

  for(i = 0; i < NumGas; i++)
    {
#ifdef NO_GAS_SELFGRAVITY
      /* This is needed otherwise the force tree will not be constructed for gas particles */
      P[i].Type = -1;
#endif /* #ifdef NO_GAS_SELFGRAVITY */
      save_masses[i] = P[i].Mass;
      P[i].Mass      = 1.0;
    }

#ifdef HIERARCHICAL_GRAVITY
  TimeBinsGravity.NActiveParticles = 0;
  for(i = 0; i < NumGas; i++)
    {
      TimeBinsGravity.ActiveParticleList[TimeBinsGravity.NActiveParticles] = i;
      TimeBinsGravity.NActiveParticles++;
    }
#endif /* #ifdef HIERARCHICAL_GRAVITY */

  construct_forcetree(1, 1, 0, 0); /* build force tree with gas particles only */

  for(i = 0; i < NumGas; i++)
    {
      no = Father[i];

      if(no < 0)
        terminate("i=%d no=%d\n", i, no);

      while(10 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
        {
          p = Nodes[no].u.d.father;

          if(p < 0)
            break;

          no = p;
        }
#ifndef TWODIMS
      SphP[i].Hsml = pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;
#else  /* #ifndef TWODIMS */
      SphP[i].Hsml = pow(1.0 / (M_PI)*All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 2) * Nodes[no].len;
#endif /* #ifndef TWODIMS #else */
#ifdef NO_GAS_SELFGRAVITY
      /* Reset the original particle type */
      P[i].Type = 0;
#endif /* #ifdef NO_GAS_SELFGRAVITY */
    }

  myfree(Father);
  myfree(Nextnode);

  myfree(Tree_Points);
  force_treefree();

  density();

  for(i = 0; i < NumGas; i++)
    P[i].Mass = save_masses[i];

  myfree(save_masses);

  for(i = 0; i < NumGas; i++)
    SphP[i].MaxDelaunayRadius = SphP[i].Hsml;

#ifdef FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES
  ngb_treefree();
  domain_free();
  domain_Decomposition();
  ngb_treeallocate();
  ngb_treebuild(NumGas);
#endif /* #ifdef FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES */
}

/*! \brief This function is used to find an initial SPH smoothing length for
 *         stars and black holes
 *
 *  It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the smoothing length is provided to the function density(), which will
 *  then iterate if needed to find the right smoothing length.
 *
 *  \return void
 */
void setup_smoothinglengths_particles(void)
{
  int i, no, p;
  int DesNgb;    /* How many neighbours within Hsml?*/
  double Hsml;   /* Smoothing length */
    
  double *save_masses = mymalloc("save_masses", NumPart * sizeof(double));

  for(i = 0; i < NumPart; i++)
    {
      save_masses[i] = P[i].Mass;
      P[i].Mass      = 1.0;
    }

#ifdef HIERARCHICAL_GRAVITY
  TimeBinsGravity.NActiveParticles = 0;
  for(i = 0; i < NumPart; i++)
    {
      TimeBinsGravity.ActiveParticleList[TimeBinsGravity.NActiveParticles] = i;
      TimeBinsGravity.NActiveParticles++;
    }
#endif /* #ifdef HIERARCHICAL_GRAVITY */

  construct_forcetree(0, 1, 0, 0); /* build force tree with all particles only */

#ifdef STARS
  TimeBinsStar.NActiveParticles = 0;
#endif /* #ifdef STARS */
#ifdef BLACKHOLES
  TimeBinsBh.NActiveParticles = 0;
#endif /* #ifdef BLACKHOLES */

  for(i = 0; i < NumPart; i++)
    {
      if (P[i].Type != 4 && P[i].Type !=5) continue; /* only stars and BH */
        
      DesNgb=All.DesNumNgb;    /* Assumes as default in case STARS/BLACKHOLES not defined*/

#ifdef STARS
      if (P[i].Type==4) DesNgb=All.StarDesNgb;
#endif
#ifdef BLACKHOLES
      if (P[i].Type==5) DesNgb=All.BhDesNgb;
#endif

      no = Father[i];

      if(no < 0)
        terminate("i=%d no=%d\n", i, no);

      while(10 * DesNgb * P[i].Mass > Nodes[no].u.d.mass)
        {
          p = Nodes[no].u.d.father;

          if(p < 0)
            break;

          no = p;
        }

#ifdef STARS
      if(P[i].Type == 4) 
      {
#ifndef TWODIMS
        Hsml = pow(3.0 / (4 * M_PI) * All.StarDesNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;
#else  /* #ifndef TWODIMS */
        Hsml = pow(1.0 / (M_PI)*All.StarDesNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 2) * Nodes[no].len;
#endif /* #ifndef TWODIMS #else */
        TimeBinsStar.ActiveParticleList[TimeBinsStar.NActiveParticles] = TimeBinsStar.NActiveParticles;
        SP[TimeBinsStar.NActiveParticles].Hsml = Hsml;
        TimeBinsStar.NActiveParticles++;  
      }
#endif /* #ifdef STARS */
#ifdef BLACKHOLES
      if(P[i].Type == 5)
      {
#ifndef TWODIMS
        Hsml = pow(3.0 / (4 * M_PI) * All.BhDesNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;
#else  /* #ifndef TWODIMS */
        Hsml = pow(1.0 / (M_PI)*All.BhDesNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 2) * Nodes[no].len;
#endif /* #ifndef TWODIMS #else */
        TimeBinsBh.ActiveParticleList[TimeBinsBh.NActiveParticles] = TimeBinsBh.NActiveParticles;
        BhP[TimeBinsBh.NActiveParticles].Hsml = Hsml;
        TimeBinsBh.NActiveParticles++;  
      }
#endif /* #ifdef BLACKHOLES */
    }

  myfree(Father);
  myfree(Nextnode);

  myfree(Tree_Points);
  force_treefree();

#ifdef STARS
  star_density();
#endif

#ifdef BLACKHOLES
  bh_density();
#endif
    
  for(i = 0; i < NumPart; i++)
    P[i].Mass = save_masses[i];

  myfree(save_masses);
}

/*! \brief This function checks for unique particle IDs.
 *
 *  The particle IDs are copied to an array and then sorted among all tasks.
 *  This array is then checked for duplicates. In that case the code
 *  terminates.
 *
 *  \return void
 */
void test_id_uniqueness(void)
{
  int i;
  double t0, t1;
  MyIDType *ids, *ids_first;

  mpi_printf("INIT: Testing ID uniqueness...\n");

  if(NumPart == 0)
    terminate("need at least one particle per cpu\n");

  t0 = second();

  ids       = (MyIDType *)mymalloc("ids", NumPart * sizeof(MyIDType));
  ids_first = (MyIDType *)mymalloc("ids_first", NTask * sizeof(MyIDType));

  for(i = 0; i < NumPart; i++)
    ids[i] = P[i].ID;

  parallel_sort(ids, NumPart, sizeof(MyIDType), compare_IDs);

  for(i = 1; i < NumPart; i++)
    {
      if(ids[i] == ids[i - 1])
        terminate("non-unique ID=%lld found on task=%d (i=%d NumPart=%d)\n", (long long)ids[i], ThisTask, i, NumPart);
    }
  MPI_Allgather(&ids[0], sizeof(MyIDType), MPI_BYTE, ids_first, sizeof(MyIDType), MPI_BYTE, MPI_COMM_WORLD);

  if(ThisTask < NTask - 1)
    {
      if(ids[NumPart - 1] == ids_first[ThisTask + 1])
        terminate("non-unique ID=%lld found on task=%d\n", (long long)ids[NumPart - 1], ThisTask);
    }
  myfree(ids_first);
  myfree(ids);

  t1 = second();

  mpi_printf("INIT: success.  took=%g sec\n", timediff(t0, t1));
}

/*! \brief Calculates global maximum of the IDs of all particles.
 *
 *  This is needed for REFINEMENT_SPLIT_CELLS.
 *
 *  \return void
 */
void calculate_maxid(void)
{
  /* determine maximum ID */
  MyIDType maxid, *tmp;
  int i;

  for(i = 0, maxid = 0; i < NumPart; i++)
    if(P[i].ID > maxid)
      {
        maxid = P[i].ID;
      }

  tmp = mymalloc("tmp", NTask * sizeof(MyIDType));

  MPI_Allgather(&maxid, sizeof(MyIDType), MPI_BYTE, tmp, sizeof(MyIDType), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < NTask; i++)
    if(tmp[i] > maxid)
      maxid = tmp[i];

#if defined(REFINEMENT_SPLIT_CELLS) || defined(USE_SFR)
  All.MaxID = maxid;
#endif /* #if defined(REFINEMENT_SPLIT_CELLS) || defined(USE_SFR) */

  myfree(tmp);
}

/*! \brief Comparison function for two MyIDType objects.
 *
 *  Used as sorting-kernel for id_uniqueness check.
 *
 *  \return (-1,0,1), -1 if a<b, 0 if a==b, 1 if a>b
 */
int compare_IDs(const void *a, const void *b)
{
  if(*((MyIDType *)a) < *((MyIDType *)b))
    return -1;

  if(*((MyIDType *)a) > *((MyIDType *)b))
    return +1;

  return 0;
}
