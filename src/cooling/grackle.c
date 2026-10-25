#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

/* Array-based interface to Grackle.
 *
 * Grackle's field-data struct describes a 1D block of N cells (grid_dimension[0] = N), and solve_chemistry() /
 * calculate_*() loop over all of them in a single call. The routines below
 *
 *     grackle_gather()  -> copy N cells from SphP/P into the (persistent) Grackle field arrays
 *     Grackle call      -> one call for the whole block
 *     grackle_scatter() -> copy the evolved abundances back into SphP
 *
 * are used both by GrackleCoolActive() (batched over all active cells, one block per time bin because solve_chemistry() takes a
 * single scalar dt) and by CallGrackle() (a block of N=1, kept for callers that handle one cell at a time).
 */

/* persistent field arrays; (re)allocated only when a larger block is needed */
static grackle_field_data GF;
static int GFCapacity = 0;

#define GF_FIELDS_BASE(X) X(x_velocity) X(y_velocity) X(z_velocity) X(density) X(internal_energy) X(metal_density)

#if (GRACKLE_CHEMISTRY >= 1)
#define GF_FIELDS_CHEM(X)                                                                                                         \
  X(e_density) X(HI_density) X(HII_density) X(HeI_density) X(HeII_density) X(HeIII_density) X(H2I_density) X(H2II_density)       \
  X(HM_density) X(DI_density) X(DII_density) X(HDI_density) X(volumetric_heating_rate) X(RT_H2_dissociation_rate)                  \
  X(RT_HI_heating_rate) X(RT_HeI_heating_rate) X(RT_HeII_heating_rate) X(RT_HI_ionization_rate) X(RT_HeI_ionization_rate)         \
  X(RT_HeII_ionization_rate) X(specific_heating_rate)
#else
#define GF_FIELDS_CHEM(X)
#endif

#define GF_ALLOC(f) GF.f = malloc(newcap * sizeof(gr_float));
#define GF_FREE(f)  free(GF.f);

static void grackle_free_fields(void)
{
  if(GFCapacity == 0)
    return;

  free(GF.grid_dimension);
  free(GF.grid_start);
  free(GF.grid_end);
  GF_FIELDS_BASE(GF_FREE)
  GF_FIELDS_CHEM(GF_FREE)
  GFCapacity = 0;
}

/* make sure the field arrays can hold n cells and describe a 1D block of exactly n cells */
static void grackle_prepare_block(int n)
{
  if(n > GFCapacity)
    {
      int newcap = (n == 1) ? 1 : (int)(1.2 * n) + 16;
      grackle_free_fields();

      gr_initialize_field_data(&GF);
      GF.grid_dimension = malloc(3 * sizeof(int));
      GF.grid_start     = malloc(3 * sizeof(int));
      GF.grid_end       = malloc(3 * sizeof(int));
      GF_FIELDS_BASE(GF_ALLOC)
      GF_FIELDS_CHEM(GF_ALLOC)
      GFCapacity = newcap;
    }

  GF.grid_rank = 3;
  for(int k = 0; k < 3; k++)
    {
      GF.grid_dimension[k] = 1;
      GF.grid_start[k]     = 0;
      GF.grid_end[k]       = 0;
    }
  GF.grid_dimension[0] = n;
  GF.grid_end[0]       = n - 1;

  /* set the right scale factor */
  All.GrackleUnits.a_value = All.cf_atime;
}

/* Fill the field arrays for the n cells target[0..n-1], using the thermal energies u[] and the (proper) densities rho[]. If
 * consume_rates is set, the accumulated radiation rates stored in SphP are zeroed once they have been passed to Grackle; this must
 * only be done for the call that actually integrates the cell (solve_chemistry), not for diagnostic queries. */
static void grackle_gather(int n, const int *target, const double *u, const double *rho, int consume_rates)
{
  grackle_prepare_block(n);

  for(int k = 0; k < n; k++)
    {
      int i = target[k];

      GF.x_velocity[k]      = P[i].Vel[0];
      GF.y_velocity[k]      = P[i].Vel[1];
      GF.z_velocity[k]      = P[i].Vel[2];
      GF.density[k]         = rho[k];
      GF.internal_energy[k] = u[k];

#ifdef METALS
      double Metallicity = SphP[i].GasMetallicity;
#else
      double Metallicity = GRACKLE_TINY;
#endif
      GF.metal_density[k] = Metallicity * GF.density[k];

#if (GRACKLE_CHEMISTRY >= 1)
      /* H and He species */
      GF.HI_density[k]    = SphP[i].GrackleSpecies(GRACKLE_HI) * GF.density[k];
      GF.HII_density[k]   = SphP[i].GrackleSpecies(GRACKLE_HII) * GF.density[k];
      GF.HeI_density[k]   = SphP[i].GrackleSpecies(GRACKLE_HeI) * GF.density[k];
      GF.HeII_density[k]  = SphP[i].GrackleSpecies(GRACKLE_HeII) * GF.density[k];
      GF.HeIII_density[k] = SphP[i].GrackleSpecies(GRACKLE_HeIII) * GF.density[k];

      /* molecular H species */
#if (GRACKLE_CHEMISTRY >= 2)
      GF.H2I_density[k]  = SphP[i].GrackleSpecies(GRACKLE_H2I) * GF.density[k];
      GF.H2II_density[k] = SphP[i].GrackleSpecies(GRACKLE_H2II) * GF.density[k];
      GF.HM_density[k]   = SphP[i].GrackleSpecies(GRACKLE_HM) * GF.density[k];
#else
      GF.H2I_density[k]  = GRACKLE_TINY * GF.density[k];
      GF.H2II_density[k] = GRACKLE_TINY * GF.density[k];
      GF.HM_density[k]   = GRACKLE_TINY * GF.density[k];
#endif

      /* deuterium species */
#if (GRACKLE_CHEMISTRY >= 3)
      GF.DI_density[k]  = SphP[i].GrackleSpecies(GRACKLE_DI) * GF.density[k];
      GF.DII_density[k] = SphP[i].GrackleSpecies(GRACKLE_DII) * GF.density[k];
      GF.HDI_density[k] = SphP[i].GrackleSpecies(GRACKLE_HDI) * GF.density[k];
#else
      GF.DI_density[k]  = GRACKLE_TINY * GF.density[k];
      GF.DII_density[k] = GRACKLE_TINY * GF.density[k];
      GF.HDI_density[k] = GRACKLE_TINY * GF.density[k];
#endif

      /* electron density from charge balance */
      double e_density = GF.HII_density[k] + GF.HeII_density[k] / 4.0 + GF.HeIII_density[k] / 2.0;
#if (GRACKLE_CHEMISTRY >= 2)
      e_density += GF.H2II_density[k] / 2.0 - GF.HM_density[k];
#endif
#if (GRACKLE_CHEMISTRY >= 3)
      e_density += GF.DII_density[k] / 2.0;
#endif
      GF.e_density[k] = e_density;

      /* Radiation */
      /* volumetric heating rate in [erg s^-1 cm^-3]; H2 dissociation rate in [1/time_units]; RT heating rates in
       * [erg s^-1 cm^-3] / n; ionization rates in [1/time_units]; specific heating rate in [erg s^-1 g^-1] */
#ifdef PHOTOELECTRIC_HEATING
      GF.volumetric_heating_rate[k] = (gr_float)(SphP[i].PE_VolHeatingRate);
      if(consume_rates)
        SphP[i].PE_VolHeatingRate = 0;
#else
      GF.volumetric_heating_rate[k] = 0;
#endif

#ifdef DISSOCIATION
      GF.RT_H2_dissociation_rate[k] = (gr_float)(SphP[i].H2_DissociationRate);
      if(consume_rates)
        SphP[i].H2_DissociationRate = 0;
#else
      GF.RT_H2_dissociation_rate[k] = 0;
#endif

#ifdef PHOTOIONIZATION
      GF.RT_HI_heating_rate[k]      = (gr_float)(SphP[i].HI_HeatingRate);
      GF.RT_HeI_heating_rate[k]     = (gr_float)(SphP[i].HeI_HeatingRate);
      GF.RT_HeII_heating_rate[k]    = (gr_float)(SphP[i].HeII_HeatingRate);
      GF.RT_HI_ionization_rate[k]   = (gr_float)(SphP[i].HI_IonizationRate);
      GF.RT_HeI_ionization_rate[k]  = (gr_float)(SphP[i].HeI_IonizationRate);
      GF.RT_HeII_ionization_rate[k] = (gr_float)(SphP[i].HeII_IonizationRate);
      if(consume_rates)
        SphP[i].HI_HeatingRate = SphP[i].HeI_HeatingRate = SphP[i].HeII_HeatingRate = SphP[i].HI_IonizationRate =
            SphP[i].HeI_IonizationRate = SphP[i].HeII_IonizationRate = 0;
#else
      GF.RT_HI_heating_rate[k]      = 0;
      GF.RT_HeI_heating_rate[k]     = 0;
      GF.RT_HeII_heating_rate[k]    = 0;
      GF.RT_HI_ionization_rate[k]   = 0;
      GF.RT_HeI_ionization_rate[k]  = 0;
      GF.RT_HeII_ionization_rate[k] = 0;
#endif
      GF.specific_heating_rate[k] = 0.0;
#endif /* GRACKLE_CHEMISTRY >= 1 */
    }
}

/* Copy the abundances evolved by solve_chemistry() back into SphP (non-eq. chemistry only). Must be called with the same target[]
 * as the preceding grackle_gather(), and before SphP species are modified. */
static void grackle_scatter_species(int n, const int *target)
{
#if (GRACKLE_CHEMISTRY >= 1)
  for(int k = 0; k < n; k++)
    {
      int i = target[k];

      /* Grackle assumes the non-metal portion of the gas always retains its original primordial abundance ratio.
       * https://arxiv.org/abs/2604.00100v1 (Appendix A); see make_consistent_g in solve_rate_cool_g.F */
#ifdef METALS
      double Metallicity = SphP[i].GasMetallicity;
#else
      double Metallicity = GRACKLE_TINY;
#endif
      /* H-bearing mass fraction before the solve; these have not yet been overwritten */
      double X_H = SphP[i].GrackleSpecies(GRACKLE_HI) + SphP[i].GrackleSpecies(GRACKLE_HII);
#if (GRACKLE_CHEMISTRY >= 2)
      X_H += SphP[i].GrackleSpecies(GRACKLE_H2I) + SphP[i].GrackleSpecies(GRACKLE_H2II) + SphP[i].GrackleSpecies(GRACKLE_HM);
#endif
#if (GRACKLE_CHEMISTRY >= 3)
      X_H += SphP[i].GrackleSpecies(GRACKLE_DI) + SphP[i].GrackleSpecies(GRACKLE_DII) + SphP[i].GrackleSpecies(GRACKLE_HDI);
#endif
      double Y_He = 1 - X_H - Metallicity;

      gr_float HHemassfrac = 1.0 - Metallicity;           /* X_H + Y_He */
      gr_float fh_correct  = (X_H / HHemassfrac) / 0.76;  /* X/(X+Y) over the value Grackle assumes (0.76) */
      gr_float fhe_correct = (Y_He / HHemassfrac) / 0.24; /* Y/(X+Y) over the value Grackle assumes (0.24) */

      double rho = GF.density[k];

      SphP[i].GrackleSpecies(GRACKLE_HI)    = GF.HI_density[k] / rho * fh_correct;
      SphP[i].GrackleSpecies(GRACKLE_HII)   = GF.HII_density[k] / rho * fh_correct;
      SphP[i].GrackleSpecies(GRACKLE_HeI)   = GF.HeI_density[k] / rho * fhe_correct;
      SphP[i].GrackleSpecies(GRACKLE_HeII)  = GF.HeII_density[k] / rho * fhe_correct;
      SphP[i].GrackleSpecies(GRACKLE_HeIII) = GF.HeIII_density[k] / rho * fhe_correct;
      sync_conserved_from_primitive(i, GRACKLE_HI);
      sync_conserved_from_primitive(i, GRACKLE_HII);
      sync_conserved_from_primitive(i, GRACKLE_HeI);
      sync_conserved_from_primitive(i, GRACKLE_HeII);
      sync_conserved_from_primitive(i, GRACKLE_HeIII);

#if (GRACKLE_CHEMISTRY >= 2)
      SphP[i].GrackleSpecies(GRACKLE_H2I)  = GF.H2I_density[k] / rho * fh_correct;
      SphP[i].GrackleSpecies(GRACKLE_H2II) = GF.H2II_density[k] / rho * fh_correct;
      SphP[i].GrackleSpecies(GRACKLE_HM)   = GF.HM_density[k] / rho * fh_correct;
      sync_conserved_from_primitive(i, GRACKLE_H2I);
      sync_conserved_from_primitive(i, GRACKLE_H2II);
      sync_conserved_from_primitive(i, GRACKLE_HM);
#endif
#if (GRACKLE_CHEMISTRY >= 3)
      SphP[i].GrackleSpecies(GRACKLE_DI)  = GF.DI_density[k] / rho * fh_correct;
      SphP[i].GrackleSpecies(GRACKLE_DII) = GF.DII_density[k] / rho * fh_correct;
      SphP[i].GrackleSpecies(GRACKLE_HDI) = GF.HDI_density[k] / rho * fh_correct;
      sync_conserved_from_primitive(i, GRACKLE_DI);
      sync_conserved_from_primitive(i, GRACKLE_DII);
      sync_conserved_from_primitive(i, GRACKLE_HDI);
#endif
    }
#endif /* GRACKLE_CHEMISTRY >= 1 */
}

/* 'mode' -- tells the routine what to do
 *
 *     0 == solve chemistry and assign new abundances
 *     1 == calculate and return cooling time
 *     2 == calculate and return temperature
 *     3 == calculate and return pressure
 *     4 == calculate and return gamma (only valid when GRACKLE_CHEMISTRY>0, also it makes sense only when you have molecules)
 *
 * Single-cell wrapper; see GrackleCoolActive() for the batched version.
 */
double CallGrackle(double u_old, double rho, double dt, int target, int mode)
{
  gr_float returnval = 0.0;

  grackle_gather(1, &target, &u_old, &rho, mode == 0);

  switch(mode)
    {
      case 0: /* returns the new internal energy; in non-eq run it evolves the abundances */
        if(solve_chemistry(&All.GrackleUnits, &GF, dt) == 0)
          terminate("GRACKLE: Error in solve_chemistry.\n");
        grackle_scatter_species(1, &target);
        returnval = GF.internal_energy[0];
        break;

      case 1: /* returns the instantaneous cooling time */
        {
          gr_float cooling_time;
          if(calculate_cooling_time(&All.GrackleUnits, &GF, &cooling_time) == 0)
            terminate("GRACKLE: Error in calculate_cooling_time.\n");
          returnval = cooling_time;
          break;
        }
      case 2: /* returns the gas temperature */
        {
          gr_float temperature;
          if(calculate_temperature(&All.GrackleUnits, &GF, &temperature) == 0)
            terminate("GRACKLE: Error in calculate_temperature.\n");
          returnval = temperature;
          break;
        }
      case 3: /* returns the gas pressure */
        {
          gr_float pressure;
          if(calculate_pressure(&All.GrackleUnits, &GF, &pressure) == 0)
            terminate("GRACKLE: Error in calculate_pressure.\n");
          returnval = pressure;
          break;
        }
      case 4: /* returns gamma (effective adiabatic index) - useful only if molecules (H2,HD) are included */
        {
          gr_float gamma;
          if(calculate_gamma(&All.GrackleUnits, &GF, &gamma) == 0)
            terminate("GRACKLE: Error in calculate_gamma.\n");
          returnval = gamma;
          break;
        }
    }

  return returnval;
}

/*! \brief Cooling times for n cells at given thermal energies and proper densities, in one Grackle call.
 *
 *  Batched equivalent of calling GetCoolingTime() for each cell: tcool[k] is 0 for heating, otherwise Grackle's cooling time
 *  (negative = cooling; same convention as GetCoolingTime()). The radiation rates in SphP are left untouched.
 */
void GrackleCoolingTimes(int n, const int *cells, const double *u, const double *rho, double *tcool)
{
  if(n == 0)
    return;

  grackle_gather(n, cells, u, rho, 0);

  gr_float *t = mymalloc("grackle_tcool", n * sizeof(gr_float));
  if(calculate_cooling_time(&All.GrackleUnits, &GF, t) == 0)
    terminate("GRACKLE: Error in calculate_cooling_time.\n");

  for(int k = 0; k < n; k++)
    tcool[k] = (t[k] >= 0) ? 0.0 : t[k];

  myfree(t);
}

/*! \brief Apply Grackle cooling to a list of gas cells, batched per time bin.
 *
 *  Replaces the per-cell chain cool_cell() -> DoCooling() -> CallGrackle(). solve_chemistry() takes one scalar dt, and the cooling
 *  step of a cell is set by its hydro time bin, so the cells are grouped by TimeBinHydro and each non-empty bin is passed to
 *  Grackle as a single array. Per cell, the update of Utherm/Energy/pressure is identical to cool_cell().
 *
 *  Cells with Mass == 0 (swallowed/eliminated) are skipped. The thermal energy floor All.MinEgySpec is applied to the Grackle input
 *  exactly as in cool_cell().
 *
 *  \param[in] ncells number of entries in list
 *  \param[in] list   indices of the gas cells to cool
 */
void GrackleCoolCells(int ncells, const int *list)
{
  /* counting sort of the cells by hydro time bin */
  int count[TIMEBINS + 1];
  memset(count, 0, sizeof(count));

  for(int k = 0; k < ncells; k++)
    {
      int i = list[k];
      if(i >= 0 && P[i].Mass != 0)
        count[P[i].TimeBinHydro + 1]++;
    }
  for(int b = 0; b < TIMEBINS; b++)
    count[b + 1] += count[b];

  int nsel = count[TIMEBINS];
  if(nsel == 0)
    return;

  int *cell   = mymalloc("grackle_cell", nsel * sizeof(int));
  double *u   = mymalloc("grackle_u", nsel * sizeof(double));
  double *rho = mymalloc("grackle_rho", nsel * sizeof(double));
  int *fill   = mymalloc("grackle_fill", (TIMEBINS + 1) * sizeof(int));
  memcpy(fill, count, (TIMEBINS + 1) * sizeof(int));

  for(int k = 0; k < ncells; k++)
    {
      int i = list[k];
      if(i >= 0 && P[i].Mass != 0)
        cell[fill[P[i].TimeBinHydro]++] = i;
    }
  myfree(fill);

  for(int bin = 0; bin < TIMEBINS; bin++)
    {
      int n = count[bin + 1] - count[bin];
      if(n == 0)
        continue;

      int *c = cell + count[bin];
      double *ub = u + count[bin], *rb = rho + count[bin];

      double dt    = (bin ? (((integertime)1) << bin) : 0) * All.Timebase_interval;
      double dtime = All.cf_atime * dt / All.cf_time_hubble_a;

      for(int k = 0; k < n; k++)
        {
          ub[k] = dmax(All.MinEgySpec, SphP[c[k]].Utherm);
          rb[k] = SphP[c[k]].Density * All.cf_a3inv;
        }

      grackle_gather(n, c, ub, rb, 1);

#ifdef COOLING_TIMESTEP_LIMITER
      /* evaluated from the pre-cooling state */
      gr_float *tcool = mymalloc("grackle_tcool", n * sizeof(gr_float));
      if(calculate_cooling_time(&All.GrackleUnits, &GF, tcool) == 0)
        terminate("GRACKLE: Error in calculate_cooling_time.\n");
      for(int k = 0; k < n; k++)
        {
          /* Grackle's cooling time is in its code time_units (negative = cooling) */
          double t                   = (tcool[k] >= 0) ? 0.0 : tcool[k];
          SphP[c[k]].CoolingTimestep = (t != 0) ? All.CoolingTimestepFactor * fabs(t) : 0;
        }
      myfree(tcool);
#endif

      /* dt == 0 (very first step): nothing to integrate, Grackle would complain */
      if(dtime > 0)
        {
          if(solve_chemistry(&All.GrackleUnits, &GF, dtime) == 0)
            terminate("GRACKLE: Error in solve_chemistry.\n");
          grackle_scatter_species(n, c);
        }

      for(int k = 0; k < n; k++)
        {
          int i       = c[k];
          double unew = (dtime > 0) ? GF.internal_energy[k] : ub[k];

          if(unew < 0)
            terminate("invalid temperature: Thistask=%d i=%d unew=%g\n", ThisTask, i, unew);

          double du = unew - SphP[i].Utherm;

          if(unew < All.MinEgySpec)
            du = All.MinEgySpec - SphP[i].Utherm;

          SphP[i].Utherm += du;
          SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;

#ifdef OUTPUT_COOLHEAT
          if(dtime > 0)
            SphP[i].CoolHeat = du * P[i].Mass / dtime;
#endif
          set_pressure_of_cell(i);
        }
    }

  myfree(rho);
  myfree(u);
  myfree(cell);
}

/*! \brief Apply Grackle cooling to all active gas cells; see GrackleCoolCells(). */
void GrackleCoolActive(void)
{
  GrackleCoolCells(TimeBinsHydro.NActiveParticles, TimeBinsHydro.ActiveParticleList);
}

/* Function that initialises Grackle */
void InitGrackle(void)
{
  // Check the consistency
  if(gr_check_consistency() != GR_SUCCESS)
    {
      terminate("GRACKLE: Error in gr_check_consistency.\n");
    }

  int grackle_verbose = 1;
  // Enable output
  if(ThisTask == 0)
    grackle_verbose = 1;

  // First, set up the units system.
  // These are conversions from code units to cgs.
  /* 1 for cosmo, 0 if not; you can put it equal to All.ComovingIntegrationOn in order to set it automatically.
   * Beware though that for a non-cosmo sim a_value has to be kept constant.
   * IMPORTANT: at the moment the density is already converted into proper frame before calling the cooling routine,
   * so it's better to set this = 0 and only change a_value at every timestep.
   */

  double h = All.HubbleParam;

  All.GrackleUnits.comoving_coordinates = 0;  // All.ComovingIntegrationOn; // 1 if cosmological sim, 0 if not
  All.GrackleUnits.density_units        = All.UnitDensity_in_cgs / h / h;
  All.GrackleUnits.length_units         = All.UnitLength_in_cm / h;
  All.GrackleUnits.time_units           = All.UnitTime_in_s / h;
  All.GrackleUnits.a_units              = 1.0;  // units for the expansion factor; NOTE: Should be 1 always
  set_velocity_units(&All.GrackleUnits);

  if(ThisTask == 0)
    if(fabs(All.GrackleUnits.velocity_units - All.UnitVelocity_in_cm_per_s) > 1e-9)
      {
        terminate("GRACKLE: Velocity units doesn't match; Check the units!");
      }

  // Set initial expansion factor (for internal units).
  // Set expansion factor to 1 for non-cosmological simulation.
  /* Set initial expansion factor (for internal units).
   * Set expansion factor to 1 (or to the scale factor that corresponds to the right redshift) for non-cosmological simulation.
   * This is the only parameter that we need to change at every timestep in the GrackleUnits struct.
   * All the others are constant throughout the simulation, since the density is already converted into proper frame before calling the
   * cooling routine.
   */
  double a_value = 1.0;
  if(All.ComovingIntegrationOn)
    a_value = All.TimeBegin;

  All.GrackleUnits.a_value = a_value;

  // Second, create a chemistry object for parameters and rate data.
  chemistry_data* my_grackle_data;
  my_grackle_data = malloc(sizeof(chemistry_data));

  if(set_default_chemistry_parameters(my_grackle_data) == 0)
    {
      terminate("\nGRACKLE: Error in set_default_chemistry_parameters.\n");
    }

  /* Third, fill this chemistry object to set parameter values for chemistry & cooling. */

  /* main flags */

  /* Flag to activate the grackle machinery */
  my_grackle_data->use_grackle = 1;
  /* Path to the data file containing the metal cooling and UV background tables (for non-eq mode) and metal and primordial
   * cooling/heating (for equilibrium mode) */
  my_grackle_data->grackle_data_file = All.GrackleDataFile;
  /* Flag to include radiative cooling and actually update the thermal energy during the
   * chemistry solver. If off, the chemistry species will still be updated. The most
   * common reason to set this to off is to iterate the chemistry network to an equilibrium state. Default: 1.
   */
  my_grackle_data->with_radiative_cooling = 1;

  /* The ratio of specific heats for an ideal gas. A direct calculation for the molecular component is used if primordial_chemistry
   * > 1. Default: 5/3. */
  my_grackle_data->Gamma = GAMMA; /* our eos set in Config.sh */
  /* Flag to control which primordial chemistry network is used (set by Config.sh) */

#ifndef GRACKLE_CHEMISTRY
  my_grackle_data->primordial_chemistry = 0; /* if nothing is set assume fully tabulated cooling */
#else
  my_grackle_data->primordial_chemistry = GRACKLE_CHEMISTRY;
#endif

#ifdef METALS
  /* Flag to enable metal cooling using the Cloudy tables. If enabled, the cooling table to be used must be specified with the
   * grackle_data_file parameter. Default: 0. */
  my_grackle_data->metal_cooling = 1;
  /* Flag to enable H2 formation on dust grains, dust cooling, and dust-gas heat transfer follow Omukai (2000). This assumes that the
   * dust to gas ratio scales with the metallicity. Default: 0. */
  my_grackle_data->h2_on_dust = 0;
  /* Flag to enable a spatially uniform heating term approximating photo-electric heating from dust from Tasker & Bryan (2008).
   * Default: 0. If photoelectric_heating enabled, photoelectric_heating_rate is the heating rate in units of erg cm-3 s-1.
   * Default: 8.5e-26. This is not adjusted to local background. (Caution: this tends to heat gas even at extremely high densities to
   * ~3000 K, when it should be entirely self-shielding) Another comment: this heats also all the low-density gas to 10^4 K, so it's
   * better not to use it.
   */
  my_grackle_data->photoelectric_heating      = 0; /* read above but DO NOT USE */
  my_grackle_data->photoelectric_heating_rate = 8.5e-26;
#else
  my_grackle_data->metal_cooling              = 0;
  my_grackle_data->h2_on_dust                 = 0;
  my_grackle_data->photoelectric_heating      = 0;
  my_grackle_data->photoelectric_heating_rate = 8.5e-26;
#endif

  /* Flag to control which three-body H2 formation rate is used.
   *    0: Abel, Bryan & Norman (2002),
   *    1: Palla, Salpeter & Stahler (1983),
   *    2: Cohen & Westberg (1983),
   *    3: Flower & Harris (2007),
   *    4: Glover (2008).
   *    These are discussed in Turk et. al. (2011). Default: 0.
   */
  my_grackle_data->three_body_rate = 0;

  /* Flag to enable an effective CMB temperature floor.
   * This is implemented by subtracting the value of the cooling rate at TCMB from the total METAL cooling rate. Default: 1.
   * Beware! You could still have Tgas<TCBM because it imposes a temperature floor only for the metal cooling.
   */
  my_grackle_data->cmb_temperature_floor = 1;

  /* Flag to enable a UV background.
   * If enabled, the cooling table to be used must be specified with the grackle_data_file parameter. Default: 0.
   */
  my_grackle_data->UVbackground = 1;
  /* The following flags are related to the UVB, but they are automatically set to the right values, so do not need to use. These
   * numbers are the correct ones for FG2011 UVB. my_grackle_data->UVbackground_redshift_on       = 10.6;
   * my_grackle_data->UVbackground_redshift_off      = 0;
   * my_grackle_data->UVbackground_redshift_fullon   = 10.6;
   * my_grackle_data->UVbackground_redshift_drop     = 0;
   */

  /* Flag to enable Compton heating from an X-ray background following Madau & Efstathiou (1999). Default: 0.
   * Beware: this flag is just broken, it causes a runaway ionisation and heating as soon as the UVB kicks in
   * and you start having some electrons. DO NOT use it.
   */
  my_grackle_data->Compton_xray_heating = 0; /* read above but DO NOT USE */

  /* Flag to enable H2 collision-induced emission cooling from Ripamonti & Abel (2004). Default: 0. */
  my_grackle_data->cie_cooling = 0;
  /* Flag to enable H2 cooling attenuation from Ripamonti & Abel (2004). Default: 0. */
  my_grackle_data->h2_optical_depth_approximation = 0;

  /* Rad Intensity of a constant Lyman-Werner H2 photo-dissociating radiation field,
   * in units of 10-21 erg s-1 cm-2 Hz-1 sr-1. Default: 0.
   */
  my_grackle_data->LWbackground_intensity = 0;
  /* Flag to enable suppression of Lyman-Werner flux due to Lyman-series absorption
   * (giving a sawtooth pattern), taken from Haiman & Abel, & Rees (2000). Default: 0.
   */
  my_grackle_data->LWbackground_sawtooth_suppression = 0;

  /* For both UV bkgd and RT; options for length scale:
   *     1: Sobolev-like (from WG11)
   *     2: array of lengths
   *     3: local Jeans length
   * Default: 0.
   */
  my_grackle_data->H2_self_shielding = 3;

  /* Flag for self-shielding from UV bkgd. Default: 0. */
  my_grackle_data->self_shielding_method = 0;

#ifdef STAR_RADIATION_ACTIVE
  /* flag to include RT */
  my_grackle_data->use_radiative_transfer      = 1;
  my_grackle_data->use_volumetric_heating_rate = 1;
#else
  my_grackle_data->use_radiative_transfer      = 0;
  my_grackle_data->use_volumetric_heating_rate = 0;
#endif

  /* Finally, initialize the chemistry object. This has to be the last step of the initialisation. */
  if(initialize_chemistry_data(&All.GrackleUnits) == 0)
    {
      terminate("Error in initialize_chemistry_data.\n");
    }

  if(ThisTask == 0)
    printf("GRACKLE: Grackle Initialized\n");
}

double compute_mu(int i)
{
  /* Level 1: atomic H and He */
#if GRACKLE_CHEMISTRY >= 1
  double XHI = SphP[i].GrackleSpeciesConserved(GRACKLE_HI) / P[i].Mass;
  double XHII = SphP[i].GrackleSpeciesConserved(GRACKLE_HII) / P[i].Mass;
  double XHeI = SphP[i].GrackleSpeciesConserved(GRACKLE_HeI) / P[i].Mass;
  double XHeII = SphP[i].GrackleSpeciesConserved(GRACKLE_HeII) / P[i].Mass;
  double XHeIII = SphP[i].GrackleSpeciesConserved(GRACKLE_HeIII) / P[i].Mass;
#else
  /* Fall back to fully neutral cosmic abundances */
  double XHI = HYDROGEN_MASSFRAC;
  double XHII = 0.0;
  double XHeI = 1.0 - HYDROGEN_MASSFRAC;
  double XHeII = 0.0;
  double XHeIII = 0.0;
#endif

  /* Level 2: molecular H and H- */
#if GRACKLE_CHEMISTRY >= 2
  double XH2I = SphP[i].GrackleSpeciesConserved(GRACKLE_H2I) / P[i].Mass;
  double XH2II = SphP[i].GrackleSpeciesConserved(GRACKLE_H2II) / P[i].Mass;
  double XHM = SphP[i].GrackleSpeciesConserved(GRACKLE_HM) / P[i].Mass;
#else
  double XH2I = 0.0;
  double XH2II = 0.0;
  double XHM = 0.0;
#endif

  /* Level 3: deuterium species */
#if GRACKLE_CHEMISTRY >= 3
  double XDI = SphP[i].GrackleSpeciesConserved(GRACKLE_DI) / P[i].Mass;
  double XDII = SphP[i].GrackleSpeciesConserved(GRACKLE_DII) / P[i].Mass;
  double XHDI = SphP[i].GrackleSpeciesConserved(GRACKLE_HDI) / P[i].Mass;
#else
  double XDI = 0.0;
  double XDII = 0.0;
  double XHDI = 0.0;
#endif

  double Xe = XHII + XHeII / 4.0 + XHeIII / 2.0 + XH2II / 2.0 - XHM + XDII / 2.0;

  /* Assemble grouped mass fractions */
  double XH = XHI + XHII + XHM; /* m_H */
  double XH2 = XH2I + XH2II; /* 2 m_H */
  double XD = XDI + XDII; /* 2 m_H */
  double XHD = XHDI; /* 3 m_H */
  double XHe = XHeI + XHeII + XHeIII; /* 4 m_H */

/* Metals, approximated as 16 m_H */
#ifdef METALS
  double Z = SphP[i].GasMetals / P[i].Mass;
#else
  double Z = 0; 
#endif

  /* mu = 1 / sum_s (X_s / A_s), where A_s is the atomic mass in units of m_H */
  return 1.0 / (Xe + XH + XH2 / 2.0 + XD / 2.0 + XHD / 3.0 + XHe / 4.0 + Z / 16.0);
}