#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

/* 'mode' -- tells the routine what to do
 *
 *     0 == solve chemistry and assign new abundances
 *     1 == calculate and return cooling time
 *     2 == calculate and return temperature
 *     3 == calculate and return pressure
 *     4 == calculate and return gamma (only valid when GRACKLE_CHEMISTRY>0, also it makes sense only when you have molecules)
 */
double CallGrackle(double u_old, double rho, double dt, int target, int mode)
{
  gr_initialize_field_data(&All.GrackleFieldData);

  /* this is the variable that is returned */
  gr_float returnval = 0.0;

  int field_size = 1;

  /* set the right scale factor */
  All.GrackleUnits.a_value = All.cf_atime;

  /* grid setup */
  All.GrackleFieldData.grid_rank      = 3;
  All.GrackleFieldData.grid_dimension = malloc(3 * sizeof(int));
  All.GrackleFieldData.grid_start     = malloc(3 * sizeof(int));
  All.GrackleFieldData.grid_end       = malloc(3 * sizeof(int));

  for(int i = 0; i < 3; i++)
    {
      All.GrackleFieldData.grid_dimension[i] = 1;
      All.GrackleFieldData.grid_start[i]     = 0;
      All.GrackleFieldData.grid_end[i]       = 0;
    }

  All.GrackleFieldData.grid_dimension[0] = field_size;
  All.GrackleFieldData.grid_end[0]       = field_size - 1;

  /* basic allocations */
  All.GrackleFieldData.x_velocity      = malloc(sizeof(gr_float));
  All.GrackleFieldData.y_velocity      = malloc(sizeof(gr_float));
  All.GrackleFieldData.z_velocity      = malloc(sizeof(gr_float));
  All.GrackleFieldData.density         = malloc(sizeof(gr_float));
  All.GrackleFieldData.internal_energy = malloc(sizeof(gr_float));
  All.GrackleFieldData.metal_density   = malloc(sizeof(gr_float));

#if (GRACKLE_CHEMISTRY >= 1)
  All.GrackleFieldData.e_density     = malloc(sizeof(gr_float));
  All.GrackleFieldData.HI_density    = malloc(sizeof(gr_float));
  All.GrackleFieldData.HII_density   = malloc(sizeof(gr_float));
  All.GrackleFieldData.HeI_density   = malloc(sizeof(gr_float));
  All.GrackleFieldData.HeII_density  = malloc(sizeof(gr_float));
  All.GrackleFieldData.HeIII_density = malloc(sizeof(gr_float));
  All.GrackleFieldData.H2I_density   = malloc(sizeof(gr_float));
  All.GrackleFieldData.H2II_density  = malloc(sizeof(gr_float));
  All.GrackleFieldData.HM_density    = malloc(sizeof(gr_float));
  All.GrackleFieldData.DI_density    = malloc(sizeof(gr_float));
  All.GrackleFieldData.DII_density   = malloc(sizeof(gr_float));
  All.GrackleFieldData.HDI_density   = malloc(sizeof(gr_float));

  // Volumetric heating rate (provide in units [erg s^-1 cm^-3])
  All.GrackleFieldData.volumetric_heating_rate = malloc(sizeof(gr_float));
  // H2 dissociation rate from radiative transfer calculations (provide in units of [1/time_units])
  All.GrackleFieldData.RT_H2_dissociation_rate = malloc(sizeof(gr_float));
  // Heating rate from radiative transfer calculations (provide in units [erg s^-1 cm^-3] / n)
  All.GrackleFieldData.RT_HI_heating_rate = malloc(sizeof(gr_float));
  // Heating rate from radiative transfer calculations (provide in units [erg s^-1 cm^-3] / n)
  All.GrackleFieldData.RT_HeI_heating_rate = malloc(sizeof(gr_float));
  // Heating rate from radiative transfer calculations (provide in units [erg s^-1 cm^-3] / n)
  All.GrackleFieldData.RT_HeII_heating_rate = malloc(sizeof(gr_float));
  // HI ionization rate from radiative transfer calculations (provide in units of [1/time_units])
  All.GrackleFieldData.RT_HI_ionization_rate = malloc(sizeof(gr_float));
  // HeI ionization rate from radiative transfer calculations (provide in units of [1/time_units])
  All.GrackleFieldData.RT_HeI_ionization_rate = malloc(sizeof(gr_float));
  // HeII ionization rate from radiative transfer calculations (provide in units of [1/time_units])
  All.GrackleFieldData.RT_HeII_ionization_rate = malloc(sizeof(gr_float));

  // specific heating rate (provide in units [egs s^-1 g^-1]
  All.GrackleFieldData.specific_heating_rate = malloc(sizeof(gr_float));
#endif

  /* basic values */
  *All.GrackleFieldData.x_velocity      = P[target].Vel[0];
  *All.GrackleFieldData.y_velocity      = P[target].Vel[1];
  *All.GrackleFieldData.z_velocity      = P[target].Vel[2];
  *All.GrackleFieldData.density         = rho;
  *All.GrackleFieldData.internal_energy = u_old;

#ifdef METALS 
  double Metallicity = SphP[target].GasMetallicity;
#else
  double Metallicity = GRACKLE_TINY;
#endif

  *All.GrackleFieldData.metal_density = Metallicity * *All.GrackleFieldData.density;

  /* non-eq. chemistry values */
#if (GRACKLE_CHEMISTRY >= 1)

  // Let's get the abundances before we call grackle
  double X_H = 0, Y_He = 0;

  /* electron density */
  //*All.GrackleFieldData.e_density = SphP[target].Ne * *All.GrackleFieldData.density;

  /* H and He species */
  *All.GrackleFieldData.HI_density = SphP[target].GrackleSpecies(GRACKLE_HI) * *All.GrackleFieldData.density;
  *All.GrackleFieldData.HII_density = SphP[target].GrackleSpecies(GRACKLE_HII) * *All.GrackleFieldData.density;
  *All.GrackleFieldData.HeI_density = SphP[target].GrackleSpecies(GRACKLE_HeI) * *All.GrackleFieldData.density;
  *All.GrackleFieldData.HeII_density = SphP[target].GrackleSpecies(GRACKLE_HeII) * *All.GrackleFieldData.density;
  *All.GrackleFieldData.HeIII_density = SphP[target].GrackleSpecies(GRACKLE_HeIII) * *All.GrackleFieldData.density;

  X_H += SphP[target].GrackleSpecies(GRACKLE_HI) + SphP[target].GrackleSpecies(GRACKLE_HII);

  /* molecular H species */
#if (GRACKLE_CHEMISTRY >= 2)
  *All.GrackleFieldData.H2I_density = SphP[target].GrackleSpecies(GRACKLE_H2I) * *All.GrackleFieldData.density;
  *All.GrackleFieldData.H2II_density = SphP[target].GrackleSpecies(GRACKLE_H2II) * *All.GrackleFieldData.density;
  *All.GrackleFieldData.HM_density = SphP[target].GrackleSpecies(GRACKLE_HM) * *All.GrackleFieldData.density;

  X_H += SphP[target].GrackleSpecies(GRACKLE_H2I) + SphP[target].GrackleSpecies(GRACKLE_H2II) + SphP[target].GrackleSpecies(GRACKLE_HM);

#else
  *All.GrackleFieldData.H2I_density = GRACKLE_TINY * *All.GrackleFieldData.density;
  *All.GrackleFieldData.H2II_density = GRACKLE_TINY * *All.GrackleFieldData.density;
  *All.GrackleFieldData.HM_density = GRACKLE_TINY * *All.GrackleFieldData.density;
#endif

  /* deuterium species */
#if (GRACKLE_CHEMISTRY >= 3)
  *All.GrackleFieldData.DI_density = SphP[target].GrackleSpecies(GRACKLE_DI) * *All.GrackleFieldData.density;
  *All.GrackleFieldData.DII_density = SphP[target].GrackleSpecies(GRACKLE_DII) * *All.GrackleFieldData.density;
  *All.GrackleFieldData.HDI_density = SphP[target].GrackleSpecies(GRACKLE_HDI) * *All.GrackleFieldData.density;

  X_H += SphP[target].GrackleSpecies(GRACKLE_DI) + SphP[target].GrackleSpecies(GRACKLE_DII) + SphP[target].GrackleSpecies(GRACKLE_HDI);

#else
  *All.GrackleFieldData.DI_density = GRACKLE_TINY * *All.GrackleFieldData.density;
  *All.GrackleFieldData.DII_density = GRACKLE_TINY * *All.GrackleFieldData.density;
  *All.GrackleFieldData.HDI_density = GRACKLE_TINY * *All.GrackleFieldData.density;
#endif

  Y_He = 1 - X_H - Metallicity;

  double e_density = 0;

  e_density += *All.GrackleFieldData.HII_density;
  e_density += *All.GrackleFieldData.HeII_density / 4.0;
  e_density += *All.GrackleFieldData.HeIII_density / 2.0;

#if (GRACKLE_CHEMISTRY >= 2)
  e_density += *All.GrackleFieldData.H2II_density / 2.0;
  e_density -= *All.GrackleFieldData.HM_density;
#endif

#if (GRACKLE_CHEMISTRY >= 3)
  e_density += *All.GrackleFieldData.DII_density / 2.0;
#endif

  *All.GrackleFieldData.e_density = e_density;

  /* Radiation */
#ifdef PHOTOELECTRIC_HEATING
  *All.GrackleFieldData.volumetric_heating_rate = (gr_float)(SphP[target].PE_VolHeatingRate);
  
  SphP[target].PE_VolHeatingRate = 0;
#else
  *All.GrackleFieldData.volumetric_heating_rate = 0;
#endif

#ifdef DISSOCIATION
  *All.GrackleFieldData.RT_H2_dissociation_rate = (gr_float)(SphP[target].H2_DissociationRate);
  
  SphP[target].H2_DissociationRate = 0;
#else
  *All.GrackleFieldData.RT_H2_dissociation_rate = 0;
#endif

#ifdef PHOTOIONIZATION
  *All.GrackleFieldData.RT_HI_heating_rate = (gr_float)(SphP[target].HI_HeatingRate);
  *All.GrackleFieldData.RT_HeI_heating_rate = (gr_float)(SphP[target].HeI_HeatingRate);
  *All.GrackleFieldData.RT_HeII_heating_rate = (gr_float)(SphP[target].HeII_HeatingRate);
  *All.GrackleFieldData.RT_HI_ionization_rate = (gr_float)(SphP[target].HI_IonizationRate);
  *All.GrackleFieldData.RT_HeI_ionization_rate = (gr_float)(SphP[target].HeI_IonizationRate);
  *All.GrackleFieldData.RT_HeII_ionization_rate = (gr_float)(SphP[target].HeII_IonizationRate);

  SphP[target].HI_HeatingRate = SphP[target].HeI_HeatingRate = SphP[target].HeII_HeatingRate 
  = SphP[target].HI_IonizationRate = SphP[target].HeI_IonizationRate = SphP[target].HeII_IonizationRate = 0;
#else
  *All.GrackleFieldData.RT_HI_heating_rate = 0;
  *All.GrackleFieldData.RT_HeI_heating_rate = 0;
  *All.GrackleFieldData.RT_HeII_heating_rate = 0;
  *All.GrackleFieldData.RT_HI_ionization_rate = 0;
  *All.GrackleFieldData.RT_HeI_ionization_rate = 0;
  *All.GrackleFieldData.RT_HeII_ionization_rate = 0;
#endif

  *All.GrackleFieldData.specific_heating_rate = 0.0;

#endif /* GRACKLE_CHEMISTRY >= 1 */

  /* call to the Grackle functions
   * remember: Grackle3 does not distinguish between non-equilibrium and tabulated version
   */
  switch(mode)
    {
      case 0: /* returns the new internal energy; in non-eq run it evolves the abundances */
        {
          if(solve_chemistry(&All.GrackleUnits, &All.GrackleFieldData, dt) == 0)
            {
              terminate("GRACKLE: Error in solve_chemistry.\n");
            }

          // Grackle assumes non-metal portion of the gas always retains its original primordial abundances ratio.
          // https://arxiv.org/abs/2604.00100v1 (Appendix A)
          // check make_consistent_g in solve_rate_cool_g.F
          gr_float HHemassfrac = 1.0 - Metallicity;  // X_H+Y_He
          // For H: X/(X+Y) divided value used in grackle 0.76
          gr_float fh_correct = (X_H / HHemassfrac) / HYDROGEN_MASSFRAC;
          // For He: Y/(X+Y) divided value in grackle 0.24
          gr_float fhe_correct = (Y_He / HHemassfrac) / (1 - HYDROGEN_MASSFRAC);

          /* if non-eq chemistry assign abundances back */
#if (GRACKLE_CHEMISTRY >= 1)
          // We balance the charges to get the proper Ne once all the abundances are assigned.
          // SphP[target].Ne                            = *All.GrackleFieldData.e_density / *All.GrackleFieldData.density;
          SphP[target].GrackleSpecies(GRACKLE_HI) = *All.GrackleFieldData.HI_density / *All.GrackleFieldData.density;
          SphP[target].GrackleSpecies(GRACKLE_HII) = *All.GrackleFieldData.HII_density / *All.GrackleFieldData.density;
          SphP[target].GrackleSpecies(GRACKLE_HeI) = *All.GrackleFieldData.HeI_density / *All.GrackleFieldData.density;
          SphP[target].GrackleSpecies(GRACKLE_HeII) = *All.GrackleFieldData.HeII_density / *All.GrackleFieldData.density;
          SphP[target].GrackleSpecies(GRACKLE_HeIII) = *All.GrackleFieldData.HeIII_density / *All.GrackleFieldData.density;

          SphP[target].GrackleSpecies(GRACKLE_HI) *= fh_correct;
          SphP[target].GrackleSpecies(GRACKLE_HII) *= fh_correct;
          SphP[target].GrackleSpecies(GRACKLE_HeI) *= fhe_correct;
          SphP[target].GrackleSpecies(GRACKLE_HeII) *= fhe_correct;
          SphP[target].GrackleSpecies(GRACKLE_HeIII) *= fhe_correct;

          //SphP[target].Ne = SphP[target].GrackleSpecies(GRACKLE_HII) + SphP[target].GrackleSpecies(GRACKLE_HeII) / 4. 
          //+ SphP[target].GrackleSpecies(GRACKLE_HeIII) / 2.;

          sync_conserved_from_primitive(target, GRACKLE_HI);
          sync_conserved_from_primitive(target, GRACKLE_HII);
          sync_conserved_from_primitive(target, GRACKLE_HeI);
          sync_conserved_from_primitive(target, GRACKLE_HeII);
          sync_conserved_from_primitive(target, GRACKLE_HeIII);
#endif

#if (GRACKLE_CHEMISTRY >= 2)
          SphP[target].GrackleSpecies(GRACKLE_H2I) = *All.GrackleFieldData.H2I_density / *All.GrackleFieldData.density;
          SphP[target].GrackleSpecies(GRACKLE_H2II) = *All.GrackleFieldData.H2II_density / *All.GrackleFieldData.density;
          SphP[target].GrackleSpecies(GRACKLE_HM) = *All.GrackleFieldData.HM_density / *All.GrackleFieldData.density;

          SphP[target].GrackleSpecies(GRACKLE_H2I) *= fh_correct;
          SphP[target].GrackleSpecies(GRACKLE_H2II) *= fh_correct;
          SphP[target].GrackleSpecies(GRACKLE_HM) *= fh_correct;

          //SphP[target].Ne += SphP[target].GrackleSpecies(GRACKLE_H2II) / 2. - SphP[target].GrackleSpecies(GRACKLE_HM);

          sync_conserved_from_primitive(target, GRACKLE_H2I);
          sync_conserved_from_primitive(target, GRACKLE_H2II);
          sync_conserved_from_primitive(target, GRACKLE_HM);    
#endif

#if (GRACKLE_CHEMISTRY >= 3)
          SphP[target].GrackleSpecies(GRACKLE_DI) = *All.GrackleFieldData.DI_density / *All.GrackleFieldData.density;
          SphP[target].GrackleSpecies(GRACKLE_DII) = *All.GrackleFieldData.DII_density / *All.GrackleFieldData.density;
          SphP[target].GrackleSpecies(GRACKLE_HDI) = *All.GrackleFieldData.HDI_density / *All.GrackleFieldData.density;

          SphP[target].GrackleSpecies(GRACKLE_DI) *= fh_correct;
          SphP[target].GrackleSpecies(GRACKLE_DII) *= fh_correct;
          SphP[target].GrackleSpecies(GRACKLE_HDI) *= fh_correct;

          //SphP[target].Ne += SphP[target].GrackleSpecies(GRACKLE_DII) / 2.;

          sync_conserved_from_primitive(target, GRACKLE_DI);
          sync_conserved_from_primitive(target, GRACKLE_DII);
          sync_conserved_from_primitive(target, GRACKLE_HDI);
#endif

          returnval = *All.GrackleFieldData.internal_energy;
          break;
        }

      case 1: /* returns the instantaneous cooling time */
        {
          gr_float cooling_time;
          if(calculate_cooling_time(&All.GrackleUnits, &All.GrackleFieldData, &cooling_time) == 0)
            {
              terminate("GRACKLE: Error in calculate_cooling_time.\n");
            }
          returnval = cooling_time;
          break;
        }

      case 2: /* returns the gas temperature */
        {
          gr_float temperature;
          if(calculate_temperature(&All.GrackleUnits, &All.GrackleFieldData, &temperature) == 0)
            {
              terminate("GRACKLE: Error in calculate_temperature.\n");
            }
          returnval = temperature;
          break;
        }

      case 3: /* returns the gas pressure */
        {
          gr_float pressure;
          if(calculate_pressure(&All.GrackleUnits, &All.GrackleFieldData, &pressure) == 0)
            {
              terminate("GRACKLE: Error in calculate_pressure.\n");
            }
          returnval = pressure;
          break;
        }

      case 4: /* returns gamma (effective adiabatic index) - useful only if molecules (H2,HD) are included in the chemical network */
        {
          gr_float gamma;
          if(calculate_gamma(&All.GrackleUnits, &All.GrackleFieldData, &gamma) == 0)
            {
              terminate("GRACKLE: Error in calculate_gamma.\n");
            }
          returnval = gamma;
          break;
        }
    } /* end switch */

  /* free the memory */
  free(All.GrackleFieldData.grid_dimension);
  free(All.GrackleFieldData.grid_start);
  free(All.GrackleFieldData.grid_end);
  free(All.GrackleFieldData.x_velocity);
  free(All.GrackleFieldData.y_velocity);
  free(All.GrackleFieldData.z_velocity);
  free(All.GrackleFieldData.density);
  free(All.GrackleFieldData.internal_energy);
  free(All.GrackleFieldData.metal_density);

#if (GRACKLE_CHEMISTRY >= 1)
  free(All.GrackleFieldData.e_density);
  free(All.GrackleFieldData.HI_density);
  free(All.GrackleFieldData.HII_density);
  free(All.GrackleFieldData.HeI_density);
  free(All.GrackleFieldData.HeII_density);
  free(All.GrackleFieldData.HeIII_density);
  free(All.GrackleFieldData.H2I_density);
  free(All.GrackleFieldData.H2II_density);
  free(All.GrackleFieldData.HM_density);
  free(All.GrackleFieldData.DI_density);
  free(All.GrackleFieldData.DII_density);
  free(All.GrackleFieldData.HDI_density);
  free(All.GrackleFieldData.volumetric_heating_rate);
  free(All.GrackleFieldData.RT_H2_dissociation_rate);
  free(All.GrackleFieldData.RT_HI_heating_rate);
  free(All.GrackleFieldData.RT_HeI_heating_rate);
  free(All.GrackleFieldData.RT_HeII_heating_rate);
  free(All.GrackleFieldData.RT_HI_ionization_rate);
  free(All.GrackleFieldData.RT_HeI_ionization_rate);
  free(All.GrackleFieldData.RT_HeII_ionization_rate);

  free(All.GrackleFieldData.specific_heating_rate);
#endif

  return returnval;
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

  my_grackle_data->HydrogenFractionByMass = HYDROGEN_MASSFRAC;    

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
/* Metals, approximated as 16 m_H */
#ifdef METALS
  double Z = SphP[i].GasMetals / P[i].Mass;
#else
  double Z = 0; 
#endif

  /* Level 1: atomic H and He */
#if GRACKLE_CHEMISTRY >= 1
  double XHI = SphP[i].GrackleSpeciesConserved(GRACKLE_HI) / P[i].Mass;
  double XHII = SphP[i].GrackleSpeciesConserved(GRACKLE_HII) / P[i].Mass;
  double XHeI = SphP[i].GrackleSpeciesConserved(GRACKLE_HeI) / P[i].Mass;
  double XHeII = SphP[i].GrackleSpeciesConserved(GRACKLE_HeII) / P[i].Mass;
  double XHeIII = SphP[i].GrackleSpeciesConserved(GRACKLE_HeIII) / P[i].Mass;
#else

  /* Fall back to fully neutral cosmic abundances */
  double XHI = (1.0 - Z) * HYDROGEN_MASSFRAC;
  double XHII = 0.0;
  double XHeI = (1.0 - Z) * (1.0 - HYDROGEN_MASSFRAC);
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

  /* mu = 1 / sum_s (X_s / A_s), where A_s is the atomic mass in units of m_H */
  return 1.0 / (Xe + XH + XH2 / 2.0 + XD / 2.0 + XHD / 3.0 + XHe / 4.0 + Z / 16.0);
}