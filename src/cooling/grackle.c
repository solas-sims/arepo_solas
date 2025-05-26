#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef USE_GRACKLE
#include <grackle.h>

/* 'mode' -- tells the routine what to do
 *
 *     0 == solve chemistry and assign new abundances
 *     1 == calculate and return cooling time
 *     2 == calculate and return temperature
 *     3 == calculate and return pressure
 *     4 == calculate and return gamma (only valid when GRACKLE_CHEMISTRY>0, also it makes sense only when you have molecules)
 */
double CallGrackle(double u_old, double rho, double dt, double *ne_guess, int target, int mode)
{
    gr_initialize_field_data(&All.GrackleFieldData);
    
    /* this is the variable that is returned */
    gr_float returnval                                  = 0.0;

    int field_size                                      = 1;
    gr_float tiny                                       = 1.0e-20;    /* just a very low number good for initialisation -> you may want to use i.e. initial values from Galli & Palla */

    /* set the right scale factor */
    All.GrackleUnits.a_value                            = All.cf_atime;

    All.GrackleFieldData.grid_rank                      = 3;
    All.GrackleFieldData.grid_dimension                 = malloc(3*sizeof(int));
    All.GrackleFieldData.grid_start                     = malloc(3*sizeof(int));
    All.GrackleFieldData.grid_end                       = malloc(3*sizeof(int));
    int i;
    for (i = 0;i < 3;i++) {
        All.GrackleFieldData.grid_dimension[i]          = 1; /* the active dimension not including ghost zones; maybe this just means that the sim is 3D */
        All.GrackleFieldData.grid_start[i]              = 0;
        All.GrackleFieldData.grid_end[i]                = 0;
    }
    All.GrackleFieldData.grid_dimension[0]              = field_size;
    All.GrackleFieldData.grid_end[0]                    = field_size - 1;

    All.GrackleFieldData.x_velocity                     = malloc(sizeof(gr_float));
    All.GrackleFieldData.y_velocity                     = malloc(sizeof(gr_float));
    All.GrackleFieldData.z_velocity                     = malloc(sizeof(gr_float));
    All.GrackleFieldData.density                        = malloc(sizeof(gr_float));
    All.GrackleFieldData.internal_energy                = malloc(sizeof(gr_float));
    All.GrackleFieldData.metal_density                  = malloc(sizeof(gr_float));

    *All.GrackleFieldData.x_velocity                    = P[target].Vel[0];     /* currently not used */
    *All.GrackleFieldData.y_velocity                    = P[target].Vel[1];     /* currently not used */
    *All.GrackleFieldData.z_velocity                    = P[target].Vel[2];     /* currently not used */
    *All.GrackleFieldData.density                       = rho;
    *All.GrackleFieldData.internal_energy               = u_old;

#ifdef METALS //TODO:
    *All.GrackleFieldData.metal_density                 = *All.GrackleFieldData.density * SphP[target].Metals;
#else
    *All.GrackleFieldData.metal_density                 = *All.GrackleFieldData.density * SOLAR_ABUNDANCE;
#endif

    /* UNDEFINED MEMBERS OF THE STRUCT (USEFUL ONLY FOR UNUSED FLAGS)
     * beware! they should be created with malloc(sizeof(..))
     * also, most of them would be used only in non-eq run
     *
     * gr_float grid_dx;
     *
     * gr_float *volumetric_heating_rate;
     * gr_float *specific_heating_rate;
     *
     * gr_float *RT_heating_rate;
     * gr_float *RT_HI_ionization_rate;
     * gr_float *RT_HeI_ionization_rate;
     * gr_float *RT_HeII_ionization_rate;
     * gr_float *RT_H2_dissociation_rate;
     *
     * gr_float *H2_self_shielding_length;
     */


/* non-equilibrium chemistry */ /* TODO: */
#if (GRACKLE_CHEMISTRY >= 1)   /* atomic chemical network only */

    All.GrackleFieldData.e_density                      = malloc(sizeof(gr_float));
    All.GrackleFieldData.HI_density                     = malloc(sizeof(gr_float));
    All.GrackleFieldData.HII_density                    = malloc(sizeof(gr_float));
    All.GrackleFieldData.HeI_density                    = malloc(sizeof(gr_float));
    All.GrackleFieldData.HeII_density                   = malloc(sizeof(gr_float));
    All.GrackleFieldData.HeIII_density                  = malloc(sizeof(gr_float));
    All.GrackleFieldData.H2I_density                    = malloc(sizeof(gr_float));
    All.GrackleFieldData.H2II_density                   = malloc(sizeof(gr_float));
    All.GrackleFieldData.HM_density                     = malloc(sizeof(gr_float));
    All.GrackleFieldData.DI_density                     = malloc(sizeof(gr_float));
    All.GrackleFieldData.DII_density                    = malloc(sizeof(gr_float));
    All.GrackleFieldData.HDI_density                    = malloc(sizeof(gr_float));

    /* electron density */
    *All.GrackleFieldData.e_density                     = *All.GrackleFieldData.density * ne_guess;
    /* H and He species */
    *All.GrackleFieldData.HI_density                    = *All.GrackleFieldData.density * SphP[target].grHI;       /* initialized with HYDROGEN_MASSFRAC */
    *All.GrackleFieldData.HII_density                   = *All.GrackleFieldData.density * SphP[target].grHII;
    *All.GrackleFieldData.HeI_density                   = *All.GrackleFieldData.density * SphP[target].grHeI;
    *All.GrackleFieldData.HeII_density                  = *All.GrackleFieldData.density * SphP[target].grHeII;
    *All.GrackleFieldData.HeIII_density                 = *All.GrackleFieldData.density * SphP[target].grHeIII;
    /* all the other species are set to tiny, I don't think Grackle uses them anyway */
    *All.GrackleFieldData.H2I_density                   = *All.GrackleFieldData.density * tiny;
    *All.GrackleFieldData.H2II_density                  = *All.GrackleFieldData.density * tiny;
    *All.GrackleFieldData.HM_density                    = *All.GrackleFieldData.density * tiny;
    *All.GrackleFieldData.DI_density                    = *All.GrackleFieldData.density * tiny;
    *All.GrackleFieldData.DII_density                   = *All.GrackleFieldData.density * tiny;
    *All.GrackleFieldData.HDI_density                   = *All.GrackleFieldData.density * tiny;
    All.GrackleFieldData.RT_heating_rate                = malloc(sizeof(gr_float));
    All.GrackleFieldData.RT_HI_ionization_rate          = malloc(sizeof(gr_float));
    All.GrackleFieldData.RT_HeI_ionization_rate         = malloc(sizeof(gr_float));
    All.GrackleFieldData.RT_HeII_ionization_rate        = malloc(sizeof(gr_float));
    All.GrackleFieldData.RT_H2_dissociation_rate        = malloc(sizeof(gr_float));
    All.GrackleFieldData.RT_HM_dissociation_rate        = malloc(sizeof(gr_float));
    // Transform rates to code units
    *All.GrackleFieldData.RT_heating_rate               = 0.0;
    *All.GrackleFieldData.RT_HI_ionization_rate         = 0.0;
    *All.GrackleFieldData.RT_HeI_ionization_rate        = 0.0;
    *All.GrackleFieldData.RT_HeII_ionization_rate       = 0.0;
    *All.GrackleFieldData.RT_H2_dissociation_rate       = 0.0; 
    *All.GrackleFieldData.RT_HM_dissociation_rate       = 0.0; 

// #if (GRACKLE_CHEMISTRY >= 2)      /* Atomic+(H2I+H2II+HM) */
//     *All.GrackleFieldData.H2I_density                   = *All.GrackleFieldData.density * SphP[target].grH2I;
//     *All.GrackleFieldData.H2II_density                  = *All.GrackleFieldData.density * SphP[target].grH2II;
//     *All.GrackleFieldData.HM_density                    = *All.GrackleFieldData.density * SphP[target].grHM;
// #if (GRACKLE_CHEMISTRY >= 3)      /* Atomic+(H2I+H2II+HM)+(DI+DII+HD) */
//     *All.GrackleFieldData.DI_density                    = *All.GrackleFieldData.density * SphP[target].grDI;
//     *All.GrackleFieldData.DII_density                   = *All.GrackleFieldData.density * SphP[target].grDII;
//     *All.GrackleFieldData.HDI_density                   = *All.GrackleFieldData.density * SphP[target].grHDI;
// #endif
// #endif
#endif /* non-eq chemistry */

    /* call to the Grackle functions
     * remember: Grackle3 does not distinguish between non-equilibrium and tabulated version
     */
    switch(mode) {
        case 0:  /* returns the new internal energy; in non-eq run it evolves the abundances */
        {
    	    if(solve_chemistry(&All.GrackleUnits,
    				&All.GrackleFieldData,
    				dt) == 0)	{
    		    terminate("GRACKLE: Error in solve_chemistry.\n");
            }

/* if non-eq chemistry assign abundances back */ /*TODO:  */
#if (GRACKLE_CHEMISTRY >= 1)
            SphP[target].Ne                             = *All.GrackleFieldData.e_density     / *All.GrackleFieldData.density;
            SphP[target].grHI                           = *All.GrackleFieldData.HI_density    / *All.GrackleFieldData.density;
            SphP[target].grHII                          = *All.GrackleFieldData.HII_density   / *All.GrackleFieldData.density;
            SphP[target].grHeI                          = *All.GrackleFieldData.HeI_density   / *All.GrackleFieldData.density;
            SphP[target].grHeII                         = *All.GrackleFieldData.HeII_density  / *All.GrackleFieldData.density;
            SphP[target].grHeIII                        = *All.GrackleFieldData.HeIII_density / *All.GrackleFieldData.density;
// #if (COOL_GRACKLE_CHEMISTRY >= 2)
//             SphP[target].grH2I                          = *All.GrackleFieldData.H2I_density   / *All.GrackleFieldData.density;
//             SphP[target].grH2II                         = *All.GrackleFieldData.H2II_density  / *All.GrackleFieldData.density;
//             SphP[target].grHM                           = *All.GrackleFieldData.HM_density    / *All.GrackleFieldData.density;
// #if (COOL_GRACKLE_CHEMISTRY >= 3)
//             SphP[target].grDI                           = *All.GrackleFieldData.DI_density    / *All.GrackleFieldData.density;
//             SphP[target].grDII                          = *All.GrackleFieldData.DII_density   / *All.GrackleFieldData.density;
//             SphP[target].grHDI                          = *All.GrackleFieldData.HDI_density   / *All.GrackleFieldData.density;
// #endif
// #endif
#endif

          returnval                                   = *All.GrackleFieldData.internal_energy;
          break;
        }
            
        case 1:  /* returns the instantaneous cooling time */
	      {
	        gr_float cooling_time;
	        if(calculate_cooling_time(&All.GrackleUnits, &All.GrackleFieldData, &cooling_time) == 0) {
            terminate("GRACKLE: Error in calculate_cooling_time.\n");
          }
          returnval = cooling_time;
          break;
	      }

        case 2:  /* returns the gas temperature */
	      {
	        gr_float temperature;
          if(calculate_temperature(&All.GrackleUnits, &All.GrackleFieldData, &temperature) == 0) {
            terminate("GRACKLE: Error in calculate_temperature.\n");
          }
          returnval = temperature;
          break;
	      }
        
        case 3:  /* returns the gas pressure */
	      {
    	    gr_float pressure;
          if(calculate_pressure(&All.GrackleUnits, &All.GrackleFieldData, &pressure) == 0) {
            terminate("GRACKLE: Error in calculate_pressure.\n");
          }
          returnval = pressure;
          break;
	      }

        case 4:  /* returns gamma (effective adiabatic index) - useful only if molecules (H2,HD) are included in the chemical network */

	      {
	        gr_float gamma;
          if(calculate_gamma(&All.GrackleUnits, &All.GrackleFieldData, &gamma) == 0) {
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
    free(All.GrackleFieldData.RT_heating_rate);
    free(All.GrackleFieldData.RT_HI_ionization_rate);
    free(All.GrackleFieldData.RT_HeI_ionization_rate);
    free(All.GrackleFieldData.RT_HeII_ionization_rate);
    free(All.GrackleFieldData.RT_H2_dissociation_rate);
    free(All.GrackleFieldData.RT_HM_dissociation_rate);
#endif

    return returnval;
}


/* Function that initialises Grackle */
void InitGrackle(void)
{
    // Check the consistency
    if (gr_check_consistency() != GR_SUCCESS) {
      terminate("GRACKLE: Error in gr_check_consistency.\n");
    }
    
    int grackle_verbose = 1;
    // Enable output
    if(ThisTask == 0) grackle_verbose = 1;
    
    // First, set up the units system.
    // These are conversions from code units to cgs.
    /* 1 for cosmo, 0 if not; you can put it equal to All.ComovingIntegrationOn in order to set it automatically.
     * Beware though that for a non-cosmo sim a_value has to be kept constant.
     * IMPORTANT: at the moment the density is already converted into proper frame before calling the cooling routine,
     * so it's better to set this = 0 and only change a_value at every timestep.
     */
    All.GrackleUnits.comoving_coordinates = 0; //All.ComovingIntegrationOn; // 1 if cosmological sim, 0 if not
    All.GrackleUnits.density_units        = All.UnitDensity_in_cgs;
    All.GrackleUnits.length_units         = All.UnitLength_in_cm;
    All.GrackleUnits.time_units           = All.UnitTime_in_s;
    All.GrackleUnits.a_units              = 1.0; // units for the expansion factor; NOTE: Should be 1 always
    set_velocity_units(&All.GrackleUnits);

    if(ThisTask == 0)
        if (fabs(All.GrackleUnits.velocity_units - All.UnitVelocity_in_cm_per_s) > 1e-9) { 
        terminate("GRACKLE: Velocity units doesn't match; Check the units!");
    }
    
    // Set initial expansion factor (for internal units).
    // Set expansion factor to 1 for non-cosmological simulation.
    /* Set initial expansion factor (for internal units).
     * Set expansion factor to 1 (or to the scale factor that corresponds to the right redshift) for non-cosmological simulation.
     * This is the only parameter that we need to change at every timestep in the GrackleUnits struct.
     * All the others are constant throughout the simulation, since the density is already converted into proper frame before calling the cooling routine.
     */
    double a_value = 1.0; 
    if(All.ComovingIntegrationOn) 
      {a_value = All.TimeBegin;}
    All.GrackleUnits.a_value = a_value;
    
    // Second, create a chemistry object for parameters and rate data.
    chemistry_data *my_grackle_data;
    my_grackle_data = malloc(sizeof(chemistry_data));
    
    if (set_default_chemistry_parameters(my_grackle_data) == 0) {
      terminate("\nGRACKLE: Error in set_default_chemistry_parameters.\n");
    }
    
    /* Third, fill this chemistry object to set parameter values for chemistry & cooling. */    
    
    /* Flag to control which three-body H2 formation rate is used.
     *    0: Abel, Bryan & Norman (2002),
     *    1: Palla, Salpeter & Stahler (1983),
     *    2: Cohen & Westberg (1983),
     *    3: Flower & Harris (2007),
     *    4: Glover (2008).
     *    These are discussed in Turk et. al. (2011). Default: 0.
     */
    my_grackle_data->three_body_rate        = 0;
   
#ifdef METALS // TODO: 
    /* Flag to enable metal cooling using the Cloudy tables. If enabled, the cooling table to be used must be specified with the grackle_data_file parameter. Default: 0. */
    my_grackle_data->metal_cooling                     = 1;
    /* Flag to enable H2 formation on dust grains, dust cooling, and dust-gas heat transfer follow Omukai (2000). This assumes that the dust to gas ratio scales with the metallicity. Default: 0. */
    my_grackle_data->h2_on_dust                        = 0;
    /* Flag to enable a spatially uniform heating term approximating photo-electric heating from dust from Tasker & Bryan (2008). Default: 0.
     * If photoelectric_heating enabled, photoelectric_heating_rate is the heating rate in units of erg cm-3 s-1. Default: 8.5e-26.
     * This is not adjusted to local background.
     * (Caution: this tends to heat gas even at extremely high densities to ~3000 K, when it should be entirely self-shielding)
     * Another comment: this heats also all the low-density gas to 10^4 K, so it's better not to use it.
     */
    my_grackle_data->photoelectric_heating             = 0;  /* read above but DO NOT USE */
    my_grackle_data->photoelectric_heating_rate        = 8.5e-26;
#else
    my_grackle_data->metal_cooling                     = 1;
    my_grackle_data->h2_on_dust                        = 0;
    my_grackle_data->photoelectric_heating             = 0;
    my_grackle_data->photoelectric_heating_rate        = 8.5e-26;
#endif

    /* Flags that are useful only if you provide an array of heating rates, so DO NOT USE
     * my_grackle_data->use_specific_heating_rate      = 0;
     * my_grackle_data->use_volumetric_heating_rate    = 0;
     */

    /* Flag to enable an effective CMB temperature floor.
     * This is implemented by subtracting the value of the cooling rate at TCMB from the total METAL cooling rate. Default: 1.
     * Beware! You could still have Tgas<TCBM because it imposes a temperature floor only for the metal cooling.
     */
    my_grackle_data->cmb_temperature_floor             = 1;

    /* Flag to enable a UV background.
     * If enabled, the cooling table to be used must be specified with the grackle_data_file parameter. Default: 0.
     */
    my_grackle_data->UVbackground                      = 1;
    /* The following flags are related to the UVB, but they are automatically set to the right values, so do not need to use. These numbers are the correct ones for FG2011 UVB.
     * my_grackle_data->UVbackground_redshift_on       = 10.6;
     * my_grackle_data->UVbackground_redshift_off      = 0;
     * my_grackle_data->UVbackground_redshift_fullon   = 10.6;
     * my_grackle_data->UVbackground_redshift_drop     = 0;
     */

    /* Flag to enable Compton heating from an X-ray background following Madau & Efstathiou (1999). Default: 0.
     * Beware: this flag is just broken, it causes a runaway ionisation and heating as soon as the UVB kicks in
     * and you start having some electrons. DO NOT use it.
     */
    my_grackle_data->Compton_xray_heating              = 0;  /* read above but DO NOT USE */
    
    
    /* Flag to enable H2 collision-induced emission cooling from Ripamonti & Abel (2004). Default: 0. */
    my_grackle_data->cie_cooling                       = 0;
    /* Flag to enable H2 cooling attenuation from Ripamonti & Abel (2004). Default: 0. */
    my_grackle_data->h2_optical_depth_approximation    = 0;
    
    /* Rad Intensity of a constant Lyman-Werner H2 photo-dissociating radiation field,
     * in units of 10-21 erg s-1 cm-2 Hz-1 sr-1. Default: 0.
     */
    my_grackle_data->LWbackground_intensity            = 0;
    /* Flag to enable suppression of Lyman-Werner flux due to Lyman-series absorption
     * (giving a sawtooth pattern), taken from Haiman & Abel, & Rees (2000). Default: 0.
     */
    my_grackle_data->LWbackground_sawtooth_suppression = 0;


    /* flag to include RT */
    /* This enables H2I, HM, H2II rates corrections 
     * in grackle per particle but also allows
     * the possibility of using a background radiation
     * */
    my_grackle_data->use_radiative_transfer            = 0;

    /* For both UV bkgd and RT; options for length scale:
     *     1: Sobolev-like (from WG11)
     *     2: array of lengths
     *     3: local Jeans length
     * Default: 0.
     */
    my_grackle_data->H2_self_shielding                 = 3;

    /* Flag for self-shielding from UV bkgd. Default: 0. */
    my_grackle_data->self_shielding_method             = 0;


    /* main flags */
    
    /* Flag to activate the grackle machinery */
    my_grackle_data->use_grackle                       = 1;
    /* Path to the data file containing the metal cooling and UV background tables (for non-eq mode) and metal and primordial cooling/heating (for equilibrium mode) */
    my_grackle_data->grackle_data_file                  = All.GrackleDataFile;
    /* Flag to include radiative cooling and actually update the thermal energy during the
     * chemistry solver. If off, the chemistry species will still be updated. The most
     * common reason to set this to off is to iterate the chemistry network to an equilibrium state. Default: 1.
     */
    my_grackle_data->with_radiative_cooling            = 1;

    /* The ratio of specific heats for an ideal gas. A direct calculation for the molecular component is used if primordial_chemistry > 1. Default: 5/3. */
    my_grackle_data->Gamma                             = GAMMA;              /* our eos set in Config.sh */
    /* Flag to control which primordial chemistry network is used (set by Config.sh) */
#ifndef GRACKLE_CHEMISTRY
    my_grackle_data->primordial_chemistry              = 0;                          /* if nothing is set assume fully tabulated cooling */
#else
    my_grackle_data->primordial_chemistry              = GRACKLE_CHEMISTRY;
#endif

    
    /* Finally, initialize the chemistry object. This has to be the last step of the initialisation. */
    if (initialize_chemistry_data(&All.GrackleUnits) == 0) {
        terminate("Error in initialize_chemistry_data.\n");
    }
    
    if(ThisTask == 0)
        printf("GRACKLE: Grackle Initialized\n");
}

#endif  /* USE_GRACKLE */
