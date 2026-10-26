#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "../main/allvars.h"
#include "../main/proto.h"


static int int_compare(const void *a, const void *b)
{
  if(*((int *)a) < *((int *)b))
    return -1;

  if(*((int *)a) > *((int *)b))
    return +1;

  return 0;
}

double gaussian_weight(double r, double h)
{
  double sigma = h / 2.0; 
  double x = r / sigma;
  return exp(-0.5 * x * x);
}

/* Sph loop kernel function */
/*void star_kernel(double u, double hinv3, double hinv4, double *wk, double *dwk)
{
  double K_norm = 8.0 / M_PI;

  if(u < 0.5)
    {
      *dwk = u * (18.0 * u - 12.0);
      
      *wk = (1.0 + 6.0 * (u - 1.0) * u * u);
    }
  else
    {
      double t1 = (1.0 - u);
      double t2 = t1 * t1;
      
      *dwk = -6.0 * t2;
      
      *wk = 2.0 * t2 * t1;
    }

  *dwk *= K_norm * hinv4;
  
  *wk  *= K_norm * hinv3;
}*/

void star_in(void)
{
  feedback_allocate(&MechanicalFeedbackEvents, NumStars);
}

void star_exit(void)
{
  feedback_free(&MechanicalFeedbackEvents);
}

void feedback_init(struct Mechanical_Feedback_Events *MFEvents)
{
  MFEvents->NumEvents = 0;
  MFEvents->MaxEvents = 0;

  MFEvents->MechanicalFeedbackData = NULL;
}

void feedback_allocate(struct Mechanical_Feedback_Events *MFEvents, int MaxEvents)
{
  MFEvents->NumEvents = 0;
  MFEvents->MaxEvents = MaxEvents > 0 ? MaxEvents : ALLOC_STAR_ROOM;
  MFEvents->MechanicalFeedbackData = malloc(MFEvents->MaxEvents * sizeof(Mechanical_Feedback_Data));
}

void feedback_reallocate(struct Mechanical_Feedback_Events *MFEvents, int NewMaxEvents)
{
  MFEvents->MaxEvents = NewMaxEvents;

  MFEvents->MechanicalFeedbackData = realloc(MFEvents->MechanicalFeedbackData, NewMaxEvents * sizeof(Mechanical_Feedback_Data));
}

void feedback_free(struct Mechanical_Feedback_Events *MFEvents)
{
  free(MFEvents->MechanicalFeedbackData);

  feedback_init(MFEvents);
}

integertime star_timestep(int i)
{ 
  /* Host hydro bin */
  double dt_host = (SP[i].HostHydroBin ? (((integertime)1) << SP[i].HostHydroBin) : 0) * All.Timebase_interval;
  
  double dt;
  
  if(dt_host != 0)
    dt = dt_host;
  else 
    dt = TIMEBASE * All.Timebase_interval;

  /* Set a maximum star timestep at 0.01 Myr */
  double dt_star = 1.0e4 / All.cf_UnitTime_in_yr;

  if(dt_star < dt)
    dt = dt_star;

  /* Park dead or low mass stars */
  if(SP[i].Active < 0)
    dt = TIMEBASE * All.Timebase_interval;
    
  integertime ti_step = (integertime)(dt / All.Timebase_interval);
  
  return ti_step;
}

void star_update_timesteps(void)
{
  int idx, i;

  for(idx = 0; idx < TimeBinsStar.NActiveParticles; idx++)
    {
      i = TimeBinsStar.ActiveParticleList[idx];

      /* Park dead or low mass stars */
      if(SP[i].Active < 0)
        {
          SP[i].TimeBinStar = TIMEBINS;
          continue;
        }

#if defined(SELFGRAVITY) ||  defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE)
      SP[i].TimeBinStar = PPS(i).TimeBinGrav;
#else
      int bin;
      timebins_get_bin_and_do_validity_checks(star_timestep(i), &bin, SP[i].TimeBinStar);
      SP[i].TimeBinStar = bin;
#endif
    }
    
  star_reconstruct_timebins();
  star_update_list_of_active_particles();
}

/* Call this function as the reconstruct_timebins() star version */
void star_reconstruct_timebins(void)
{
  int i, bin;

  for(bin = 0; bin < TIMEBINS; bin++)
    {
      TimeBinsStar.TimeBinCount[bin] = 0;
      TimeBinsStar.FirstInTimeBin[bin] = -1;
      TimeBinsStar.LastInTimeBin[bin] = -1;
    }
  
  for(i = 0; i < NumStars; i++)
    {
      
      bin = SP[i].TimeBinStar;
      if(bin >= TIMEBINS)
        continue;

      if(TimeBinsStar.TimeBinCount[bin] > 0)
        {
          TimeBinsStar.PrevInTimeBin[i] = TimeBinsStar.LastInTimeBin[bin];
          TimeBinsStar.NextInTimeBin[i] = -1;
          TimeBinsStar.NextInTimeBin[TimeBinsStar.LastInTimeBin[bin]] = i;
          TimeBinsStar.LastInTimeBin[bin] = i;
        }
      else
        {
          TimeBinsStar.FirstInTimeBin[bin] = TimeBinsStar.LastInTimeBin[bin] = i;
          TimeBinsStar.PrevInTimeBin[i] = TimeBinsStar.NextInTimeBin[i] = -1;
        }
      TimeBinsStar.TimeBinCount[bin]++;
    }
}

/* Call this function after updating the star-timebin */
void star_update_list_of_active_particles(void)
{
  int i, n;
  TimeBinsStar.NActiveParticles = 0;

  for(n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinSynchronized[n]) 
        {
          for(i = TimeBinsStar.FirstInTimeBin[n]; i >= 0; i = TimeBinsStar.NextInTimeBin[i])
            {
              TimeBinsStar.ActiveParticleList[TimeBinsStar.NActiveParticles] = i;
              TimeBinsStar.NActiveParticles++;  
            }
        }
    }

  mysort(TimeBinsStar.ActiveParticleList, TimeBinsStar.NActiveParticles, sizeof(int), int_compare);

  sumup_large_ints(1, &TimeBinsStar.NActiveParticles, &TimeBinsStar.GlobalNActiveParticles);
}

static inline void deactivate_star(int i)
{
#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  SP[i].TimeToSN = MAX_REAL_NUMBER;
  SP[i].NextSNEnergy = 0.0;
#endif
  
  SP[i].Active = STAR_INACTIVE; 
  SP[i].WithFeedback = 0; 
  
  SP[i].Hsml = 0.0;
  SP[i].DensityFlag = -1;
  SP[i].NgbsMass = 0.0;
  SP[i].NgbsVolume = 0.0;
  SP[i].HostHydroBin = TIMEBINS;
  SP[i].TimeBinStar = TIMEBINS;
}

/* Compute feedback properties of active stars */
void star_prep(void)
{
  TIMER_START(CPU_STARS_PREP);

  int idx, i;
  
  for(idx = 0; idx < TimeBinsStar.NActiveParticles; idx++)
    {
      i = TimeBinsStar.ActiveParticleList[idx];

      /* Put newly formed stars on the main sequence */
      if(SP[i].Active == STAR_UNBORN)
        {
          SP[i].MassOfStar = PPS(i).Mass;
          SP[i].Active = STAR_ACTIVE;
          
          SP[i].Age = 0.0;
          SP[i].Birthtime = All.Time;
        }
      
      /* Advance timestep and age */
      MyDouble star_timestep = (SP[i].TimeBinStar ? (((integertime)1) << SP[i].TimeBinStar) : 0) * All.Timebase_interval;

      /* Convert properties to yr and msun. Birthtime/All.Time are scale factors in comoving runs,
       * so the age must come from the cosmic-time difference rather than All.Time - Birthtime;
       * get_time_difference_in_Gyr() handles both cases (and includes the 1/h conversion) */
      MyDouble star_mass_msun = SP[i].MassOfStar * All.cf_UnitMass_in_Msun;
      MyDouble star_timestep_yr = star_timestep * All.cf_UnitTime_in_yr;
      MyDouble star_age_yr = get_time_difference_in_Gyr(SP[i].Birthtime, All.Time) * 1.0e9;

      /* Age is kept in the same (dln a, if comoving) units as star_timestep */
      SP[i].Age = star_age_yr / All.cf_UnitTime_in_yr;

#ifdef METALS 
      MyDouble star_metallicity = SP[i].Metallicity;
#else 
      MyDouble star_metallicity = 0;
#endif

      Star_Feedback StarFeedback;

      /* Call the interpolation functions */
#if defined(STAR_PARTICLES) && STAR_PARTICLES < 2
      StarFeedback = units_for_feedback(star_particle_feedback(i, star_timestep_yr, star_metallicity, star_age_yr));
#elif STAR_PARTICLES == 2     
      StarFeedback = units_for_feedback(star_feedback_compute(star_timestep_yr, star_metallicity, star_mass_msun, star_age_yr));
#endif

      /* Deactivate dead or low mass stars */
      if(StarFeedback.State == -1)
        {
          deactivate_star(i);
          continue;
        }

      /* Assign stellar feedback variables */  
#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
      SP[i].TimeToSN = StarFeedback.TimeToSN;
      SP[i].NextSNEnergy = StarFeedback.NextSNEnergy;
#endif

      for(int k = 0; k < 3; k++)
        {
          SP[i].MechanicalFeedback.StarPosition[k] = PPS(i).Pos[k];
          SP[i].MechanicalFeedback.StarVelocity[k] = PPS(i).Vel[k];
        } 

#ifdef WINDS
      SP[i].MechanicalFeedback.MassLoss = StarFeedback.MassLoss;
#if GRACKLE_CHEMISTRY >= 1
      SP[i].MechanicalFeedback.HLoss = StarFeedback.HLoss;
      SP[i].MechanicalFeedback.HeLoss = StarFeedback.HeLoss;
#endif
#ifdef METALS
      SP[i].MechanicalFeedback.MetalsLoss = StarFeedback.MetalsLoss;
#endif
      SP[i].MechanicalFeedback.WindMomentum = StarFeedback.WindMomentum;
#endif

#ifdef STAR_RADIATION_ACTIVE      
      for(int w = 0; w < WAVEBANDS; w++)
        {
          SP[i].MechanicalFeedback.Radiated[w].Photons = StarFeedback.Radiated[w].Photons;
          SP[i].MechanicalFeedback.Radiated[w].Energy = StarFeedback.Radiated[w].Energy;
        }
#endif

#ifdef SUPERNOVAE
      SP[i].MechanicalFeedback.SN_MassLoss = StarFeedback.SN_MassLoss;
#if GRACKLE_CHEMISTRY >= 1
      SP[i].MechanicalFeedback.SN_HLoss = StarFeedback.SN_HLoss;
      SP[i].MechanicalFeedback.SN_HeLoss = StarFeedback.SN_HeLoss;
#endif
#ifdef METALS
      SP[i].MechanicalFeedback.SN_MetalsLoss = StarFeedback.SN_MetalsLoss;
#endif
      SP[i].MechanicalFeedback.SN_EnergyInject = StarFeedback.SN_EnergyInject;
#endif

      /* Determine if star provides feedback */
      /* If with_feedback == 0 the star is skipped in the feedback functions */
      int with_feedback = 0;

#ifdef WINDS 
      if(SP[i].MechanicalFeedback.MassLoss > 0)
        with_feedback++;
#endif

#ifdef STAR_RADIATION_ACTIVE
      for(int w = 0; w < WAVEBANDS; w++)
        {
          if(SP[i].MechanicalFeedback.Radiated[w].Photons > 0 || SP[i].MechanicalFeedback.Radiated[w].Energy > 0)
            { 
              with_feedback++; 
            }
        }  
#endif

#ifdef SUPERNOVAE
      if(SP[i].MechanicalFeedback.SN_MassLoss > 0 || SP[i].MechanicalFeedback.SN_EnergyInject > 0)
        with_feedback++;
#endif

      SP[i].WithFeedback = with_feedback;
    }

  TIMER_STOP(CPU_STARS_PREP);
}

/*! Fraction of a host cell's mass that must survive the net removal by SN host sweeping (matches the 0.1 floor used when the
 *  removal is sized in star_feedback.c). */
#define STAR_FEEDBACK_MIN_HOST_MASS_FRACTION 0.1

/*! Cumulative host-removal limiter statistics (this task): [0] cells limited, [1] mass, [2] energy not removed. */
static double FeedbackLimiter[3] = {0, 0, 0};

/* The host-cell removal of an SN event (star_feedback.c) is sized from the cell's mass at the time of the event, but without
 * STAR_FEEDBACK_IMMEDIATE_DEPOSIT it only reaches the cell when the cell is next active. A cell on a long time step exchanges
 * fluxes with its short-step neighbours in between and can have lost mass, so the pending net removal can exceed what is left
 * (and momentum and energy go negative with it). Limit the net removal so that at least STAR_FEEDBACK_MIN_HOST_MASS_FRACTION of
 * the cell survives, scaling all pending components together so the surviving gas keeps its specific momentum and energy.
 * This is not conservative -- the neighbours have already received the full amounts -- so the mass and energy that are not
 * removed are accumulated and reported in the STARS log lines. */
static void limit_host_removal(int i)
{
  if(SphP[i].StarMassFeed < 0 && P[i].Mass + SphP[i].StarMassFeed < STAR_FEEDBACK_MIN_HOST_MASS_FRACTION * P[i].Mass)
    {
      double s = (1.0 - STAR_FEEDBACK_MIN_HOST_MASS_FRACTION) * P[i].Mass / (-SphP[i].StarMassFeed);

      FeedbackLimiter[0] += 1;
      FeedbackLimiter[1] += (1.0 - s) * (-SphP[i].StarMassFeed);
      FeedbackLimiter[2] += (1.0 - s) * (-SphP[i].StarEnergyFeed);

      SphP[i].StarMassFeed *= s;
      SphP[i].StarEnergyFeed *= s;
      for(int k = 0; k < 3; k++)
        SphP[i].StarMomentumFeed[k] *= s;
#ifdef METALS
      SphP[i].StarMetalsFeed *= s;
#endif
#if GRACKLE_CHEMISTRY >= 1
      for(int sp = 0; sp < GRACKLE_SPECIES_NUMBER; sp++)
        SphP[i].StarChemFeed[sp] *= s;
#endif
    }
}

#ifdef STAR_FEEDBACK_IMMEDIATE_DEPOSIT
/*! \brief Apply the pending wind/SN feed of the given cells to their conserved variables right away.
 *
 *  Called at the end of star_feedback(), for every cell that received a kick. Unlike the end-of-step application in
 *  star_perform_end_of_step_physics() (which only reaches active cells), this updates the conserved variables (Mass, Momentum,
 *  Energy, metals, chemical species) of the cell whether or not it is active, exactly as the hydro flux exchange does. Primitive
 *  variables are left alone: active cells get them from update_primitive_variables(), inactive ones when they are next active.
 *  Removals are therefore evaluated against the same cell content they were sized from and cannot overshoot.
 *
 *  \param[in] n    number of entries in list
 *  \param[in] list local gas cell indices (duplicates are harmless)
 */
void star_flush_pending_feed(int n, const int *list)
{
  for(int k = 0; k < n; k++)
    {
      int i = list[k];
      if(i < 0 || i >= NumGas)
        continue;

      int pending = (SphP[i].StarMassFeed != 0 || SphP[i].StarEnergyFeed != 0 || SphP[i].StarMomentumFeed[0] != 0 ||
                     SphP[i].StarMomentumFeed[1] != 0 || SphP[i].StarMomentumFeed[2] != 0);
#ifdef METALS
      pending |= (SphP[i].StarMetalsFeed != 0);
#endif
      if(!pending)
        continue;

      limit_host_removal(i);

      /* Add mass */
      double old_mass_before_feed = P[i].Mass;
      P[i].Mass += SphP[i].StarMassFeed;

      if(P[i].Mass <= 0)
        {
          double floor_mass = fmax(1.0e-3 * All.TargetGasMass, 1.0e-3 * old_mass_before_feed);
          if(floor_mass <= 0)
            floor_mass = 1.0e-3 * All.TargetGasMass;

          printf("WARNING: STAR_FEEDBACK: cell ID=%lld (task=%d) driven to non-positive mass by immediate feedback deposition (old "
                 "Mass=%g, StarMassFeed=%g, resulting Mass=%g) -- clamping to %g\n",
                 (long long)P[i].ID, ThisTask, old_mass_before_feed, SphP[i].StarMassFeed, P[i].Mass, floor_mass);
          myflush(stdout);

          P[i].Mass = floor_mass;
        }

      All.StarFeedbackLocal[3] += SphP[i].StarMassFeed;
      SphP[i].StarMassFeed = 0;

#if GRACKLE_CHEMISTRY >= 1
      for(int sp = 0; sp < GRACKLE_SPECIES_NUMBER; sp++)
        {
          SphP[i].GrackleSpeciesConserved(GRACKLE_SPECIES_INDEX + sp) += SphP[i].StarChemFeed[sp];
          SphP[i].StarChemFeed[sp] = 0;
        }
#endif
#ifdef METALS
      SphP[i].GasMetals += SphP[i].StarMetalsFeed;
      All.StarFeedbackLocal[4] += SphP[i].StarMetalsFeed;
      SphP[i].StarMetalsFeed = 0;
#endif

      for(int d = 0; d < 3; d++)
        {
          SphP[i].Momentum[d] += SphP[i].StarMomentumFeed[d];
          SphP[i].StarMomentumFeed[d] = 0;
        }

      SphP[i].Energy += SphP[i].StarEnergyFeed;
      All.StarFeedbackLocal[5] += SphP[i].StarEnergyFeed;
      SphP[i].StarEnergyFeed = 0;

      /* OldMass bounds the mass a face may transport (face_limit_fluxes); keep it consistent with a lowered mass, as
       * done in the primitive-variable exchange */
      if(P[i].Mass < SphP[i].OldMass)
        SphP[i].OldMass = P[i].Mass;
    }
}
#endif /* #ifdef STAR_FEEDBACK_IMMEDIATE_DEPOSIT */

void star_perform_end_of_step_physics(void)
{
  int idx, i;

  struct pv_update_data pvd;
  if(All.ComovingIntegrationOn)
    {
      pvd.atime = All.Time;
      pvd.hubble_a = hubble_function(All.Time);
      pvd.a3inv = 1 / (All.Time * All.Time * All.Time);
    }
  else
    pvd.atime = pvd.hubble_a = pvd.a3inv = 1.0;

  /* Subtract massloss from stars */
  for(idx = 0; idx < TimeBinsStar.NActiveParticles; idx++)
    {
      i = TimeBinsStar.ActiveParticleList[idx];

#ifdef WINDS
      PPS(i).Mass -= SP[i].MechanicalFeedback.MassLoss;
#endif

#ifdef SUPERNOVAE
      PPS(i).Mass -= SP[i].MechanicalFeedback.SN_MassLoss;
#endif
    }

  /* Dump star injected mass, momentum, and energy into gas */  
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

#if defined(WINDS) || defined(SUPERNOVAE)
      limit_host_removal(i);

      /* Add mass */
      double old_mass_before_feed = P[i].Mass;
      P[i].Mass += SphP[i].StarMassFeed;

      if(P[i].Mass <= 0)
        {
          double floor_mass = fmax(1.0e-3 * All.TargetGasMass, 1.0e-3 * old_mass_before_feed);
          if(floor_mass <= 0)
            floor_mass = 1.0e-3 * All.TargetGasMass;

          printf(
              "WARNING: STAR_FEEDBACK: cell ID=%lld (task=%d, pos=%g|%g|%g) driven to non-positive mass by "
              "feedback deposition (old Mass=%g, StarMassFeed=%g, resulting Mass=%g) -- clamping to %g\n",
              (long long)P[i].ID, ThisTask, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], old_mass_before_feed,
              SphP[i].StarMassFeed, P[i].Mass, floor_mass);
          /* diagnostic: state that is left behind when only the mass is clamped (momentum and energy are not rescaled) */
          {
            double p_after[3], mom2 = 0;
            for(int k = 0; k < 3; k++)
              {
                p_after[k] = SphP[i].Momentum[k] + SphP[i].StarMomentumFeed[k];
                mom2 += p_after[k] * p_after[k];
              }
            double E_after = SphP[i].Energy + SphP[i].StarEnergyFeed;
            printf(
                "WARNING: STAR_FEEDBACK:   cell ID=%lld TimeBinHydro=%d Volume=%g  before feed: Momentum=%g|%g|%g Energy=%g  pending feed: "
                "Momentum=%g|%g|%g Energy=%g  => after clamp |v|=%g specific energy=%g (old: |v|=%g, E/m=%g)\n",
                (long long)P[i].ID, P[i].TimeBinHydro, SphP[i].Volume, SphP[i].Momentum[0], SphP[i].Momentum[1], SphP[i].Momentum[2],
                SphP[i].Energy, SphP[i].StarMomentumFeed[0], SphP[i].StarMomentumFeed[1], SphP[i].StarMomentumFeed[2],
                SphP[i].StarEnergyFeed, sqrt(mom2) / floor_mass, E_after / floor_mass,
                old_mass_before_feed > 0 ? sqrt(SphP[i].Momentum[0] * SphP[i].Momentum[0] + SphP[i].Momentum[1] * SphP[i].Momentum[1] +
                                                SphP[i].Momentum[2] * SphP[i].Momentum[2]) / old_mass_before_feed
                                         : 0.0,
                old_mass_before_feed > 0 ? SphP[i].Energy / old_mass_before_feed : 0.0);
          }
          myflush(stdout);

          P[i].Mass = floor_mass;
        }

      All.StarFeedbackLocal[3] += SphP[i].StarMassFeed;

      SphP[i].StarMassFeed = 0;
#if GRACKLE_CHEMISTRY >= 1
      for(int s = 0; s < GRACKLE_SPECIES_NUMBER; s++)
        {    
          SphP[i].GrackleSpeciesConserved(GRACKLE_SPECIES_INDEX + s) += SphP[i].StarChemFeed[s];
          
          sync_primitive_from_conserved(i, GRACKLE_SPECIES_INDEX + s);

          SphP[i].StarChemFeed[s] = 0;
        }
#endif
#ifdef METALS
      /* Add metals */
      SphP[i].GasMetals += SphP[i].StarMetalsFeed;
      All.StarFeedbackLocal[4] += SphP[i].StarMetalsFeed;
      
      sync_primitive_from_conserved(i, METALS_INDEX);

      SphP[i].StarMetalsFeed = 0;
#endif
#endif
            
#if defined(WINDS) || defined(RADIATION_PRESSURE) || defined(SUPERNOVAE) 
      /* Update momentum */ 
      SphP[i].Momentum[0] += SphP[i].StarMomentumFeed[0];
      SphP[i].Momentum[1] += SphP[i].StarMomentumFeed[1];
      SphP[i].Momentum[2] += SphP[i].StarMomentumFeed[2];
      
      /* Update velocities */ 
      update_primitive_variables_single(P, SphP, i, &pvd);
      
      /* Set feed flags to zero */
      SphP[i].StarMomentumFeed[0] = SphP[i].StarMomentumFeed[1] = SphP[i].StarMomentumFeed[2] = 0;
#endif

      /* Update total energy */
      SphP[i].Energy += SphP[i].StarEnergyFeed;
      All.StarFeedbackLocal[5] += SphP[i].StarEnergyFeed;
      
      /* Update internal energy */ 
      update_internal_energy(P, SphP, i, &pvd);
      /* Update pressure */
      set_pressure_of_cell_internal(P, SphP, i);

      /* Set feed flags to zero */
      SphP[i].StarEnergyFeed = 0;
    } // for(idx...

    MPI_Allreduce(All.StarFeedbackLocal, All.StarFeedbackGlobal, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    mpi_printf("STARS: Number of stars = %lld, active stars = %lld, feedback events = %lld \n", 
    All.TotNumStars, TimeBinsStar.GlobalNActiveParticles, MechanicalFeedbackEvents.TotEvents);
    mpi_printf("STARS: Mass given by StarParts = %e, Mass taken up by gas particles = %e \n",
    All.StarFeedbackGlobal[0], All.StarFeedbackGlobal[3]);
    mpi_printf("STARS: Metals given by StarParts = %e, Metals taken up by gas particles = %e \n",
    All.StarFeedbackGlobal[1], All.StarFeedbackGlobal[4]);
    mpi_printf("STARS: Energy given by StarParts = %e, Energy taken up by gas particles = %e \n",
    All.StarFeedbackGlobal[2], All.StarFeedbackGlobal[5]);

    double limiter_global[3];
    MPI_Allreduce(FeedbackLimiter, limiter_global, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if(limiter_global[0] > 0)
      mpi_printf("STARS: host-removal limiter (cumulative): %.0f cells limited, mass not removed = %e, energy not removed = %e \n",
                 limiter_global[0], limiter_global[1], limiter_global[2]);
} 