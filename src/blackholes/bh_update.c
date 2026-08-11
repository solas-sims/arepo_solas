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

/* Sph loop kernel function */
void bh_kernel(double u, double hinv3, double hinv4, double *wk, double *dwk)
{
  // Cubic spline
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
}

integertime bh_timestep(int i)
{ 
  /* Neighbours minimum bin */
  double dt_ngbmin = (BhP[i].NgbsMinBin ? (((integertime)1) << BhP[i].NgbsMinBin) : 0) * All.Timebase_interval;
  
  double dt;
  
  if(dt_ngbmin != 0)
    dt = dt_ngbmin;
  else 
    dt = TIMEBASE * All.Timebase_interval;

#ifdef BH_ACCRETION_ACTIVE
  /* Accretion timescale */
  double bh_timestep = (BhP[i].TimeBinBh ? (((integertime)1) << BhP[i].TimeBinBh) : 0) * All.Timebase_interval;
  
  double dt_bh;
  
  if(BhP[i].Accretion != 0.0 && bh_timestep != 0.0)
    dt_bh = PPB(i).Mass / (BhP[i].Accretion / bh_timestep);
  else
    dt_bh = TIMEBASE * All.Timebase_interval;

  if(dt_bh < dt)
    dt = dt_bh;
#endif

  integertime ti_step = (integertime)(dt / All.Timebase_interval);
  
  return ti_step;
}

void bh_update_timesteps(void)
{
  int idx, i;

  for(idx = 0; idx < TimeBinsBh.NActiveParticles; idx++)
    {
      i = TimeBinsBh.ActiveParticleList[idx];
    
#if defined(SELFGRAVITY) ||  defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE)
      BhP[i].TimeBinBh = PPB(i).TimeBinGrav;
#else
      int bin;
      timebins_get_bin_and_do_validity_checks(bh_timestep(i), &bin, BhP[i].TimeBinBh);
      BhP[i].TimeBinBh = bin;
#endif
    }
    
  bh_reconstruct_timebins();
  bh_update_list_of_active_particles();
}

/* Call this function as the reconstruct_timebins() bh version */
void bh_reconstruct_timebins(void)
{
  int i, bin;

  for(bin = 0; bin < TIMEBINS; bin++)
    {
      TimeBinsBh.TimeBinCount[bin]   = 0;
      TimeBinsBh.FirstInTimeBin[bin] = -1;
      TimeBinsBh.LastInTimeBin[bin]  = -1;
    }
  
  for(i = 0; i < NumBhs; i++)
    {
      bin = BhP[i].TimeBinBh;
      if(bin >= TIMEBINS)
        continue;

      if(TimeBinsBh.TimeBinCount[bin] > 0)
        {
          TimeBinsBh.PrevInTimeBin[i]                                  = TimeBinsBh.LastInTimeBin[bin];
          TimeBinsBh.NextInTimeBin[i]                                  = -1;
          TimeBinsBh.NextInTimeBin[TimeBinsBh.LastInTimeBin[bin]]      = i;
          TimeBinsBh.LastInTimeBin[bin]                                = i;
        }
      else
        {
          TimeBinsBh.FirstInTimeBin[bin] = TimeBinsBh.LastInTimeBin[bin] = i;
          TimeBinsBh.PrevInTimeBin[i] = TimeBinsBh.NextInTimeBin[i] = -1;
        }
      TimeBinsBh.TimeBinCount[bin]++;
    }
}

/* Call this function after updating the bh-timebin to the ngb condition */
void bh_update_list_of_active_particles(void)
{
  int i, n;
  TimeBinsBh.NActiveParticles = 0;

  for(n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinSynchronized[n]) 
        {
          for(i = TimeBinsBh.FirstInTimeBin[n]; i >= 0; i = TimeBinsBh.NextInTimeBin[i])
            {
              TimeBinsBh.ActiveParticleList[TimeBinsBh.NActiveParticles] = i;
              TimeBinsBh.NActiveParticles++;  
            }
        }
    }

  mysort(TimeBinsBh.ActiveParticleList, TimeBinsBh.NActiveParticles, sizeof(int), int_compare);

  sumup_large_ints(1, &TimeBinsBh.NActiveParticles, &TimeBinsBh.GlobalNActiveParticles);
}

void bh_perform_end_of_step_physics(void)
{
  int idx, i;
    
#ifdef BH_ACCRETION_ACTIVE
  /* Accrete mass and angular momentum onto the bh */
  for(idx = 0; idx < TimeBinsBh.NActiveParticles; idx++)
    {
      i = TimeBinsBh.ActiveParticleList[idx];

      PPB(i).Mass += (1 - All.Epsilon_r) * BhP[i].Accretion;

      /* Do we want to track bh spin? */
      //BhP[i].AngularMomentum[0] += BhP[i].Accretion * BhP[i].VelocityGasCircular[0];
      //BhP[i].AngularMomentum[1] += BhP[i].Accretion * BhP[i].VelocityGasCircular[1];
      //BhP[i].AngularMomentum[2] += BhP[i].Accretion * BhP[i].VelocityGasCircular[2];
    }
      
  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Type != 0 || P[i].Mass == 0 || P[i].ID == 0)
        continue;

      double original_mass = P[i].Mass;

      P[i].Mass -= SphP[i].BhMassDrain;

      double factor = P[i].Mass / original_mass;

      SphP[i].Energy *= factor;
      SphP[i].Momentum[0] *= factor;
      SphP[i].Momentum[1] *= factor;
      SphP[i].Momentum[2] *= factor;

#ifdef MAXSCALARS
      for(int s = 0; s < N_Scalar; s++)
        *(MyFloat *)(((char *)(&SphP[i])) + scalar_elements[s].offset_mass) *= factor;
#endif

      SphP[i].BhMassDrain = 0;
    }
#endif

#ifdef BH_FEEDBACK_ACTIVE   
  struct pv_update_data pvd;
  if(All.ComovingIntegrationOn)
    {
      pvd.atime    = All.Time;
      pvd.hubble_a = hubble_function(All.Time);
      pvd.a3inv    = 1 / (All.Time * All.Time * All.Time);
    }
  else
      pvd.atime = pvd.hubble_a = pvd.a3inv = 1.0;
    
  /* Inject feedback to ngb cells */ 

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
            
      if(i < 0)
        continue;
                
      // Dump mass, momentum and energy injected by bh 
      if(SphP[i].BhThermalFeed > 0)
        {
          // Add mass
          //P[i].Mass += SphP[i].BhMassFeed;
                    
          // Add kinetic energy */
          //SphP[i].Energy += SphP[i].KineticFeed;
                    
          // Calculate momentum feed exactly so energy is conserved */
          //-> we need to do this here so that particle properties don't change between loading the buffer and emptying it*/
          //kick_vector[0] = SphP[i].BhKickVector[0];
          //kick_vector[1] = SphP[i].BhKickVector[1];
          //kick_vector[2] = SphP[i].BhKickVector[2];
                    
          //p0 = sqrt(pow(SphP[i].Momentum[0], 2) + pow(SphP[i].Momentum[1], 2) + pow(SphP[i].Momentum[2], 2));
                    
          //pj = sqrt(2 * P[i].Mass * (SphP[i].Energy - (P[i].Mass-SphP[i].MassLoading)*SphP[i].Utherm*All.cf_atime*All.cf_atime)) - p0;

          // Update total energy 
          SphP[i].Energy += SphP[i].BhThermalFeed;
          All.BhFeedbackLocal[1] += SphP[i].BhThermalFeed;
          // Update momentum 
          //SphP[i].Momentum[0] += kick_vector[0] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2));
          //SphP[i].Momentum[1] += kick_vector[1] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2));
          //SphP[i].Momentum[2] += kick_vector[2] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2));
          // Update velocities
          //update_primitive_variables_single(P, SphP, i, &pvd);
          // Update internal energy 
          update_internal_energy(P, SphP, i, &pvd);
          // Update pressure
          set_pressure_of_cell_internal(P, SphP, i);
          // Set feed flags to zero 
          SphP[i].BhMassFeed = SphP[i].BhThermalFeed = SphP[i].BhKineticFeed = 0;

#ifdef JET_TRACER
          // Tracer field advected passively 
          SphP[i].PScalars[JET_INDEX] = 1;
          sync_conserved_from_primitive(i, JET_INDEX);
#endif
        }
    }
#ifdef BURST_MODE
      All.FeedbackFlag = -1;
#endif
        
    MPI_Allreduce(&All.BhFeedbackLocal, &All.BhFeedbackGlobal, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
    mpi_printf("BLACKHOLES: Number of active blackholes = %lld \n", TimeBinsBh.GlobalNActiveParticles);
    mpi_printf("BLACKHOLES: Energy given by BHs = %e, Energy taken up by gas particles = %e \n",
    All.BhFeedbackGlobal[0], All.BhFeedbackGlobal[1]);

#ifdef BURST_MODE
    if(All.EnergyExchangeTot[0] - All.EnergyExchangeTot[1] > 10)
        All.FeedbackFlag = 1;
#endif

#endif
} // perform_end_of_step_physics(void)
