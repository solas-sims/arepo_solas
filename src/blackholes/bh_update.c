#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef BH_ACTIVE

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
      /* skip BHs killed (e.g. by bh_merger()) since the last rebuild; their P[]/BhP[]
       * slots persist until the next domain decomposition sweeps them up */
      if(PPB(i).Mass == 0 && PPB(i).ID == 0)
        continue;

      bin = BhP[i].TimeBinBh;
      if(bin < 0 || bin >= TIMEBINS)
        {
          /* TIMEBINS is the established "inactive/parked" sentinel elsewhere in this
           * codebase (see star_reconstruct_timebins()/deactivate_star()); a negative bin
           * has no such meaning and indicates genuinely unexpected state (ASan caught this
           * as a global-buffer-overflow read into TimeBinsBh.TimeBinCount[] via a stale
           * BhP[i].TimeBinBh). Confirmed via production runs to be state baked into the
           * restart checkpoint itself (predates any swap this run performs), not something
           * introduced live -- so *skipping* this BH here would orphan it from the timebin
           * linked list permanently: it would never again be picked up by
           * bh_update_list_of_active_particles(), silently freezing its gravity/accretion
           * out for the rest of the run. Reset to bin 0 (the same default a freshly seeded
           * BH gets, see spawn_black_hole_from_cell()) and fall through to the normal
           * insertion below instead, so it resumes being tracked -- self-healing rather than
           * permanently dropping the particle. Flagged loudly since the underlying cause
           * (how TimeBinBh went negative in the first place) is still unexplained. */
          printf("WARNING: BH_UPDATE: BhP[%d] (ID=%llu) has invalid TimeBinBh=%d (valid range [0,%d)) -- resetting to bin 0\n", i,
                 (unsigned long long)PPB(i).ID, bin, TIMEBINS);
          myflush(stdout);
          bin             = 0;
          BhP[i].TimeBinBh = 0;
        }

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

#ifdef BH_JET_FEEDBACK
/*! \brief Applies the jet buffer of active cell i (see bh_jet_feedback.c for the three injection modes).
 *
 *  Energy and momentum are the code's comoving quantities (Energy = a^2 M u + 1/2 M Vel^2, Momentum = M Vel),
 *  so the kinetic energy of a cell is 1/2 |Momentum|^2 / Mass and everything else in Energy (thermal, and
 *  magnetic if present) is left untouched by the momentum changes below.
 */
static void bh_jet_apply_cell(int i, struct pv_update_data *pvd)
{
  double dE = SphP[i].BhJetEnergy;
  double dM = SphP[i].BhJetMass;
  double dp[3] = {SphP[i].BhJetMomentum[0], SphP[i].BhJetMomentum[1], SphP[i].BhJetMomentum[2]};

  SphP[i].BhJetEnergy = SphP[i].BhJetMass = 0;
  SphP[i].BhJetMomentum[0] = SphP[i].BhJetMomentum[1] = SphP[i].BhJetMomentum[2] = 0;

  if(P[i].Mass <= 0)
    return; /* cell was dissolved or drained since the energy was buffered */

  /* Jet mass created in the cell (BH_JET_FIXED_POWER) is added at rest and with no thermal energy, as in arepobh:
   * Energy and Momentum are left alone, so only the jet energy dE changes the cell's total energy. Scaling them with
   * the mass instead would give the new mass the host cell's specific energy, which for a cell that has already been
   * kicked is an uncontrolled source of kinetic energy (measured ~2% of dE in a test with Mach ~50 cones). */
  double mass_old = P[i].Mass;

  P[i].Mass += dM;

  double mass    = P[i].Mass;
  double *p      = SphP[i].Momentum;
  double energy0 = SphP[i].Energy;

  (void)mass_old;
  (void)mass;
  (void)dE;

#if defined(MOMENTUM_JET)
  /* add momentum only; the non-kinetic part of Energy is unchanged */
  double p0sq = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];

  SphP[i].Energy -= 0.5 * p0sq / mass_old;
  for(int k = 0; k < 3; k++)
    p[k] += dp[k];
  SphP[i].Energy += 0.5 * (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) / mass;
#elif defined(THERMAL_JET)
  /* add momentum and the full jet energy; what the momentum does not take up becomes thermal */
  for(int k = 0; k < 3; k++)
    p[k] += dp[k];
  SphP[i].Energy += dE;
#else /* ENERGY_JET */
  double dnorm = sqrt(dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2]);

  SphP[i].Energy += dE;

  if(dnorm > 0)
    {
      double n[3]  = {dp[0] / dnorm, dp[1] / dnorm, dp[2] / dnorm};
      double pn    = p[0] * n[0] + p[1] * n[1] + p[2] * n[2];
      double p0sq  = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];

      /* kick s n such that 1/2 |p + s n|^2 / M = 1/2 |p|^2 / M_old + dE, i.e. all of dE becomes kinetic and the
       * non-kinetic energy of the cell is untouched:
       * s = -p.n + sqrt((p.n)^2 - |p|^2 + 2 M (1/2 |p|^2 / M_old + dE)), the positive root, exact for any angle between p and n */
      double s = -pn + sqrt(pn * pn - p0sq + 2.0 * mass * (0.5 * p0sq / mass_old + dE));

      for(int k = 0; k < 3; k++)
        p[k] += s * n[k];
    }
  /* else: no usable direction (opposing contributions cancelled); dE stays thermal */
#endif

  All.BhFeedbackLocal[1] += SphP[i].Energy - energy0;

  update_primitive_variables_single(P, SphP, i, pvd);
  update_internal_energy(P, SphP, i, pvd);
  set_pressure_of_cell_internal(P, SphP, i);

#ifdef JET_TRACER
  SphP[i].PScalars[JET_INDEX] = 1;
  sync_conserved_from_primitive(i, JET_INDEX);
#endif
}
#endif /* #ifdef BH_JET_FEEDBACK */

void bh_perform_end_of_step_physics(void)
{
  int idx, i;
    
#ifdef BH_ACCRETION_ACTIVE
  /* Accrete mass and angular momentum onto the bh */
  for(idx = 0; idx < TimeBinsBh.NActiveParticles; idx++)
    {
      i = TimeBinsBh.ActiveParticleList[idx];

      PPB(i).Mass += (1 - All.Epsilon_r) * BhP[i].Accretion;

      /* BH spin direction (BH_JET_SPIN) is tracked in bh_jet_prepare() */
    }
      
  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Type != 0 || P[i].Mass == 0 || P[i].ID == 0)
        continue;

      double original_mass = P[i].Mass;

      /* each BH caps its own contribution to BhMassDrain at 0.9x the cell's mass, but
       * multiple BHs sharing this cell as a neighbour add to the same accumulator without
       * seeing each other's contributions, so the sum can still exceed the cell's mass.
       * Guard against P[i].Mass itself already being non-positive (corrupted upstream, e.g.
       * by the mesh flux update): the 0.9x-of-current-mass clamp below only makes sense for
       * a positive mass -- applied to a negative one it flips sign and UNDER-clamps, letting
       * the drain shrink the negative value's magnitude instead of being rejected outright. */
      if(P[i].Mass > 0)
        {
          if(SphP[i].BhMassDrain > 0.9 * P[i].Mass)
            SphP[i].BhMassDrain = 0.9 * P[i].Mass;
        }
      else
        {
          if(SphP[i].BhMassDrain != 0)
            printf("BH_ACCRETION: WARNING Task %d: cell ID=%llu already has Mass=%g <= 0 before BhMassDrain=%g is applied "
                   "-- corruption upstream of bh_perform_end_of_step_physics(), not from accretion draining\n",
                   ThisTask, (unsigned long long)P[i].ID, P[i].Mass, SphP[i].BhMassDrain);

          SphP[i].BhMassDrain = 0;
        }

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
          // Update total energy 
          SphP[i].Energy += SphP[i].BhThermalFeed;
          All.BhFeedbackLocal[1] += SphP[i].BhThermalFeed;
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

#ifdef BH_JET_FEEDBACK
      if(SphP[i].BhJetEnergy > 0)
        bh_jet_apply_cell(i, &pvd);
#endif
    }
        
  MPI_Allreduce(&All.BhFeedbackLocal, &All.BhFeedbackGlobal, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  mpi_printf("BLACKHOLES: Number of active blackholes = %lld \n", TimeBinsBh.GlobalNActiveParticles);
  mpi_printf("BLACKHOLES: Energy given by BHs = %e, Energy taken up by gas particles = %e \n",
             All.BhFeedbackGlobal[0], All.BhFeedbackGlobal[1]);

#endif
} // perform_end_of_step_physics(void)

#endif /* #ifdef BH_ACTIVE */
