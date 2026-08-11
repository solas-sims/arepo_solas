#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "../main/allvars.h"
#include "../main/proto.h"


/* Feedback tables interpolation */
static inline double linear_interpolation(double x, double x0, double x1, double y0, double y1);
static inline double star_lifetime(int z_idx, double m_val);
static double lifetime(double z_val, double m_val);

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
static double next_SN_time(double tau, double z_val, double m_val, double a);
#endif

#if defined(WINDS) || defined(STAR_RADIATION_ACTIVE)
static inline Star_Interpolate interpolate_age(int z_idx, int m_idx, double a);
static Star_Interpolate interpolate_mass(int z_idx, double m_val, double a);
static Star_Interpolate interpolate_metallicity(double z_val, double m_val, double a);
#endif

#ifdef SUPERNOVAE
static inline Star_Interpolate SN_interpolate_mass(int z_idx, double m_val);
static Star_Interpolate SN_interpolate_metallicity(double z_val, double m_val);
#endif

/* Linear interpolation helper function */
static inline double linear_interpolation(double x, double x0, double x1, double y0, double y1) 
{
  /* Avoid divide by zero */
  if(x1 == x0) return y0;

  return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

static inline double star_lifetime(int z_idx, double m_val)
{
  /* Lifetime for a given metallicity and mass is just the last entry in the corresponding mass-loss table */
  
  /* Handle mass interpolation */
  if(m_val <= M_VALUES[0]) 
    {
      int N_last = N[z_idx][0];
      return Age[z_idx][0][N_last - 1];
    }
    
  if(m_val >= M_VALUES[M_COUNT - 1]) 
    {
      int N_last = N[z_idx][M_COUNT - 1];
      return Age[z_idx][M_COUNT - 1][N_last - 1];
    }

  const double lm = log10(m_val);

  for(int m = 0; m < M_COUNT - 1; m++) 
    {
      double lm0 = logM_VALUES[m];
      double lm1 = logM_VALUES[m + 1];
      if(lm >= lm0 && lm <= lm1) 
        {
          int N0 = N[z_idx][m];
          int N1 = N[z_idx][m + 1];
          double t0 = Age[z_idx][m][N0 - 1];
          double t1 = Age[z_idx][m + 1][N1 - 1];
          double lt = linear_interpolation(lm, lm0, lm1, log10(t0), log10(t1));

          return pow(10, lt);
        }
    }
  
  terminate("star_lifetime: failed to bracket mass");
}

static double lifetime(double z_val, double m_val)
{
  if(z_val <= Z_VALUES[0])
    return star_lifetime(0, m_val);

  if(z_val >= Z_VALUES[Z_COUNT - 1])
    return star_lifetime(Z_COUNT - 1, m_val);

  const double lz = log10(z_val); 

  for(int z = 0; z < Z_COUNT - 1; z++) 
    {
      double lz0 = logZ_VALUES[z];
      double lz1 = logZ_VALUES[z + 1];
      if(lz >= lz0 && lz <= lz1) 
        {
          double t0 = star_lifetime(z, m_val);
          double t1 = star_lifetime(z + 1, m_val);
          double lt = linear_interpolation(lz, lz0, lz1, log10(t0), log10(t1));

          return pow(10, lt);
        }
    }
  
  terminate("lifetime: failed to bracket metallicity");
}

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
static double next_SN_time(double tau, double z_val, double m_val, double a)
{ 
  Star_Interpolate SN_Feedback = SN_interpolate_metallicity(z_val, m_val);
      
  /* Real SN */
  if(SN_Feedback.SN_MassLoss > 0.0)
    return tau - a;

  /* Failed SN or direct-collapse BH */
  return MAX_REAL_NUMBER;   
}
#endif

#if defined(WINDS) || defined(STAR_RADIATION_ACTIVE)
/* Linear interpolation in age */
static inline Star_Interpolate interpolate_age(int z_idx, int m_idx, double a) 
{
  Star_Interpolate Feedback = {0};
  
  const double *age = Age[z_idx][m_idx];
  const double *radius = Radius[z_idx][m_idx];
  const double *temperature = Temperature[z_idx][m_idx];
  
#ifdef WINDS
  const double *masslossrate  = MassLossRate[z_idx][m_idx];
#if GRACKLE_CHEMISTRY >= 1
  const double *Hlossrate  = HLossRate[z_idx][m_idx];
  const double *Helossrate  = HeLossRate[z_idx][m_idx];
#endif
#ifdef METALS
  const double *metalslossrate = MetalsLossRate[z_idx][m_idx];
#endif
  const double *windvelocity = WindVelocity[z_idx][m_idx];
#endif

#ifdef STAR_RADIATION_ACTIVE
  const WavebandData *flux[WAVEBANDS];
  for(int w = 0; w < WAVEBANDS; w++)
    flux[w] = Flux[w][z_idx][m_idx];
#endif

  int n = N[z_idx][m_idx];

  if(a <= age[0])
    {
      Feedback.Radius = radius[0];
      Feedback.Temperature = temperature[0];

#ifdef WINDS
      Feedback.MassLossRate = masslossrate[0];
#if GRACKLE_CHEMISTRY >= 1
      Feedback.HLossRate = Hlossrate[0];
      Feedback.HeLossRate = Helossrate[0];
#endif
#ifdef METALS
      Feedback.MetalsLossRate = metalslossrate[0];
#endif
      Feedback.WindVelocity = windvelocity[0];
#endif

#ifdef STAR_RADIATION_ACTIVE
      for(int w = 0; w < WAVEBANDS; w++)
        Feedback.Flux[w] = flux[w][0];
#endif
      
      return Feedback;
    }       
    
  if(a >= age[n - 1])  
    {
      Feedback.Radius = radius[n - 1];
      Feedback.Temperature = temperature[n - 1];

#ifdef WINDS
      Feedback.MassLossRate = masslossrate[n - 1];
#if GRACKLE_CHEMISTRY >= 1
      Feedback.HLossRate = Hlossrate[n - 1];
      Feedback.HeLossRate = Helossrate[n - 1];
#endif 
#ifdef METALS
      Feedback.MetalsLossRate = metalslossrate[n - 1];
#endif
      Feedback.WindVelocity = windvelocity[n - 1];
#endif

#ifdef STAR_RADIATION_ACTIVE
      for(int w = 0; w < WAVEBANDS; w++)
        Feedback.Flux[w] = flux[w][n - 1];
#endif
      
      return Feedback;
    } 
  
  for(int i = 0; i < n - 1; i++)
    {
      if(a >= age[i] && a <= age[i + 1])
        {
          Feedback.Radius = linear_interpolation(a, age[i], age[i + 1], radius[i], radius[i + 1]);
          Feedback.Temperature = linear_interpolation(a, age[i], age[i + 1], temperature[i], temperature[i + 1]);

#ifdef WINDS
          Feedback.MassLossRate = linear_interpolation(a, age[i], age[i + 1], masslossrate[i], masslossrate[i + 1]);
#if GRACKLE_CHEMISTRY >= 1
          Feedback.HLossRate = linear_interpolation(a, age[i], age[i + 1], Hlossrate[i], Hlossrate[i + 1]);
          Feedback.HeLossRate = linear_interpolation(a, age[i], age[i + 1], Helossrate[i], Helossrate[i + 1]);
#endif
#ifdef METALS
          Feedback.MetalsLossRate = linear_interpolation(a, age[i], age[i + 1], metalslossrate[i], metalslossrate[i + 1]);
#endif
          Feedback.WindVelocity = linear_interpolation(a, age[i], age[i + 1], windvelocity[i], windvelocity[i + 1]);
#endif

#ifdef STAR_RADIATION_ACTIVE
          for(int w = 0; w < WAVEBANDS; w++)
            {
              Feedback.Flux[w].Energy  = linear_interpolation(a, age[i], age[i+1], flux[w][i].Energy,  flux[w][i+1].Energy);
              Feedback.Flux[w].Photons = linear_interpolation(a, age[i], age[i+1], flux[w][i].Photons, flux[w][i+1].Photons);
            }
#endif

          return Feedback;
        } 
    }
  
  terminate("interpolate_age: failed to bracket age");
}

/* Linear interpolation in mass */
static Star_Interpolate interpolate_mass(int z_idx, double m_val, double a) 
{
  if(m_val <= M_VALUES[0])
    return interpolate_age(z_idx, 0, a);

  if(m_val >= M_VALUES[M_COUNT - 1])
    return interpolate_age(z_idx, M_COUNT - 1, a);

  const double lm = log10(m_val);

  for(int m = 0; m < M_COUNT - 1; m++)
    {
      double lm0 = logM_VALUES[m];
      double lm1 = logM_VALUES[m + 1];
      if(lm >= lm0 && lm <= lm1)
        {
          Star_Interpolate Feedback0 = interpolate_age(z_idx, m, a);
          Star_Interpolate Feedback1 = interpolate_age(z_idx, m + 1, a);
          Star_Interpolate Feedback = {0};

          Feedback.Radius = linear_interpolation(lm, lm0, lm1, Feedback0.Radius, Feedback1.Radius);
          Feedback.Temperature = linear_interpolation(lm, lm0, lm1, Feedback0.Temperature, Feedback1.Temperature);

#ifdef WINDS
          Feedback.MassLossRate = linear_interpolation(lm, lm0, lm1, Feedback0.MassLossRate, Feedback1.MassLossRate);
#if GRACKLE_CHEMISTRY >= 1
          Feedback.HLossRate = linear_interpolation(lm, lm0, lm1, Feedback0.HLossRate, Feedback1.HLossRate);
          Feedback.HeLossRate = linear_interpolation(lm, lm0, lm1, Feedback0.HeLossRate, Feedback1.HeLossRate);
#endif
#ifdef METALS
          Feedback.MetalsLossRate = linear_interpolation(lm, lm0, lm1, Feedback0.MetalsLossRate, Feedback1.MetalsLossRate);
#endif
          Feedback.WindVelocity = linear_interpolation(lm, lm0, lm1, Feedback0.WindVelocity, Feedback1.WindVelocity);
#endif

#ifdef STAR_RADIATION_ACTIVE
          for(int w = 0; w < WAVEBANDS; w++)
            {
              Feedback.Flux[w].Energy = linear_interpolation(lm, lm0, lm1, Feedback0.Flux[w].Energy, Feedback1.Flux[w].Energy);
              Feedback.Flux[w].Photons = linear_interpolation(lm, lm0, lm1, Feedback0.Flux[w].Photons, Feedback1.Flux[w].Photons);
            }
#endif

          return Feedback;
        }
    }
  
  terminate("interpolate_mass: failed to bracket mass");
}

/* Linear interpolation in metallicity */
static Star_Interpolate interpolate_metallicity(double z_val, double m_val, double a)
{
  if(z_val <= Z_VALUES[0])
    return interpolate_mass(0, m_val, a);

  if(z_val >= Z_VALUES[Z_COUNT - 1])
    return interpolate_mass(Z_COUNT - 1, m_val, a);

  const double lz = log10(z_val);

  for(int z = 0; z < Z_COUNT - 1; z++)
    {
      double lz0 = logZ_VALUES[z];
      double lz1 = logZ_VALUES[z + 1];
      if(lz >= lz0 && lz <= lz1)
        {
          Star_Interpolate Feedback0 = interpolate_mass(z, m_val, a);
          Star_Interpolate Feedback1 = interpolate_mass(z + 1, m_val, a);
          Star_Interpolate Feedback = {0};

          Feedback.Radius = linear_interpolation(lz, lz0, lz1, Feedback0.Radius, Feedback1.Radius);
          Feedback.Temperature = linear_interpolation(lz, lz0, lz1, Feedback0.Temperature, Feedback1.Temperature);

#ifdef WINDS
          Feedback.MassLossRate = linear_interpolation(lz, lz0, lz1, Feedback0.MassLossRate, Feedback1.MassLossRate);
#if GRACKLE_CHEMISTRY >= 1
          Feedback.HLossRate = linear_interpolation(lz, lz0, lz1, Feedback0.HLossRate, Feedback1.HLossRate);
          Feedback.HeLossRate = linear_interpolation(lz, lz0, lz1, Feedback0.HeLossRate, Feedback1.HeLossRate);
#endif
#ifdef METALS
          Feedback.MetalsLossRate = linear_interpolation(lz, lz0, lz1, Feedback0.MetalsLossRate, Feedback1.MetalsLossRate);
#endif
          Feedback.WindVelocity = linear_interpolation(lz, lz0, lz1, Feedback0.WindVelocity, Feedback1.WindVelocity);
#endif

#ifdef STAR_RADIATION_ACTIVE
          for(int w = 0; w < WAVEBANDS; w++)
            {
              Feedback.Flux[w].Energy = linear_interpolation(lz, lz0, lz1, Feedback0.Flux[w].Energy, Feedback1.Flux[w].Energy);
              Feedback.Flux[w].Photons = linear_interpolation(lz, lz0, lz1, Feedback0.Flux[w].Photons, Feedback1.Flux[w].Photons);
            }
#endif

          return Feedback;
        }
    }
  
  terminate("interpolate_metallicity: failed to bracket metallicity");
}
#endif

#ifdef SUPERNOVAE
/* Linear interpolation in mass */
static inline Star_Interpolate SN_interpolate_mass(int z_idx, double m_val) 
{
  Star_Interpolate SN_Feedback = {0};
  
  const double *SN_massloss  = SN_MassLoss[z_idx];
#if GRACKLE_CHEMISTRY >= 1
  const double *SN_Hloss  = SN_HLoss[z_idx];
  const double *SN_Heloss  = SN_HeLoss[z_idx];
#endif
#ifdef METALS
  const double *SN_metalsloss = SN_MetalsLoss[z_idx];
#endif

  if(m_val <= M_VALUES[0])
    {
      SN_Feedback.SN_MassLoss = SN_massloss[0];
#if GRACKLE_CHEMISTRY >= 1
      SN_Feedback.SN_HLoss = SN_Hloss[0];
      SN_Feedback.SN_HeLoss = SN_Heloss[0];
#endif
#ifdef METALS
      SN_Feedback.SN_MetalsLoss = SN_metalsloss[0];
#endif
      SN_Feedback.SN_EnergyInject = (SN_Feedback.SN_MassLoss > 0.0) ? SN_ENERGY : 0.0;
      
      return SN_Feedback;
    }       
    
  if(m_val >= M_VALUES[M_COUNT - 1])
    {
      SN_Feedback.SN_MassLoss = SN_massloss[M_COUNT - 1];
#if GRACKLE_CHEMISTRY >= 1
      SN_Feedback.SN_HLoss = SN_Hloss[M_COUNT - 1];
      SN_Feedback.SN_HeLoss = SN_Heloss[M_COUNT - 1];
#endif
#ifdef METALS
      SN_Feedback.SN_MetalsLoss = SN_metalsloss[M_COUNT - 1];
#endif
      SN_Feedback.SN_EnergyInject = (SN_Feedback.SN_MassLoss > 0.0) ? SN_ENERGY : 0.0;
      
      return SN_Feedback;
    } 

  const double lm = log10(m_val);
  
  for(int m = 0; m < M_COUNT - 1; m++)
    {
      double lm0 = logM_VALUES[m];
      double lm1 = logM_VALUES[m + 1];
      if(lm >= lm0 && lm <= lm1)
        {
          /* Both SN -> interpolate */
          if(SN_massloss[m] > 0.0 && SN_massloss[m + 1] > 0.0)
            {
              SN_Feedback.SN_MassLoss = linear_interpolation(lm, lm0, lm1, SN_massloss[m], SN_massloss[m + 1]);
#if GRACKLE_CHEMISTRY >= 1
              SN_Feedback.SN_HLoss = linear_interpolation(lm, lm0, lm1, SN_Hloss[m], SN_Hloss[m + 1]);
              SN_Feedback.SN_HeLoss = linear_interpolation(lm, lm0, lm1, SN_Heloss[m], SN_Heloss[m + 1]);
#endif
#ifdef METALS
              SN_Feedback.SN_MetalsLoss = linear_interpolation(lm, lm0, lm1, SN_metalsloss[m], SN_metalsloss[m + 1]);
#endif
              SN_Feedback.SN_EnergyInject = SN_ENERGY;
            }
          /* At least one failed SN or direct-collapse BH -> clamp to nearest */
          else
            {
              if(lm - lm0 < lm1 - lm)
                {
                  SN_Feedback.SN_MassLoss = SN_massloss[m];
#if GRACKLE_CHEMISTRY >= 1
                  SN_Feedback.SN_HLoss = SN_Hloss[m];
                  SN_Feedback.SN_HeLoss = SN_Heloss[m];
#endif
#ifdef METALS
                  SN_Feedback.SN_MetalsLoss = SN_metalsloss[m];
#endif
                  SN_Feedback.SN_EnergyInject = (SN_Feedback.SN_MassLoss > 0.0) ? SN_ENERGY : 0.0;
                }
              else 
                {
                  SN_Feedback.SN_MassLoss = SN_massloss[m + 1];
#if GRACKLE_CHEMISTRY >= 1
                  SN_Feedback.SN_HLoss = SN_Hloss[m + 1];
                  SN_Feedback.SN_HeLoss = SN_Heloss[m + 1];
#endif
#ifdef METALS
                  SN_Feedback.SN_MetalsLoss = SN_metalsloss[m + 1];
#endif
                  SN_Feedback.SN_EnergyInject = (SN_Feedback.SN_MassLoss > 0.0) ? SN_ENERGY : 0.0;
                }
            }
          
          return SN_Feedback;
        } 
    }
  
  terminate("SN_interpolate_mass: failed to bracket mass");
}

/* Linear interpolation in metallicity */
static Star_Interpolate SN_interpolate_metallicity(double z_val, double m_val)
{
  if(z_val <= Z_VALUES[0])
    return SN_interpolate_mass(0, m_val);
  
  if(z_val >= Z_VALUES[Z_COUNT - 1])
    return SN_interpolate_mass(Z_COUNT - 1, m_val);

  const double lz = log10(z_val);

  for(int z = 0; z < Z_COUNT - 1; z++)
    {
      double lz0 = logZ_VALUES[z];
      double lz1 = logZ_VALUES[z + 1];
      if(lz >= lz0 && lz <= lz1)
        {
          Star_Interpolate SNfeedback0 = SN_interpolate_mass(z, m_val);
          Star_Interpolate SNfeedback1 = SN_interpolate_mass(z + 1, m_val);
          Star_Interpolate SN_Feedback = {0};

          /* Both SN -> interpolate */
          if(SNfeedback0.SN_MassLoss > 0.0 && SNfeedback1.SN_MassLoss > 0.0)
            {
              SN_Feedback.SN_MassLoss = linear_interpolation(lz, lz0, lz1, SNfeedback0.SN_MassLoss, SNfeedback1.SN_MassLoss);
#if GRACKLE_CHEMISTRY >= 1
              SN_Feedback.SN_HLoss = linear_interpolation(lz, lz0, lz1, SNfeedback0.SN_HLoss, SNfeedback1.SN_HLoss);
              SN_Feedback.SN_HeLoss = linear_interpolation(lz, lz0, lz1, SNfeedback0.SN_HeLoss, SNfeedback1.SN_HeLoss);
#endif
#ifdef METALS
              SN_Feedback.SN_MetalsLoss = linear_interpolation(lz, lz0, lz1, SNfeedback0.SN_MetalsLoss, SNfeedback1.SN_MetalsLoss);
#endif
              SN_Feedback.SN_EnergyInject = SN_ENERGY;
            }
          /* At least one failed SN or direct-collapse BH -> clamp to nearest */
          else
            {
              if(lz - lz0 < lz1 - lz)
                SN_Feedback = SNfeedback0;
              else
                SN_Feedback = SNfeedback1;
            }
          
          return SN_Feedback;
        }
    }

  terminate("SN_interpolate_metallicity: failed to bracket metallicity");
}
#endif

/* Wrapper function */
Star_Feedback star_feedback_compute(double dt, double z_val, double m_val, double a)
{
  Star_Feedback Star = {0};

  Star.State = -1;

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  Star.TimeToSN = MAX_REAL_NUMBER;
  Star.NextSNEnergy = 0.0;
#endif

  if(m_val <= LOWEST_MASS_FEEDBACK)
    return Star;

  double tau = lifetime(z_val, m_val);

  if(a > tau)
    {
      Star.Stage = STAR_POST_SN;

      return Star;
    }

  Star.State = 1; 

  if(tau < a + dt)
    {
      Star.Stage = STAR_SN;

      if(m_val < LOWEST_MASS_SN)
        return Star;

#ifdef SUPERNOVAE
      Star_Interpolate SN_Feedback = SN_interpolate_metallicity(z_val, m_val);
      
      Star.SN_MassLoss = SN_Feedback.SN_MassLoss;
#if GRACKLE_CHEMISTRY >= 1
      Star.SN_HLoss = SN_Feedback.SN_HLoss;
      Star.SN_HeLoss = SN_Feedback.SN_HeLoss;
#endif
#ifdef METALS
      Star.SN_MetalsLoss = SN_Feedback.SN_MetalsLoss;
#endif
      Star.SN_EnergyInject = SN_Feedback.SN_EnergyInject;
#endif

      return Star;
    }
  else 
    {
      Star.Stage = STAR_MS;

#if defined(WINDS) || defined(STAR_RADIATION_ACTIVE)
      Star_Interpolate Feedback = interpolate_metallicity(z_val, m_val, a);

#ifdef WINDS
      Star.MassLoss = Feedback.MassLossRate * dt;
#if GRACKLE_CHEMISTRY >= 1
      Star.HLoss = Feedback.HLossRate * dt;
      Star.HeLoss = Feedback.HeLossRate * dt;
#endif
#ifdef METALS
      Star.MetalsLoss = Feedback.MetalsLossRate * dt;
#endif
      Star.WindMomentum = Feedback.MassLossRate * dt * Feedback.WindVelocity;
#endif
      
#ifdef STAR_RADIATION_ACTIVE
      double dt_rad = dt * SEC_PER_YEAR;
      double flux_to_luminosity = 4 * M_PI * Feedback.Radius * Feedback.Radius;
        
      for(int w = 0; w < WAVEBANDS; w++)
        {
          Star.Radiated[w].Energy  = Feedback.Flux[w].Energy  * flux_to_luminosity * dt_rad;
          Star.Radiated[w].Photons = Feedback.Flux[w].Photons * flux_to_luminosity * dt_rad;
        }
#endif

#endif

      if(m_val < LOWEST_MASS_SN)
        return Star;

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
      Star.TimeToSN = next_SN_time(tau, z_val, m_val, a);
      
      if(Star.TimeToSN < MAX_REAL_NUMBER)
        Star.NextSNEnergy = SN_ENERGY;
#endif

      return Star;
    }
}

Star_Feedback units_for_feedback(Star_Feedback StarFeedback)
{

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  if(StarFeedback.TimeToSN < MAX_REAL_NUMBER)
    {
      StarFeedback.TimeToSN /= All.cf_UnitTime_in_yr;
      StarFeedback.NextSNEnergy /= All.cf_UnitEnergy_in_cgs;
    }
#endif

#ifdef WINDS
  StarFeedback.MassLoss /= All.cf_UnitMass_in_Msun;
#if GRACKLE_CHEMISTRY >= 1
  StarFeedback.HLoss /= All.cf_UnitMass_in_Msun;
  StarFeedback.HeLoss /= All.cf_UnitMass_in_Msun;
#endif
#ifdef METALS
  StarFeedback.MetalsLoss /= All.cf_UnitMass_in_Msun;
#endif
  StarFeedback.WindMomentum *= SOLAR_MASS * 1e5 / All.cf_UnitMomentum_in_cgs;
#endif

#ifdef STAR_RADIATION_ACTIVE
  for(int w = 0; w < WAVEBANDS; w++)
    StarFeedback.Radiated[w].Energy /= All.cf_UnitEnergy_in_cgs;
#endif

#ifdef SUPERNOVAE
  StarFeedback.SN_MassLoss /= All.cf_UnitMass_in_Msun;
#if GRACKLE_CHEMISTRY >= 1
  StarFeedback.SN_HLoss /= All.cf_UnitMass_in_Msun;
  StarFeedback.SN_HeLoss /= All.cf_UnitMass_in_Msun;
#endif
#ifdef METALS
  StarFeedback.SN_MetalsLoss /= All.cf_UnitMass_in_Msun;
#endif
  StarFeedback.SN_EnergyInject /= All.cf_UnitEnergy_in_cgs;
#endif

  return StarFeedback;
}

#if defined(STAR_PARTICLES) && STAR_PARTICLES < 2
Star_Feedback star_particle_feedback(int index, double dt, double z, double a)
{  
  int i, Nstars;
  double m;
  Star_Feedback StarParticle = {0};

  StarParticle.State = -1;

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  StarParticle.TimeToSN = MAX_REAL_NUMBER;
  StarParticle.NextSNEnergy = 0.0; 
#endif

  /* Add feedback contributions for each bin */ 
  for(i = 0; i < NBINS; i++) 
    {
      if(SP[index].NumOfStarsInBins[i] == 0)
        continue;

      Nstars = SP[index].NumOfStarsInBins[i];
      
      m = StarMeanMassInBins[i]; 

      Star_Feedback Star = star_feedback_compute(dt, z, m, a);

      if(Star.State < 0)
        continue;

      StarParticle.State = 1;

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
      if(Star.TimeToSN < StarParticle.TimeToSN)
        {
          StarParticle.TimeToSN = Star.TimeToSN;
          StarParticle.NextSNEnergy = Nstars * Star.NextSNEnergy;
        }
#endif

      switch(Star.Stage)
        {
          case STAR_MS:

#ifdef WINDS
          StarParticle.MassLoss += Nstars * Star.MassLoss;
#if GRACKLE_CHEMISTRY >= 1
          StarParticle.HLoss += Nstars * Star.HLoss;
          StarParticle.HeLoss += Nstars * Star.HeLoss;
#endif
#ifdef METALS
          StarParticle.MetalsLoss += Nstars * Star.MetalsLoss;
#endif
          StarParticle.WindMomentum += Nstars * Star.WindMomentum;
#endif

#ifdef STAR_RADIATION_ACTIVE
          for(int w = 0; w < WAVEBANDS; w++)
            {
              StarParticle.Radiated[w].Energy  += Nstars * Star.Radiated[w].Energy;
              StarParticle.Radiated[w].Photons += Nstars * Star.Radiated[w].Photons;
            }
#endif
          break;

          case STAR_SN:

#ifdef SUPERNOVAE
          StarParticle.SN_MassLoss += Nstars * Star.SN_MassLoss;
#if GRACKLE_CHEMISTRY >= 1
          StarParticle.SN_HLoss += Nstars * Star.SN_HLoss;
          StarParticle.SN_HeLoss += Nstars * Star.SN_HeLoss;
#endif
#ifdef METALS
          StarParticle.SN_MetalsLoss += Nstars * Star.SN_MetalsLoss;
#endif
          StarParticle.SN_EnergyInject += Nstars * Star.SN_EnergyInject;   
#endif

          break;

          default:
            terminate("star_particle_feedback: bad stage %d (m=%g, a=%g)", Star.Stage, m, a);
        }
    }

  return StarParticle;
}
#endif