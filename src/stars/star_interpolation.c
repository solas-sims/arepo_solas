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
  // avoid divide by zero
  if(x1 == x0) return y0;

  return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

static inline double star_lifetime(int z_idx, double m_val)
{
  // Lifetime for a given metallicity and mass is just the last entry
  // in the corresponding mass-loss table
  
  // Handle mass interpolation
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

  for(int m = 0; m < M_COUNT - 1; m++) 
    {
      double m0 = M_VALUES[m];
      double m1 = M_VALUES[m + 1];
      if(m_val >= m0 && m_val <= m1) 
        {
          int N0 = N[z_idx][m];
          int N1 = N[z_idx][m + 1];
          double t0 = Age[z_idx][m][N0 - 1];
          double t1 = Age[z_idx][m + 1][N1 - 1];
          return linear_interpolation(m_val, m0, m1, t0, t1);
        }
    }
  terminate("Star_lifetime: failed to bracket mass");
}

static double lifetime(double z_val, double m_val)
{
  if(z_val <= Z_VALUES[0])
    return star_lifetime(0, m_val);
  if(z_val >= Z_VALUES[Z_COUNT - 1])
    return star_lifetime(Z_COUNT - 1, m_val);

  for(int z = 0; z < Z_COUNT - 1; z++) 
    {
      double z0 = Z_VALUES[z];
      double z1 = Z_VALUES[z + 1];
      if(z_val >= z0 && z_val <= z1) 
        {
          double t0 = star_lifetime(z, m_val);
          double t1 = star_lifetime(z + 1, m_val);
          return linear_interpolation(z_val, z0, z1, t0, t1);
        }
    }
  terminate("Lifetime: failed to bracket metallicity");
}

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
static double next_SN_time(double tau, double z_val, double m_val, double a)
{ 
  if(m_val >= 8 && tau > a)
    {
      Star_Interpolate SN_Feedback = SN_interpolate_metallicity(z_val, m_val);
      if(SN_Feedback.SN_MassLoss > 0.0)
        return tau;
    }
  return MAX_REAL_NUMBER; /* No SN or already past SN */
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
  terminate("Interpolate_age: failed to bracket age");
}

/* Linear interpolation in mass */
static Star_Interpolate interpolate_mass(int z_idx, double m_val, double a) 
{
  if(m_val <= M_VALUES[0])
    return interpolate_age(z_idx, 0, a);
  if(m_val >= M_VALUES[M_COUNT - 1])
    return interpolate_age(z_idx, M_COUNT - 1, a);

  for(int m = 0; m < M_COUNT - 1; m++)
    {
    double m0 = M_VALUES[m];
      double m1 = M_VALUES[m + 1];
      if(m_val >= m0 && m_val <= m1)
        {
          Star_Interpolate Feedback0 = interpolate_age(z_idx, m, a);
          Star_Interpolate Feedback1 = interpolate_age(z_idx, m + 1, a);
          Star_Interpolate Feedback = {0};

          Feedback.Radius = linear_interpolation(m_val, m0, m1, Feedback0.Radius, Feedback1.Radius);
          Feedback.Temperature = linear_interpolation(m_val, m0, m1, Feedback0.Temperature, Feedback1.Temperature);

#ifdef WINDS
          Feedback.MassLossRate = linear_interpolation(m_val, m0, m1, Feedback0.MassLossRate, Feedback1.MassLossRate);
#ifdef METALS
          Feedback.MetalsLossRate = linear_interpolation(m_val, m0, m1, Feedback0.MetalsLossRate, Feedback1.MetalsLossRate);
#endif
          Feedback.WindVelocity = linear_interpolation(m_val, m0, m1, Feedback0.WindVelocity, Feedback1.WindVelocity);
#endif

#ifdef STAR_RADIATION_ACTIVE
          for(int w = 0; w < WAVEBANDS; w++)
            {
              Feedback.Flux[w].Energy = linear_interpolation(m_val, m0, m1, Feedback0.Flux[w].Energy, Feedback1.Flux[w].Energy);
              Feedback.Flux[w].Photons = linear_interpolation(m_val, m0, m1, Feedback0.Flux[w].Photons, Feedback1.Flux[w].Photons);
            }
#endif

          return Feedback;
        }
    }
  terminate("Interpolate_mass: failed to bracket mass");
}

/* Linear interpolation in metallicity */
static Star_Interpolate interpolate_metallicity(double z_val, double m_val, double a)
{
  if(z_val <= Z_VALUES[0])
    return interpolate_mass(0, m_val, a);
  if(z_val >= Z_VALUES[Z_COUNT - 1])
    return interpolate_mass(Z_COUNT - 1, m_val, a);

  for(int z = 0; z < Z_COUNT - 1; z++)
    {
      double z0 = Z_VALUES[z];
      double z1 = Z_VALUES[z + 1];
      if(z_val >= z0 && z_val <= z1)
        {
          Star_Interpolate Feedback0 = interpolate_mass(z, m_val, a);
          Star_Interpolate Feedback1 = interpolate_mass(z + 1, m_val, a);
          Star_Interpolate Feedback = {0};

          Feedback.Radius = linear_interpolation(z_val, z0, z1, Feedback0.Radius, Feedback1.Radius);
          Feedback.Temperature = linear_interpolation(z_val, z0, z1, Feedback0.Temperature, Feedback1.Temperature);

#ifdef WINDS
          Feedback.MassLossRate = linear_interpolation(z_val, z0, z1, Feedback0.MassLossRate, Feedback1.MassLossRate);
#ifdef METALS
          Feedback.MetalsLossRate = linear_interpolation(z_val, z0, z1, Feedback0.MetalsLossRate, Feedback1.MetalsLossRate);
#endif
          Feedback.WindVelocity = linear_interpolation(z_val, z0, z1, Feedback0.WindVelocity, Feedback1.WindVelocity);
#endif

#ifdef STAR_RADIATION_ACTIVE
          for(int w = 0; w < WAVEBANDS; w++)
            {
              Feedback.Flux[w].Energy = linear_interpolation(z_val, z0, z1, Feedback0.Flux[w].Energy, Feedback1.Flux[w].Energy);
              Feedback.Flux[w].Photons = linear_interpolation(z_val, z0, z1, Feedback0.Flux[w].Photons, Feedback1.Flux[w].Photons);
            }
#endif

          return Feedback;
        }
    }
  terminate("Interpolate_metallicity: failed to bracket metallicity");
}
#endif

#ifdef SUPERNOVAE
/* Linear interpolation in mass */
static inline Star_Interpolate SN_interpolate_mass(int z_idx, double m_val) 
{
  Star_Interpolate SN_Feedback = {0};
  
  const double *SN_massloss  = SN_MassLoss[z_idx];
#ifdef METALS
  const double *SN_metalsloss = SN_MetalsLoss[z_idx];
#endif

  if(m_val <= M_VALUES[0])
    {
      SN_Feedback.SN_MassLoss = SN_massloss[0];
#ifdef METALS
      SN_Feedback.SN_MetalsLoss = SN_metalsloss[0];
#endif
      SN_Feedback.SN_EnergyInject = (SN_Feedback.SN_MassLoss > 0.0) ? 1e51 : 0.0;
      
      return SN_Feedback;
    }       
    
  if(m_val >= M_VALUES[M_COUNT - 1])
    {
      SN_Feedback.SN_MassLoss = SN_massloss[M_COUNT - 1];
#ifdef METALS
      SN_Feedback.SN_MetalsLoss = SN_metalsloss[M_COUNT - 1];
#endif
      SN_Feedback.SN_EnergyInject = (SN_Feedback.SN_MassLoss > 0.0) ? 1e51 : 0.0;
      
      return SN_Feedback;
    } 
  
  for(int m = 0; m < M_COUNT - 1; m++)
    {
      double m0 = M_VALUES[m];
      double m1 = M_VALUES[m + 1];
      if(m_val >= m0 && m_val <= m1)
        {
          if(SN_massloss[m] > 0 && SN_massloss[m + 1] > 0)
            {
              SN_Feedback.SN_MassLoss = linear_interpolation(m_val, m0, m1, SN_massloss[m], SN_massloss[m + 1]);
#ifdef METALS
              SN_Feedback.SN_MetalsLoss = linear_interpolation(m_val, m0, m1, SN_metalsloss[m], SN_metalsloss[m + 1]);
#endif
              SN_Feedback.SN_EnergyInject = (SN_Feedback.SN_MassLoss > 0.0) ? 1e51 : 0.0;
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

  for(int z = 0; z < Z_COUNT - 1; z++)
    {
      double z0 = Z_VALUES[z];
      double z1 = Z_VALUES[z + 1];
      if(z_val >= z0 && z_val <= z1)
        {
          Star_Interpolate SNfeedback0 = SN_interpolate_mass(z, m_val);
          Star_Interpolate SNfeedback1 = SN_interpolate_mass(z + 1, m_val);
          Star_Interpolate SN_Feedback = {0};

          if(SNfeedback0.SN_MassLoss > 0.0 && SNfeedback1.SN_MassLoss > 0.0)
            {
              SN_Feedback.SN_MassLoss = linear_interpolation(z_val, z0, z1, SNfeedback0.SN_MassLoss, SNfeedback1.SN_MassLoss);
#ifdef METALS
              SN_Feedback.SN_MetalsLoss = linear_interpolation(z_val, z0, z1, SNfeedback0.SN_MetalsLoss, SNfeedback1.SN_MetalsLoss);
#endif
              SN_Feedback.SN_EnergyInject = (SN_Feedback.SN_MassLoss > 0.0) ? 1e51 : 0.0;
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
  double tau = lifetime(z_val, m_val);
  Star_Feedback Star = {0};

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  Star.TimeSN = MAX_REAL_NUMBER;
#endif

  if(m_val <= LOWEST_MASS_FEEDBACK)
    return Star;

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  Star.TimeSN = next_SN_time(tau, z_val, m_val, a);
#endif

  if(a < tau)
    {
      Star.Stage = 0;

#if defined(WINDS) || defined(STAR_RADIATION_ACTIVE)
      Star_Interpolate Feedback = interpolate_metallicity(z_val, m_val, a);

#ifdef WINDS
      Star.MassLoss = Feedback.MassLossRate * dt;
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

      return Star;
    }

  if(a >= tau && a < (tau+dt)) 
    {
      Star.Stage = 1;

#ifdef SUPERNOVAE
      Star_Interpolate SN_Feedback = SN_interpolate_metallicity(z_val, m_val);
      Star.SN_MassLoss = SN_Feedback.SN_MassLoss;
#ifdef METALS
      Star.SN_MetalsLoss = SN_Feedback.SN_MetalsLoss;
#endif
      Star.SN_EnergyInject = SN_Feedback.SN_EnergyInject;
#endif

      return Star;
    }
  
  if(a >= tau+dt)
    {
      Star.Stage = 2;

      return Star;
    }
  return Star;
}

Star_Feedback units_for_feedback(Star_Feedback StarFeedback)
{
#ifdef WINDS
  StarFeedback.MassLoss /= All.cf_UnitMass_in_Msun;
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

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  StarParticle.TimeSN = MAX_REAL_NUMBER; 
#endif

  // Add feedback contributions for each bin 
  for(i = 0; i < NBINS; i++) 
    {
      if(SP[index].NumOfStarsInBins[i] == 0)
        continue;

      Nstars = SP[index].NumOfStarsInBins[i];
      
      m = StarMeanMassInBins[i]; 

      Star_Feedback Star = star_feedback_compute(dt, z, m, a);

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
      if(Star.TimeSN < StarParticle.TimeSN)
      StarParticle.TimeSN = Star.TimeSN;
#endif

      switch(Star.Stage)
        {
          case 0:

#ifdef WINDS
          StarParticle.MassLoss += Nstars * Star.MassLoss;
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

          case 1:

#ifdef SUPERNOVAE
          StarParticle.SN_MassLoss += Nstars * Star.SN_MassLoss;
#ifdef METALS
          StarParticle.SN_MetalsLoss += Nstars * Star.SN_MetalsLoss;
#endif
          StarParticle.SN_EnergyInject += Nstars * Star.SN_EnergyInject;   
#endif

          break;

          case 2:

          break;
        }
    }
  return StarParticle;
}
#endif