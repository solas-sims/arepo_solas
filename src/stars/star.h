#ifndef STAR_H
#define STAR_H

#include "../main/allvars.h"

#ifdef STAR_PARTICLES
#include "../stars/star_particle.h"
#endif

#ifdef STAR_RADIATION_ACTIVE
#include "../stars/star_radiation.h"
#endif

#ifdef STAR_FEEDBACK_ACTIVE
#include "../stars/star_tables.h"
#endif


#define ALLOC_STAR_ROOM 64
extern int NumStars;

#ifdef STAR_FEEDBACK_ACTIVE
extern struct TimeBinData TimeBinsStar;
#endif

#if defined(WINDS) || defined(SUPERNOVAE)
typedef struct Mechanical_Feedback
{
  MyDouble StarPosition[3];
  MyDouble StarVelocity[3];

#ifdef WINDS 
  MyDouble MassLoss;
#ifdef METALS
  MyDouble MetalsLoss;
#endif
  MyDouble WindMomentum;
#endif

#ifdef SUPERNOVAE
  MyDouble SN_MassLoss;
#ifdef METALS
  MyDouble SN_MetalsLoss;
#endif
  MyDouble SN_EnergyInject;
#endif
} Mechanical_Feedback;

typedef struct Mechanical_Feedback_Data
{
  int StarIndex; /* local star index */
  int HostIndex; /* local gas-cell index */
  int HostTask; /* task that owns the host cell */

  Mechanical_Feedback WindsAndSN;
} Mechanical_Feedback_Data;

typedef struct Mechanical_Feedback_Pack
{
  int NumEvents;
  int MaxEvents;

  Mechanical_Feedback_Data *Data;
} Mechanical_Feedback_Pack;

extern Mechanical_Feedback_Pack MechanicalFeedbackEvents;
#endif

typedef struct Star_Particle_Data
{
  MyIDType PID;

#ifdef STAR_PARTICLES
  int PopulationType;
#endif

#ifdef METALS
  MyDouble Metallicity;
#endif

#ifdef STAR_PARTICLES
  MyDouble MassOfStar;
#endif

#if defined(STAR_PARTICLES) && STAR_PARTICLES < 2
  int NumOfStarsInBins[NBINS];
#endif

#ifdef INDIVIDUAL_STAR_BY_STAR_FORMATION
  MyDouble MassToDrain;
#endif

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  MyDouble TimeSN_yr;
#endif

#ifdef STAR_FEEDBACK_ACTIVE
  int Active;
  MyDouble Hsml;
  MyDouble NgbsMass;
  MyDouble NgbsVolume;
  int HostHydroBin;
  int DensityFlag;
  signed char TimeBinStar;
  MyDouble PhysicalAge_yr;
#endif

#if defined(WINDS) || defined(SUPERNOVAE)
  Mechanical_Feedback WindsAndSN;
#endif

#ifdef STAR_RADIATION_ACTIVE
  WavebandData Radiated[WAVEBANDS];
#endif
} Star_Particle_Data;

extern Star_Particle_Data *SP;

#define SPP(i) SP[P[i].SID]
#define PPS(i) P[SP[i].PID]

#ifdef STAR_FEEDBACK_ACTIVE
typedef struct Star_Interpolate
{
  MyDouble Radius;
  MyDouble Temperature;

#ifdef WINDS
  MyDouble MassLossRate;
#ifdef METALS
  MyDouble MetalsLossRate;
#endif
  MyDouble WindVelocity;
#endif

#ifdef STAR_RADIATION_ACTIVE
  WavebandData Flux[WAVEBANDS];
#endif

#ifdef SUPERNOVAE
  MyDouble SN_MassLoss;
#ifdef METALS
  MyDouble SN_MetalsLoss;
#endif
  MyDouble SN_EnergyInject;
#endif
} Star_Interpolate;

typedef struct Star_Feedback
{
#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  double TimeSN;
#endif
  
  int Stage; // 0:preSN, 1:SN, 2:postSN

#ifdef WINDS
  MyDouble MassLoss;
#ifdef METALS
  MyDouble MetalsLoss;
#endif
  MyDouble WindMomentum;
#endif

#ifdef STAR_RADIATION_ACTIVE
  WavebandData Radiated[WAVEBANDS];
#endif

#ifdef SUPERNOVAE
  MyDouble SN_MassLoss;
#ifdef METALS
  MyDouble SN_MetalsLoss;
#endif
  MyDouble SN_EnergyInject;
#endif
} Star_Feedback;
#endif

#endif