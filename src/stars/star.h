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
/* Lowest stellar mass (in Msolar) which contributes to feedback */
#define LOWEST_MASS_FEEDBACK 2 
/* Lowest stellar mass (in Msolar) which explodes as an SN */
#define LOWEST_MASS_SN 8

#define STAR_MS 0
#define STAR_SN 1
#define STAR_POST_SN 2

#define STAR_UNBORN 0
#define STAR_ACTIVE 1
#define STAR_INACTIVE (-1)  

#define SN_ENERGY 1.0e51

extern struct TimeBinData TimeBinsStar;

typedef struct Star_Interpolate
{
  MyDouble Radius;
  MyDouble Temperature;

#ifdef WINDS
  MyDouble MassLossRate;
#if GRACKLE_CHEMISTRY >= 1
  MyDouble HLossRate;
  MyDouble HeLossRate;
#endif
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
#if GRACKLE_CHEMISTRY >= 1
  MyDouble SN_HLoss;
  MyDouble SN_HeLoss;
#endif
#ifdef METALS
  MyDouble SN_MetalsLoss;
#endif
  MyDouble SN_EnergyInject;
#endif
} Star_Interpolate;

typedef struct Star_Feedback
{
#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  MyDouble TimeToSN;
  MyDouble NextSNEnergy;
#endif

  /* 
    STAR_MS: Main sequence, STAR_SN: Supernova (if any), STAR_POST_SN: After main sequence/Supernova
  */  
  int Stage; 

  /* 1: Active, -1: Inactive */
  int State;

#ifdef WINDS
  MyDouble MassLoss;
#if GRACKLE_CHEMISTRY >= 1
  MyDouble HLoss;
  MyDouble HeLoss;
#endif
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
#if GRACKLE_CHEMISTRY >= 1
  MyDouble SN_HLoss;
  MyDouble SN_HeLoss;
#endif
#ifdef METALS
  MyDouble SN_MetalsLoss;
#endif
  MyDouble SN_EnergyInject;
#endif
} Star_Feedback;

typedef struct Mechanical_Feedback
{
  MyDouble StarPosition[3];
  MyDouble StarVelocity[3];

#ifdef WINDS 
  MyDouble MassLoss;
#if GRACKLE_CHEMISTRY >= 1
  MyDouble HLoss;
  MyDouble HeLoss;
#endif
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
#if GRACKLE_CHEMISTRY >= 1
  MyDouble SN_HLoss;
  MyDouble SN_HeLoss;
#endif
#ifdef METALS
  MyDouble SN_MetalsLoss;
#endif
  MyDouble SN_EnergyInject;
#endif
} Mechanical_Feedback;

typedef struct Mechanical_Feedback_Data
{
  int StarIndex; /* local star index */
  int StarTask; /* task that owns the star */
  int HostIndex; /* local gas-cell index */
  int HostTask; /* task that owns the host */

  Mechanical_Feedback MechanicalFeedback;
} Mechanical_Feedback_Data;

typedef struct Mechanical_Feedback_Events
{
  int NumEvents;
  int MaxEvents;

  long long TotEvents;

  Mechanical_Feedback_Data *MechanicalFeedbackData;
} Mechanical_Feedback_Events;

extern Mechanical_Feedback_Events MechanicalFeedbackEvents;
#endif

typedef struct Star_Particle_Data
{
  MyIDType PID;

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
  MyDouble TimeToSN;
  MyDouble NextSNEnergy;
#endif

#ifdef STAR_FEEDBACK_ACTIVE
  /* 
    Permanent flag
    STAR_UNBORN: before activation, STAR_ACTIVE: post activation, STAR_INACTIVE: inactive
  */  
  int Active; 
  /* Per timestep-> 0: no feedback, 1: feedback (of any type) */
  int WithFeedback; 
  
  MyDouble Hsml;
  int DensityFlag;
  MyDouble NgbsMass;
  MyDouble NgbsVolume;
  int HostHydroBin;
  signed char TimeBinStar;
  
  MyDouble Age;
  MyDouble Birthtime;
  
  Mechanical_Feedback MechanicalFeedback;
#endif
} Star_Particle_Data;

extern Star_Particle_Data *SP;

#define SPP(i) SP[P[i].SID]
#define PPS(i) P[SP[i].PID]

#endif