#ifndef BH_H
#define BH_H

#define ALLOC_BH_ROOM 4
extern int NumBhs;

#include "../main/allvars.h"


extern FILE *FdBlackHoles; 

#ifdef BH_ACTIVE
extern struct TimeBinData TimeBinsBh;
#endif

typedef struct Bh_Particle_Data
{
  MyIDType PID;

#ifdef HALO_SEEDING
  MyDouble FormationTime;         /*!< All.Time at which this black hole was seeded */
  MyFloat FormationMetallicity;   /*!< donor cell's own metal mass fraction at seeding; -1 if METALS is off */
  int FormationChannel;           /*!< bitmask of BH_SEED_CHANNEL_* values that triggered seeding (see fof_seeding.h);
                                        int (not a smaller type) so it can be written to snapshots via MEM_INT */
  MyFloat DonorVelocity[3];       /*!< donor cell's velocity at seeding, captured before P[ibh].Vel is set */
#endif

#ifdef BH_ACTIVE
  MyDouble Hsml;
  MyDouble NgbsMass;
  MyDouble NgbsVolume;
  MyDouble AngularMomentum[3];
  int NgbsMinBin;
  int DensityFlag;
  signed char TimeBinBh;
#endif

#ifdef BH_ACCRETION_ACTIVE
  MyDouble Accretion;
  MyDouble MassToDrain;
#endif

#ifdef BONDI_ACCRETION
  MyDouble GasVelocity[3];
  MyDouble GasCircularVelocity[3];
  MyDouble GasDensity;
  MyDouble GasInternalEnergy;
#endif

#ifdef TORQUE_ACCRETION
  MyDouble TorqueMgas;
  MyDouble TorqueMstar;
  MyDouble TorqueMgasDisk;
  MyDouble TorqueMstarDisk;
  MyDouble TorqueFd;
  MyDouble GasAngularMomentum[3];
  MyDouble GasCircularVelocity[3];
#endif

#ifdef ADP_ACCRETION
  MyDouble ADP_CapturedMass;   /* mass captured since last update (code mass) */
  MyDouble ADP_ReservoirMass;  /* reservoir mass waiting to enter disc (code mass) */
  MyDouble ADP_DiscMass;       /* disc mass available to accrete (code mass) */
  MyDouble GasCircularVelocity[3];
#endif

//#ifdef INFALL_ACCRETION
//  MyDouble Accretion;
//#endif
} Bh_Particle_Data; 

extern Bh_Particle_Data *BhP;

#define BPP(i) BhP[P[i].BhID]
#define PPB(i) P[BhP[i].PID]

#endif