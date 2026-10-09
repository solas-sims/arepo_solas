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
#ifdef BH_JET_PRESSURE_TRIGGER
  MyDouble NgbsPressureVolume;    /*!< sum over kernel cells of Pressure * Volume (comoving code units) */
#endif
  int DensityFlag;
  signed char TimeBinBh;
#endif

#ifdef BH_ACCRETION_ACTIVE
  MyDouble Accretion;
  MyDouble MassToDrain;
#endif

#ifdef BH_JET_PRESSURE_TRIGGER
  MyDouble AmbientPressure;       /*!< volume-weighted mean physical gas pressure within the BH kernel [cgs] */
  MyDouble JetRamPressure;        /*!< physical ram pressure of the jet at its launch base [cgs] */
  int JetLaunchFlag;              /*!< 1 if JetRamPressure > BhJetPressureFactor * AmbientPressure this step, else 0 */
#endif

#ifdef BH_JET_FEEDBACK
  MyDouble JetEnergyReservoir;    /*!< jet energy produced but not yet released (code units); carried over while JetLaunchFlag = 0 */
  MyDouble JetEnergyRelease;      /*!< jet energy to inject on this BH step (code units); 0 if no launch */
  MyDouble JetConeMass;           /*!< total mass of the gas cells inside the jet cones within the BH kernel */
  MyDouble JetAxis[3];            /*!< unit vector along the jet axis (either lobe) used on this step */
#ifdef BH_JET_SPIN
  MyDouble JetGasL[3];            /*!< mass-weighted gas angular momentum (relative to the BH) within the kernel */
  MyDouble SpinVec[3];            /*!< spin proxy J: BH mass along the initial axis plus accreted mass times the unit
                                       gas angular-momentum direction of each accretion event; JetAxis = J/|J| */
#endif
#endif

#ifdef BH_MERGER
  int NumMergers; /*!< number of other black holes this one has absorbed via bh_merger() */
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

#ifdef BH_MERGER
/*! \brief Which length scale bh_merger() uses for its "close enough to merge" radius test.
 *  See bh_merger_radius() in bh_merger.c. Values are assigned from the BhMergerRadiusCriterion
 *  string parameter by check_parameters() in src/io/parameters.c -- do not renumber without
 *  updating that resolution logic. */
enum bh_merger_radius_criterion
{
  BH_MERGER_RADIUS_HSML = 0,           /*!< factor * max(Hsml_i, Hsml_j) -- legacy/default behaviour */
  BH_MERGER_RADIUS_SOFTENING,          /*!< factor * max(softening_i, softening_j) */
  BH_MERGER_RADIUS_MAX_HSML_SOFTENING, /*!< factor * max(Hsml_i, Hsml_j, softening_i, softening_j) */
  BH_MERGER_RADIUS_MIN_HSML_SOFTENING, /*!< factor * min(max(Hsml_i,Hsml_j), max(softening_i,softening_j)) */
  BH_MERGER_RADIUS_OTHER_PHYSICAL      /*!< reserved for a future physically-motivated criterion (e.g. mutual
                                             Bondi radius); not implemented -- selecting it terminates the run */
};
#endif /* #ifdef BH_MERGER */

#endif