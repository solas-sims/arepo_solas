/*! \file bh_jet_common.c
 *  \brief Jet bookkeeping shared by the launch trigger and the injection: jet power, energy
 *         reservoir, jet axis (fixed or following the BH spin proxy) and the release decision.
 *
 *  Energy bookkeeping (all per BH step):
 *    produced  E_new  = Epsilon_f Epsilon_r Accretion c^2          (accretion-driven, default)
 *                     = BhJetPower * dt                            (BH_JET_FIXED_POWER)
 *    reservoir R     += E_new
 *    released  E_rel  = R                                          if the launch condition holds
 *                     = n * E_burst, n = floor(R / E_burst)        (BURST_MODE, E_burst = BhJetBurstMass v_jet^2)
 *                     = 0                                          otherwise, and R is carried over
 *  The launch condition is BhP[].JetLaunchFlag (BH_JET_PRESSURE_TRIGGER) or always true without it.
 *  Energy that finds no gas in the jet cone is returned to R by bh_jet_feedback().
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#if defined(BH_JET_FEEDBACK) || defined(BH_JET_PRESSURE_TRIGGER)

double bh_jet_step_time_cgs(int i)
{
  /* All.Timebase_interval is in ln(a) for comoving runs; cf_UnitTime_in_s converts either to seconds */
  double dt_code = (BhP[i].TimeBinBh ? (((integertime)1) << BhP[i].TimeBinBh) : 0) * All.Timebase_interval;

  return dt_code * All.cf_UnitTime_in_s;
}

double bh_jet_step_energy_cgs(int i)
{
#ifdef BH_JET_FIXED_POWER
  return All.BhJetPower * bh_jet_step_time_cgs(i);
#else
  if(BhP[i].Accretion <= 0)
    return 0;

  /* Accretion is the code mass swallowed this step; the radiated fraction is Epsilon_r of it */
  return All.Epsilon_f * All.Epsilon_r * BhP[i].Accretion * All.cf_UnitMass_in_g * CLIGHT * CLIGHT;
#endif
}

#endif /* #if defined(BH_JET_FEEDBACK) || defined(BH_JET_PRESSURE_TRIGGER) */

#ifdef BH_JET_FEEDBACK

void bh_jet_prepare(void)
{
#ifdef BURST_MODE
  const double v_cgs  = All.BhJetVelocity * 1.0e5;
  const double e_burst = All.BhJetBurstMass * SOLAR_MASS * v_cgs * v_cgs / All.cf_UnitEnergy_in_cgs; /* code energy */

  if(e_burst <= 0)
    terminate("BURST_MODE needs BhJetBurstMass > 0 and BhJetVelocity > 0");
#endif

  double sums[3] = {0, 0, 0}, glob[3]; /* produced, released, reservoir (code energy) */
  long long cnt[2] = {0, 0}, cnt_glob[2]; /* active BHs, BHs releasing */

  for(int idx = 0; idx < TimeBinsBh.NActiveParticles; idx++)
    {
      int i = TimeBinsBh.ActiveParticleList[idx];

      if(i < 0)
        continue;

      cnt[0]++;

      /* --- jet axis --- */
#ifdef BH_JET_SPIN
      double *J = BhP[i].SpinVec;

      if(J[0] == 0 && J[1] == 0 && J[2] == 0)
        for(int k = 0; k < 3; k++)
          J[k] = PPB(i).Mass * All.BhJetAxis[k];

      double lnorm = sqrt(BhP[i].JetGasL[0] * BhP[i].JetGasL[0] + BhP[i].JetGasL[1] * BhP[i].JetGasL[1] +
                          BhP[i].JetGasL[2] * BhP[i].JetGasL[2]);

      if(lnorm > 0 && BhP[i].Accretion > 0)
        for(int k = 0; k < 3; k++)
          J[k] += BhP[i].Accretion * BhP[i].JetGasL[k] / lnorm;

      double jnorm = sqrt(J[0] * J[0] + J[1] * J[1] + J[2] * J[2]);

      for(int k = 0; k < 3; k++)
        BhP[i].JetAxis[k] = J[k] / jnorm;
#else
      for(int k = 0; k < 3; k++)
        BhP[i].JetAxis[k] = All.BhJetAxis[k];
#endif

      /* --- energy: produce, then decide what to release --- */
      double e_new = bh_jet_step_energy_cgs(i) / All.cf_UnitEnergy_in_cgs;

      BhP[i].JetEnergyReservoir += e_new;
      BhP[i].JetEnergyRelease = 0;
      BhP[i].JetConeMass      = 0;

      int launch = 1;
#ifdef BH_JET_PRESSURE_TRIGGER
      launch = BhP[i].JetLaunchFlag;
#endif

      double release = 0;

      if(launch && BhP[i].JetEnergyReservoir > 0)
        {
#ifdef BURST_MODE
          double nb = floor(BhP[i].JetEnergyReservoir / e_burst);
          release   = nb * e_burst;
#else
          release = BhP[i].JetEnergyReservoir;
#endif
        }

      BhP[i].JetEnergyRelease = release;
      BhP[i].JetEnergyReservoir -= release;
      if(BhP[i].JetEnergyReservoir < 0)
        BhP[i].JetEnergyReservoir = 0;

      sums[0] += e_new;
      sums[1] += release;
      sums[2] += BhP[i].JetEnergyReservoir;
      if(release > 0)
        cnt[1]++;
    }

  MPI_Allreduce(sums, glob, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(cnt, cnt_glob, 2, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

  mpi_printf("BLACKHOLES: Jet energy [code]: produced %g, released %g by %lld of %lld active BHs, held in reservoirs %g\n", glob[0], glob[1],
             cnt_glob[1], cnt_glob[0], glob[2]);
}

#endif /* #ifdef BH_JET_FEEDBACK */
