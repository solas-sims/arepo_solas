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

#ifdef BH_JET_LAUNCH_POINT

/* splitmix64: small, stateless-per-call generator so the draws depend only on (BH ID, burst number) and are
 * identical across reruns, restarts and task counts */
static double jet_uniform(unsigned long long *state)
{
  unsigned long long z = (*state += 0x9E3779B97F4A7C15ULL);

  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  z = z ^ (z >> 31);

  return (z >> 11) * (1.0 / 9007199254740992.0);
}

/*! \brief Puts a burst of energy e (code units) on its way: random direction within BhJetOpeningAngle of the jet axis,
 *         uniform in cos(theta), arriving at the launch radius after the flight time R / v_jet. */
static void bh_jet_launch_packet(int i, double e)
{
  unsigned long long state = 0x2545F4914F6CDD1DULL ^ ((unsigned long long)PPB(i).ID * 0x9E3779B97F4A7C15ULL) ^
                             ((unsigned long long)BhP[i].JetBurstCount * 0xBF58476D1CE4E5B9ULL);

  double cos_max = cos(All.BhJetOpeningAngle * M_PI / 180.0);
  double cth     = 1.0 - jet_uniform(&state) * (1.0 - cos_max);
  double sth     = sqrt(1.0 - cth * cth);
  double phi     = 2.0 * M_PI * jet_uniform(&state);

  const double *a = BhP[i].JetAxis;

  /* orthonormal frame (e1, e2, a) */
  double h[3] = {0, 0, 0};
  h[fabs(a[0]) < 0.9 ? 0 : 1] = 1.0;

  double e1[3] = {h[1] * a[2] - h[2] * a[1], h[2] * a[0] - h[0] * a[2], h[0] * a[1] - h[1] * a[0]};
  double n1    = sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);

  for(int k = 0; k < 3; k++)
    e1[k] /= n1;

  double e2[3] = {a[1] * e1[2] - a[2] * e1[1], a[2] * e1[0] - a[0] * e1[2], a[0] * e1[1] - a[1] * e1[0]};

  double dir[3];
  for(int k = 0; k < 3; k++)
    dir[k] = cth * a[k] + sth * (cos(phi) * e1[k] + sin(phi) * e2[k]);

  /* flight time in code time (ln a for comoving runs) -> integer time */
  double t_fly_s    = All.BhJetLaunchRadius * 3.0856775814913673e21 / (All.BhJetVelocity * 1.0e5);
  double t_fly_code = t_fly_s / All.cf_UnitTime_in_s;
  integertime ti_arrive = All.Ti_Current + (integertime)(t_fly_code / All.Timebase_interval);

  BhP[i].JetBurstCount++;

  int n = BhP[i].JetNPackets;

  if(n >= BH_JET_MAXPACKETS)
    {
      /* buffer full: fold the burst into the newest packet (energy adds, direction and arrival time are kept) */
      BhP[i].JetPacket[n - 1].Energy += e;
      return;
    }

  BhP[i].JetPacket[n].Energy = e;
  for(int k = 0; k < 3; k++)
    BhP[i].JetPacket[n].Dir[k] = dir[k];
  BhP[i].JetPacket[n].TiArrive = ti_arrive;
  BhP[i].JetPacket[n].New      = 1;
  BhP[i].JetNPackets           = n + 1;
}

#endif /* #ifdef BH_JET_LAUNCH_POINT */

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

      BhP[i].JetEnergyReservoir -= release;
      if(BhP[i].JetEnergyReservoir < 0)
        BhP[i].JetEnergyReservoir = 0;

#ifdef BH_JET_LAUNCH_POINT
      /* the released energy is now a packet in flight; bh_jet_points() delivers it */
      if(release > 0)
        bh_jet_launch_packet(i, release);
      BhP[i].JetEnergyRelease = 0;
#else
      BhP[i].JetEnergyRelease = release;
#endif

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
