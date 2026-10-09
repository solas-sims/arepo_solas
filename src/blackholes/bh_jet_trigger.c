/*! \file bh_jet_trigger.c
 *  \brief Ambient-pressure criterion for launching a black hole jet.
 *
 *  A jet is allowed to launch only if its ram pressure at the launch base exceeds
 *  the thermal pressure of the surrounding gas:
 *
 *      P_ram = rho_jet * v_jet^2  >  BhJetPressureFactor * P_ambient
 *
 *  Jet model (all quantities physical, evaluated in cgs):
 *    - kinetic power      L_jet   = E_step / dt_BH, with E_step from bh_jet_step_energy_cgs():
 *                                   Epsilon_f Epsilon_r Accretion c^2   (default)
 *                                   BhJetPower * dt_BH                  (BH_JET_FIXED_POWER)
 *    - launch speed       v_jet   = BhJetVelocity
 *    - mass outflow rate  Mdot_j  = 2 L_jet / v_jet^2        (purely kinetic jet)
 *    - launch area        A_jet   = pi * (R_inj * sin(theta))^2
 *      with R_inj the BH kernel radius (Hsml) and theta = BhJetOpeningAngle the half-opening angle.
 *    - hence rho_jet = Mdot_j / (A_jet v_jet) and
 *      P_ram = Mdot_j v_jet / A_jet = 2 L_jet / (v_jet A_jet).
 *
 *  P_ambient is the volume-weighted mean of the cell pressures inside the BH kernel,
 *  accumulated in bh_density_evaluate() (BhP[].NgbsPressureVolume / BhP[].NgbsVolume).
 *
 *  Result: BhP[].JetLaunchFlag, consumed by bh_jet_prepare(). With BH_JET_FEEDBACK, energy produced on
 *  a step where the flag is 0 is held in BhP[].JetEnergyReservoir and released at the next launch; the
 *  jet power on release is then larger than the L_jet tested here, so the condition remains satisfied.
 *  With BH_THERMAL_FEEDBACK alone the flag is diagnostic.
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef BH_JET_PRESSURE_TRIGGER

void bh_jet_pressure_trigger(void)
{
  const double v_cgs  = All.BhJetVelocity * 1.0e5;
  const double sin_th = sin(All.BhJetOpeningAngle * M_PI / 180.0);

  long long loc[2] = {0, 0}, glob[2]; /* active BHs, launching BHs */
  double pratio_max = 0;

  for(int idx = 0; idx < TimeBinsBh.NActiveParticles; idx++)
    {
      int i = TimeBinsBh.ActiveParticleList[idx];

      if(i < 0)
        continue;

      loc[0]++;

      BhP[i].AmbientPressure = BhP[i].JetRamPressure = 0;
      BhP[i].JetLaunchFlag = 0;

      double dt_s   = bh_jet_step_time_cgs(i);
      double e_step = bh_jet_step_energy_cgs(i);

      if(e_step <= 0 || dt_s <= 0 || BhP[i].NgbsVolume <= 0)
        continue;

      /* physical ambient pressure [cgs]: SphP.Pressure is stored in comoving code units */
      double p_amb = BhP[i].NgbsPressureVolume / BhP[i].NgbsVolume * All.cf_UnitPressure_in_cgs;

      double r_inj = BhP[i].Hsml * All.cf_UnitLength_in_cm; /* physical length */

      double l_jet = e_step / dt_s;
      double area  = M_PI * r_inj * r_inj * sin_th * sin_th;

      if(v_cgs <= 0 || area <= 0)
        continue;

      double p_ram = 2.0 * l_jet / (v_cgs * area);

      BhP[i].AmbientPressure = p_amb;
      BhP[i].JetRamPressure  = p_ram;

      if(p_ram > All.BhJetPressureFactor * p_amb)
        {
          BhP[i].JetLaunchFlag = 1;
          loc[1]++;
        }

      if(p_amb > 0 && p_ram / p_amb > pratio_max)
        pratio_max = p_ram / p_amb;
    }

  double pratio_max_glob;
  MPI_Allreduce(loc, glob, 2, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&pratio_max, &pratio_max_glob, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  mpi_printf("BLACKHOLES: Jet pressure trigger: %lld of %lld active BHs can launch (max P_ram/P_amb = %g)\n", glob[1],
             glob[0], pratio_max_glob);
}

#endif /* #ifdef BH_JET_PRESSURE_TRIGGER */
