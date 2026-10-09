/*! \file bh_jet_feedback.c
 *  \brief Bipolar jet injection from black holes into the surrounding gas cells.
 *
 *  Each active BH with BhP[].JetEnergyRelease > 0 (see bh_jet_prepare()) deposits that kinetic energy in
 *  the gas cells that lie within its kernel radius and within BhJetOpeningAngle of its jet axis (either
 *  lobe). Cell shares are proportional to cell mass (m_j / BhP[].JetConeMass). The kick on a cell is
 *  directed radially away from the BH, so the two lobes carry equal and opposite momentum.
 *
 *  This routine only fills per-cell buffers (SphP[].BhJet*). They are applied by bh_perform_end_of_step_physics()
 *  when the cell is next active, in one of three ways (compile-time, one of):
 *    MOMENTUM_JET  add momentum 2 E / v_jet; thermal energy of the cell is unchanged
 *    THERMAL_JET   add the same momentum and add E to the total energy; the kinetic part of the cell's
 *                  energy gain comes from the momentum, the remainder is thermal
 *    ENERGY_JET    add E to the total energy and kick along the cone direction by exactly the momentum
 *                  that puts all of E into the kinetic energy of the cell
 *
 *  BH_JET_FIXED_POWER also creates the jet mass 2 E / v_jet^2, at rest and cold (arepobh behaviour); otherwise the jet
 *  mass is entrained from the cone cells themselves and the net cell mass is unchanged.
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"

#ifdef BH_JET_FEEDBACK

static int bh_jet_feedback_evaluate(int target, int mode, int threadid);

typedef struct
{
  MyDouble Pos[3];
  MyDouble Axis[3];
  MyDouble Energy;
  MyDouble ConeMass;
  MyFloat Hsml;
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

static void particle2in(data_in *in, int i, int firstnode)
{
  for(int j = 0; j < 3; j++)
    {
      in->Pos[j]  = PPB(i).Pos[j];
      in->Axis[j] = BhP[i].JetAxis[j];
    }
  in->Energy   = BhP[i].JetEnergyRelease;
  in->ConeMass = BhP[i].JetConeMass;
  in->Hsml     = BhP[i].Hsml;
  in->Firstnode = firstnode;
}

typedef struct
{
} data_out;

static data_out *DataResult, *DataOut;

static void out2particle(data_out *out, int i, int mode)
{
}

#include "../utils/generic_comm_helpers2.h"

static void kernel_local(void)
{
  int i, idx;
  int j, threadid = get_thread_num();

  for(j = 0; j < NTask; j++)
    Thread[threadid].Exportflag[j] = -1;

  while(1)
    {
      if(Thread[threadid].ExportSpace < MinSpace)
        break;

      idx = NextParticle++;

      if(idx >= TimeBinsBh.NActiveParticles)
        break;

      i = TimeBinsBh.ActiveParticleList[idx];

      if(i < 0)
        continue;

      bh_jet_feedback_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
    }
}

static void kernel_imported(void)
{
  int i, cnt = 0;
  int threadid = get_thread_num();

  while(1)
    {
      i = cnt++;

      if(i >= Nimport)
        break;

      bh_jet_feedback_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
    }
}

void bh_jet_feedback(void)
{
  TIMER_START(CPU_BLACKHOLES_FEEDBACK);

  /* energy released into a cone with no gas stays in the reservoir */
  for(int idx = 0; idx < TimeBinsBh.NActiveParticles; idx++)
    {
      int i = TimeBinsBh.ActiveParticleList[idx];

      if(i < 0)
        continue;

      if(BhP[i].JetEnergyRelease > 0 && BhP[i].JetConeMass <= 0)
        {
          BhP[i].JetEnergyReservoir += BhP[i].JetEnergyRelease;
          BhP[i].JetEnergyRelease = 0;
        }
    }

  generic_set_MaxNexport();

  generic_comm_pattern(TimeBinsBh.NActiveParticles, kernel_local, kernel_imported);

  TIMER_STOP(CPU_BLACKHOLES_FEEDBACK);
}

static int bh_jet_feedback_evaluate(int target, int mode, int threadid)
{
  int j, n, numnodes, *firstnode;
  MyDouble xtmp, ytmp, ztmp;
  double d[3], r;

  data_in local, *target_data;
  data_out out;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      target_data = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  if(target_data->Energy > 0 && target_data->ConeMass > 0)
    {
      const double h      = target_data->Hsml;
      const double h2     = h * h;
      const double cos_th = cos(All.BhJetOpeningAngle * M_PI / 180.0);
      MyDouble *pos       = target_data->Pos;

      /* jet speed in cgs; momentum and mass follow from the cell's energy share */
#if !defined(ENERGY_JET) || defined(BH_JET_FIXED_POWER)
      const double v_cgs = All.BhJetVelocity * 1.0e5;
#endif
#ifndef ENERGY_JET
      const double p_per_e = 2.0 * All.cf_UnitEnergy_in_cgs / (v_cgs * All.cf_UnitMomentum_in_cgs); /* code momentum per code energy */
#endif
#ifdef BH_JET_FIXED_POWER
      const double m_per_e = 2.0 * All.cf_UnitEnergy_in_cgs / (v_cgs * v_cgs * All.cf_UnitMass_in_g);
#endif

      int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

      for(n = 0; n < nfound; n++)
        {
          j = Thread[threadid].Ngblist[n];

          if(P[j].Type != 0 || P[j].Mass == 0 || P[j].ID == 0)
            continue;

          d[0] = NEAREST_X(P[j].Pos[0] - pos[0]);
          d[1] = NEAREST_Y(P[j].Pos[1] - pos[1]);
          d[2] = NEAREST_Z(P[j].Pos[2] - pos[2]);

          if(d[0] * d[0] + d[1] * d[1] + d[2] * d[2] >= h2)
            continue;

          if(!bh_jet_in_cone(d, target_data->Axis, cos_th, &r))
            continue;

          double ej = target_data->Energy * P[j].Mass / target_data->ConeMass;

          SphP[j].BhJetEnergy += ej;

#ifdef ENERGY_JET
          /* direction accumulator, weighted by the energy each BH deposits */
          for(int k = 0; k < 3; k++)
            SphP[j].BhJetMomentum[k] += ej * d[k] / r;
#else
          for(int k = 0; k < 3; k++)
            SphP[j].BhJetMomentum[k] += ej * p_per_e * d[k] / r;
#endif
#ifdef BH_JET_FIXED_POWER
          SphP[j].BhJetMass += ej * m_per_e;
#endif
          All.BhFeedbackLocal[0] += ej;
        }
    }

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif /* #ifdef BH_JET_FEEDBACK */
