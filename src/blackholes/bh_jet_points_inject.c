/*! \file bh_jet_points_inject.c
 *  \brief Second pass of the launch-point jet (BH_JET_LAUNCH_POINT): fills the per-cell jet buffers and marks the
 *         cells for the timestep limiter. See bh_jet_points.c for the geometry and bh_jet_feedback.c for the three
 *         injection modes that consume the buffers.
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"

#ifdef BH_JET_LAUNCH_POINT

static int bh_jet_inject_evaluate(int target, int mode, int threadid);

typedef struct
{
  MyDouble Pos[3];
  MyDouble Dir[3];
  MyDouble Energy;
  MyDouble SphereMass;
  MyFloat Hsml;
  integertime MarkTiUntil;
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

static void particle2in(data_in *in, int i, int firstnode)
{
  for(int k = 0; k < 3; k++)
    {
      in->Pos[k] = BhJetPts[i].Pos[k];
      in->Dir[k] = BhJetPts[i].Dir[k];
    }
  in->Energy      = BhJetPts[i].Energy;
  in->SphereMass  = BhJetPts[i].Mass;
  in->Hsml        = BhJetPts[i].Hsml;
  in->MarkTiUntil = BhJetPts[i].MarkTiUntil;
  in->Firstnode   = firstnode;
}

typedef struct
{
} data_out;

static data_out *DataResult, *DataOut;

static void out2particle(data_out *out, int i, int mode) {}

#include "../utils/generic_comm_helpers2.h"

static void kernel_local(void)
{
  int i, j;
  int threadid = get_thread_num();

  for(j = 0; j < NTask; j++)
    Thread[threadid].Exportflag[j] = -1;

  while(1)
    {
      if(Thread[threadid].ExportSpace < MinSpace)
        break;

      i = NextParticle++;

      if(i >= NBhJetPts)
        break;

      bh_jet_inject_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

      bh_jet_inject_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
    }
}

void bh_jet_points_inject(void)
{
  generic_set_MaxNexport();

  generic_comm_pattern(NBhJetPts, kernel_local, kernel_imported);
}

static int bh_jet_inject_evaluate(int target, int mode, int threadid)
{
  int j, n, numnodes, *firstnode;
  MyDouble xtmp, ytmp, ztmp;

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

  MyDouble *pos  = target_data->Pos;
  const double h  = target_data->Hsml;
  const double h2 = h * h;

  const int deposit = target_data->Energy > 0 && target_data->SphereMass > 0;

  /* momentum per unit energy at v_jet, and mass per unit energy (code units) */
  const double v_cgs   = All.BhJetVelocity * 1.0e5;
#ifndef ENERGY_JET
  const double p_per_e = 2.0 * All.cf_UnitEnergy_in_cgs / (v_cgs * All.cf_UnitMomentum_in_cgs);
#endif
#ifdef BH_JET_FIXED_POWER
  const double m_per_e = 2.0 * All.cf_UnitEnergy_in_cgs / (v_cgs * v_cgs * All.cf_UnitMass_in_g);
#endif
  (void)v_cgs;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      if(P[j].Type != 0 || P[j].Mass == 0 || P[j].ID == 0)
        continue;

      double dx = NEAREST_X(P[j].Pos[0] - pos[0]);
      double dy = NEAREST_Y(P[j].Pos[1] - pos[1]);
      double dz = NEAREST_Z(P[j].Pos[2] - pos[2]);

      if(dx * dx + dy * dy + dz * dz >= h2)
        continue;

      if(target_data->MarkTiUntil > SphP[j].BhJetTiUntil)
        SphP[j].BhJetTiUntil = target_data->MarkTiUntil;

      if(!deposit)
        continue;

      double ej = target_data->Energy * P[j].Mass / target_data->SphereMass;

      SphP[j].BhJetEnergy += ej;

      for(int k = 0; k < 3; k++)
#ifdef ENERGY_JET
        SphP[j].BhJetMomentum[k] += ej * target_data->Dir[k]; /* direction accumulator */
#else
        SphP[j].BhJetMomentum[k] += ej * p_per_e * target_data->Dir[k];
#endif
#ifdef BH_JET_FIXED_POWER
      SphP[j].BhJetMass += ej * m_per_e;
#endif
      All.BhFeedbackLocal[0] += ej;
    }

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif /* #ifdef BH_JET_LAUNCH_POINT */
