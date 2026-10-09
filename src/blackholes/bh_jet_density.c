/*! \file bh_jet_density.c
 *  \brief Mass of the gas cells inside the jet cones of each active black hole.
 *
 *  The jet energy is distributed over the cone cells in proportion to cell mass (m_j / sum m), as in the
 *  Maragkakis et al. (2026) arepobh model, so the total cone mass has to be known before injecting.
 *  The same walk marks the cone cells for the dynamical jet timestep limiter (see get_timestep_hydro()):
 *  a cell is marked while its BH holds jet energy that is about to be, or is waiting to be, released.
 *
 *  Result: BhP[].JetConeMass.
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"

#ifdef BH_JET_FEEDBACK

static int bh_jet_density_evaluate(int target, int mode, int threadid);

typedef struct
{
  MyDouble Pos[3];
  MyDouble Axis[3];
  MyFloat Hsml;
  int Active;            /* 1 if this BH needs the cone mass (it releases energy this step) */
  integertime MarkTiUntil; /* 0: do not mark cone cells for the timestep limiter */
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
  in->Hsml   = BhP[i].Hsml;
  in->Active = BhP[i].JetEnergyRelease > 0;

  in->MarkTiUntil = 0;
  if(BhP[i].JetEnergyRelease > 0 || BhP[i].JetEnergyReservoir > 0)
    {
      integertime span = BhP[i].TimeBinBh ? (((integertime)1) << BhP[i].TimeBinBh) : 1;

      /* the BH renews the mark when it is next active; keep one extra step of margin */
      in->MarkTiUntil = All.Ti_Current + 2 * span;
    }

  in->Firstnode = firstnode;
}

typedef struct
{
  MyDouble ConeMass;
} data_out;

static data_out *DataResult, *DataOut;

static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    BhP[i].JetConeMass = out->ConeMass;
  else /* combine */
    BhP[i].JetConeMass += out->ConeMass;
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

      bh_jet_density_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

      bh_jet_density_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
    }
}

void bh_jet_density(void)
{
  TIMER_START(CPU_BLACKHOLES_FEEDBACK);

  generic_set_MaxNexport();

  generic_comm_pattern(TimeBinsBh.NActiveParticles, kernel_local, kernel_imported);

  TIMER_STOP(CPU_BLACKHOLES_FEEDBACK);
}

static int bh_jet_density_evaluate(int target, int mode, int threadid)
{
  int j, n, numnodes, *firstnode;
  MyDouble xtmp, ytmp, ztmp;
  double d[3];
  double r, conemass = 0;

  data_in local, *target_data;
  data_out out = {0};

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

  if(target_data->Active || target_data->MarkTiUntil)
    {
      const double h        = target_data->Hsml;
      const double h2       = h * h;
      const double cos_th   = cos(All.BhJetOpeningAngle * M_PI / 180.0);
      MyDouble *pos         = target_data->Pos;

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

          conemass += P[j].Mass;

          if(target_data->MarkTiUntil > SphP[j].BhJetTiUntil)
            SphP[j].BhJetTiUntil = target_data->MarkTiUntil;
        }
    }

  out.ConeMass = conemass;

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif /* #ifdef BH_JET_FEEDBACK */
