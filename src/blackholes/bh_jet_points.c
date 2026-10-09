/*! \file bh_jet_points.c
 *  \brief Launch-point jet geometry (BH_JET_LAUNCH_POINT): the arepobh model of Maragkakis et al. (2026).
 *
 *  A burst released by a BH (see bh_jet_prepare()) is a packet that flies out along +-d, d a random direction within
 *  BhJetOpeningAngle of the jet axis, and is deposited when it reaches BhJetLaunchRadius. The two lobes sit at
 *  BH.Pos +- R d and receive half of the energy each, spread over the gas inside a sphere around the lobe whose
 *  radius is chosen so it holds (BhJetDesNgb +- BhJetDesDev) * TargetGasMass of gas (as jets_improvements does with
 *  the virtual particle's smoothing length); shares are proportional to cell mass and the kick points along the lobe.
 *
 *  This file finds the sphere radii (bisection like bh_density()) and handles the packet bookkeeping; the cell
 *  buffers are filled by bh_jet_points_inject().
 *
 *  Packets launched on this step are also run through the pass with zero energy, so the cells they will hit are
 *  marked for the jet timestep limiter before the burst arrives. If either lobe's sphere is empty the packet is kept
 *  and retried on the next BH step. Spheres that cannot reach the target mass within the launch radius are accepted
 *  as they are (arepobh stops the run in that case).
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"

#ifdef BH_JET_LAUNCH_POINT

struct bh_jet_point *BhJetPts;
int NBhJetPts;

static int bh_jet_points_density_evaluate(int target, int mode, int threadid);

typedef struct
{
  MyDouble Pos[3];
  MyFloat Hsml;
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

static void particle2in(data_in *in, int i, int firstnode)
{
  for(int k = 0; k < 3; k++)
    in->Pos[k] = BhJetPts[i].Pos[k];
  in->Hsml      = BhJetPts[i].Hsml;
  in->Firstnode = firstnode;
}

typedef struct
{
  MyDouble Mass;
  MyDouble NumNgb;
} data_out;

static data_out *DataResult, *DataOut;

static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      BhJetPts[i].Mass   = out->Mass;
      BhJetPts[i].NumNgb = out->NumNgb;
    }
  else /* combine */
    {
      BhJetPts[i].Mass += out->Mass;
      BhJetPts[i].NumNgb += out->NumNgb;
    }
}

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

      if(BhJetPts[i].Done)
        continue;

      bh_jet_points_density_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

      bh_jet_points_density_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
    }
}

static int bh_jet_points_density_evaluate(int target, int mode, int threadid)
{
  int j, n, numnodes, *firstnode;
  MyDouble xtmp, ytmp, ztmp;
  double mass = 0, ngb = 0;

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

  MyDouble *pos = target_data->Pos;
  double h      = target_data->Hsml;
  double h2     = h * h;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      if(P[j].Type != 0 || P[j].Mass == 0 || P[j].ID == 0)
        continue;

      double dx = NEAREST_X(P[j].Pos[0] - pos[0]);
      double dy = NEAREST_Y(P[j].Pos[1] - pos[1]);
      double dz = NEAREST_Z(P[j].Pos[2] - pos[2]);

      if(dx * dx + dy * dy + dz * dz < h2)
        {
          mass += P[j].Mass;
          ngb += 1;
        }
    }

  out.Mass   = mass;
  out.NumNgb = ngb;

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

void bh_jet_points(void)
{
  TIMER_START(CPU_BLACKHOLES_FEEDBACK);

  MyDouble xtmp, ytmp, ztmp;

  const double r_com = All.BhJetLaunchRadius * 3.0856775814913673e21 / All.cf_UnitLength_in_cm; /* comoving code length */

  /* ---- 1. one pair of points per packet that is new or has arrived ---- */
  int npk = 0;
  for(int idx = 0; idx < TimeBinsBh.NActiveParticles; idx++)
    {
      int i = TimeBinsBh.ActiveParticleList[idx];

      if(i < 0)
        continue;

      for(int k = 0; k < BhP[i].JetNPackets; k++)
        if(BhP[i].JetPacket[k].New || BhP[i].JetPacket[k].TiArrive <= All.Ti_Current)
          npk++;
    }

  BhJetPts  = (struct bh_jet_point *)mymalloc("BhJetPts", (2 * npk > 0 ? 2 * npk : 1) * sizeof(struct bh_jet_point));
  NBhJetPts = 0;

  for(int idx = 0; idx < TimeBinsBh.NActiveParticles; idx++)
    {
      int i = TimeBinsBh.ActiveParticleList[idx];

      if(i < 0)
        continue;

      integertime span = BhP[i].TimeBinBh ? (((integertime)1) << BhP[i].TimeBinBh) : 1;

      for(int k = 0; k < BhP[i].JetNPackets; k++)
        {
          int due = BhP[i].JetPacket[k].TiArrive <= All.Ti_Current;

          if(!due && !BhP[i].JetPacket[k].New)
            continue;

          for(int s = 0; s < 2; s++)
            {
              struct bh_jet_point *pt = &BhJetPts[NBhJetPts++];
              double sign             = s == 0 ? 1.0 : -1.0;

              for(int d = 0; d < 3; d++)
                {
                  pt->Dir[d] = sign * BhP[i].JetPacket[k].Dir[d];
                  pt->Pos[d] = PPB(i).Pos[d] + r_com * pt->Dir[d];
                }
              pt->Pos[0] = WRAP_X(pt->Pos[0]);
              pt->Pos[1] = WRAP_Y(pt->Pos[1]);
              pt->Pos[2] = WRAP_Z(pt->Pos[2]);

              pt->Energy      = due ? 0.5 * BhP[i].JetPacket[k].Energy : 0;
              pt->MarkTiUntil = (due ? All.Ti_Current : BhP[i].JetPacket[k].TiArrive) + 2 * span;
              pt->Hsml        = BhP[i].Hsml > 0 ? BhP[i].Hsml : 0.5 * r_com;
              if(pt->Hsml > r_com)
                pt->Hsml = r_com;
              pt->Mass = pt->NumNgb = pt->Left = pt->Right = 0;
              pt->Bh   = i;
              pt->Pk   = k;
              pt->Due  = due;
              pt->Done = 0;
            }
        }
    }

  /* ---- 2. Hsml of every point so that its sphere holds the target gas mass ---- */
  generic_set_MaxNexport();

  int iter = 0, ncapped = 0;
  long long ntot;

  do
    {
      generic_comm_pattern(NBhJetPts, kernel_local, kernel_imported);

      int nleft = 0;

      for(int p = 0; p < NBhJetPts; p++)
        {
          struct bh_jet_point *pt = &BhJetPts[p];

          if(pt->Done)
            continue;

#ifdef REFINEMENT
          double lo = (All.BhJetDesNgb - All.BhJetDesDev) * All.TargetGasMass;
          double hi = (All.BhJetDesNgb + All.BhJetDesDev) * All.TargetGasMass;
          double v  = pt->Mass;
#else
          double lo = All.BhJetDesNgb - All.BhJetDesDev;
          double hi = All.BhJetDesNgb + All.BhJetDesDev;
          double v  = pt->NumNgb;
#endif
          if(v >= lo && v <= hi)
            {
              pt->Done = 1;
              continue;
            }

          if(pt->Left > 0 && pt->Right > 0 && (pt->Right - pt->Left) < 1.0e-3 * pt->Left)
            {
              pt->Done = 1;
              continue;
            }

          if(v < lo && pt->Hsml >= r_com)
            {
              pt->Done = 1; /* cannot grow further: accept what is there */
              ncapped++;
              continue;
            }

          if(iter > 60)
            {
              pt->Done = 1;
              continue;
            }

          if(v < lo)
            pt->Left = dmax(pt->Hsml, pt->Left);
          else if(pt->Right > 0)
            pt->Right = dmin(pt->Right, pt->Hsml);
          else
            pt->Right = pt->Hsml;

          if(pt->Left > 0 && pt->Right > 0)
            pt->Hsml = pow(0.5 * (pow(pt->Left, 3) + pow(pt->Right, 3)), 1.0 / 3);
          else if(pt->Left > 0)
            pt->Hsml *= 1.26;
          else
            pt->Hsml /= 1.26;

          if(pt->Hsml > r_com)
            pt->Hsml = r_com;

          nleft++;
        }

      sumup_large_ints(1, &nleft, &ntot);
      iter++;
    }
  while(ntot > 0);

  /* ---- 3. a packet is delivered only if both of its lobes found gas ---- */
  int ndeliver = 0, npostponed = 0;
  for(int p = 0; p + 1 < NBhJetPts; p += 2)
    if(BhJetPts[p].Due)
      {
        if(BhJetPts[p].Mass > 0 && BhJetPts[p + 1].Mass > 0)
          ndeliver++;
        else
          {
            BhJetPts[p].Energy = BhJetPts[p + 1].Energy = 0;
            npostponed++;
          }
      }

  bh_jet_points_inject();

  /* ---- 4. drop delivered packets, clear the new-packet flags ---- */
  for(int p = 0; p + 1 < NBhJetPts; p += 2)
    if(BhJetPts[p].Due && BhJetPts[p].Energy > 0)
      BhP[BhJetPts[p].Bh].JetPacket[BhJetPts[p].Pk].Energy = -1;

  long long inflight = 0;
  for(int idx = 0; idx < TimeBinsBh.NActiveParticles; idx++)
    {
      int i = TimeBinsBh.ActiveParticleList[idx];

      if(i < 0)
        continue;

      int m = 0;
      for(int k = 0; k < BhP[i].JetNPackets; k++)
        if(BhP[i].JetPacket[k].Energy >= 0)
          {
            BhP[i].JetPacket[m]     = BhP[i].JetPacket[k];
            BhP[i].JetPacket[m].New = 0;
            m++;
          }
      BhP[i].JetNPackets = m;
      inflight += m;
    }

  long long loc[4] = {ndeliver, npostponed, ncapped, inflight}, glob[4];
  MPI_Allreduce(loc, glob, 4, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

  if(glob[0] + glob[1] + glob[3] > 0)
    mpi_printf("BLACKHOLES: Jet launch points: %lld packets delivered, %lld postponed (empty sphere), %lld spheres capped at the launch radius, "
               "%lld in flight (%d Hsml iterations)\n",
               glob[0], glob[1], glob[2], glob[3], iter);

  myfree(BhJetPts);
  BhJetPts  = NULL;
  NBhJetPts = 0;

  TIMER_STOP(CPU_BLACKHOLES_FEEDBACK);
}

#endif /* #ifdef BH_JET_LAUNCH_POINT */
