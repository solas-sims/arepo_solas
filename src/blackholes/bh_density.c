#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"

#ifdef BLACKHOLES

static int bh_density_evaluate(int target, int mode, int threadid);
static int bh_density_isactive(int n);

static MyFloat *BhNumNgb;

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Pos[3];
#ifdef BONDI_ACCRETION
  MyDouble Vel[3];
#endif
#ifdef TORQUE_ACCRETION
  MyDouble Vel[3];
#endif
#ifdef ADP_ACCRETION
  MyDouble Vel[3];
#endif
  MyFloat Hsml;
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

/*! \brief Routine that fills the relevant particle/cell data into the input
 *         structure defined above. Needed by generic_comm_helpers2.
 *
 *  \param[out] in Data structure to fill.
 *  \param[in] i Index of particle in P and SphP arrays.
 *  \param[in] firstnode First note of communication.
 *
 *  \return void
 */
static void particle2in(data_in *in, int i, int firstnode)
{
  in->Pos[0]        = PPB(i).Pos[0];
  in->Pos[1]        = PPB(i).Pos[1];
  in->Pos[2]        = PPB(i).Pos[2];
#ifdef BONDI_ACCRETION
  in->Vel[0]        = PPB(i).Vel[0];
  in->Vel[1]        = PPB(i).Vel[1];
  in->Vel[2]        = PPB(i).Vel[2];
#endif
#ifdef TORQUE_ACCRETION
  in->Vel[0]        = PPB(i).Vel[0];
  in->Vel[1]        = PPB(i).Vel[1];
  in->Vel[2]        = PPB(i).Vel[2];
#endif
#ifdef ADP_ACCRETION
  in->Vel[0]        = PPB(i).Vel[0];
  in->Vel[1]        = PPB(i).Vel[1];
  in->Vel[2]        = PPB(i).Vel[2];
#endif
  in->Hsml          = BhP[i].Hsml;
  in->Firstnode     = firstnode;
}  

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Ngb;
  MyDouble Rho;
  MyDouble Mass;
  integertime NgbMinStep;
#ifdef BONDI_ACCRETION
  MyDouble VelocityGas[3];
  MyDouble VelocityGasCircular[3];
  MyDouble InternalEnergyGas;
#endif
#ifdef INFALL_ACCRETION
  MyDouble Accretion;
#endif
#ifdef TORQUE_ACCRETION
  MyDouble TorqueMgas;
  MyDouble TorqueMstar;
  MyDouble TorqueMgasDisk;      /* Disk component gas mass */
  MyDouble TorqueMstarDisk;     /* Disk component stellar mass */ 
  MyDouble TorqueR0;
  MyDouble TorqueFd;
  MyDouble VelocityGasCircular[3];
#endif
#ifdef ADP_ACCRETION
  MyDouble ADP_Racc;
  MyDouble ADP_CapturedMass;
  MyDouble ADP_ReservoirMass;
  MyDouble ADP_DiscMass;
  MyDouble VelocityGasCircular[3];
  MyDouble VelocityGas[3];
#endif
} data_out;

static data_out *DataResult, *DataOut;

/*! \brief Routine to store or combine result data. Needed by
 *         generic_comm_helpers2.
 *
 *  \param[in] out Data to be moved to appropriate variables in global
 *  particle and cell data arrays (P, SphP,...)
 *  \param[in] i Index of particle in P and SphP arrays
 *  \param[in] mode Mode of function: local particles or information that was
 *  communicated from other tasks and has to be added locally?
 *
 *  \return void
 */
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      BhNumNgb[i]                      = out->Ngb;
      BhP[i].Density                   = out->Rho;
      BhP[i].NgbMass                   = out->Mass;
      BhP[i].NgbMinStep                = out->NgbMinStep;
#ifdef BONDI_ACCRETION
      BhP[i].VelocityGas[0]            = out->VelocityGas[0];
      BhP[i].VelocityGas[1]            = out->VelocityGas[1];
      BhP[i].VelocityGas[2]            = out->VelocityGas[2];
      BhP[i].VelocityGasCircular[0]    = out->VelocityGasCircular[0];
      BhP[i].VelocityGasCircular[1]    = out->VelocityGasCircular[1];
      BhP[i].VelocityGasCircular[2]    = out->VelocityGasCircular[2];
      BhP[i].InternalEnergyGas         = out->InternalEnergyGas;
#endif
#ifdef INFALL_ACCRETION
      BhP[i].Accretion                += out->Accretion;
#endif
#ifdef TORQUE_ACCRETION
      BhP[i].TorqueMgas                = out->TorqueMgas;
      BhP[i].TorqueMstar               = out->TorqueMstar;
      BhP[i].TorqueMgasDisk            = out->TorqueMgasDisk;
      BhP[i].TorqueMstarDisk           = out->TorqueMstarDisk;
      BhP[i].TorqueR0                  = out->TorqueR0;
      BhP[i].TorqueFd                  = out->TorqueFd;
      BhP[i].VelocityGasCircular[0]    = out->VelocityGasCircular[0];
      BhP[i].VelocityGasCircular[1]    = out->VelocityGasCircular[1];
      BhP[i].VelocityGasCircular[2]    = out->VelocityGasCircular[2];
#endif
#ifdef ADP_ACCRETION
     BhP[i].ADP_Racc                   = out->ADP_Racc;
     BhP[i].ADP_CapturedMass           = out->ADP_CapturedMass;
     BhP[i].ADP_ReservoirMass          = out->ADP_ReservoirMass;
     BhP[i].ADP_DiscMass               = out->ADP_DiscMass;
     BhP[i].VelocityGasCircular[0]     = out->VelocityGasCircular[0];
     BhP[i].VelocityGasCircular[1]     = out->VelocityGasCircular[1];
     BhP[i].VelocityGasCircular[2]     = out->VelocityGasCircular[2];
     BhP[i].VelocityGas[0]             = out->VelocityGas[0];
     BhP[i].VelocityGas[1]             = out->VelocityGas[1];
     BhP[i].VelocityGas[2]             = out->VelocityGas[2];
#endif
    }
  else /* combine */
    {
      BhNumNgb[i]                      += out->Ngb;
      BhP[i].Density                   += out->Rho;
      BhP[i].NgbMass                   += out->Mass;
      if(out->NgbMinStep < BhP[i].NgbMinStep)
        BhP[i].NgbMinStep               = out->NgbMinStep;
#ifdef BONDI_ACCRETION
      BhP[i].VelocityGas[0]            += out->VelocityGas[0];
      BhP[i].VelocityGas[1]            += out->VelocityGas[1];
      BhP[i].VelocityGas[2]            += out->VelocityGas[2];
      BhP[i].VelocityGasCircular[0]    += out->VelocityGasCircular[0];
      BhP[i].VelocityGasCircular[1]    += out->VelocityGasCircular[1];
      BhP[i].VelocityGasCircular[2]    += out->VelocityGasCircular[2];
      BhP[i].InternalEnergyGas         += out->InternalEnergyGas;
#endif
#ifdef INFALL_ACCRETION
      BhP[i].Accretion                 += out->Accretion; 
#endif
#ifdef TORQUE_ACCRETION
      BhP[i].TorqueMgas                += out->TorqueMgas;
      BhP[i].TorqueMstar               += out->TorqueMstar;
      BhP[i].TorqueMgasDisk            += out->TorqueMgasDisk;
      BhP[i].TorqueMstarDisk           += out->TorqueMstarDisk;
      if(out->TorqueR0 > BhP[i].TorqueR0)
        BhP[i].TorqueR0                 = out->TorqueR0;
      BhP[i].TorqueFd                  += out->TorqueFd;
      BhP[i].VelocityGasCircular[0]    += out->VelocityGasCircular[0];
      BhP[i].VelocityGasCircular[1]    += out->VelocityGasCircular[1];
      BhP[i].VelocityGasCircular[2]    += out->VelocityGasCircular[2];
#endif
#ifdef ADP_ACCRETION
     BhP[i].ADP_CapturedMass           += out->ADP_CapturedMass;
     BhP[i].ADP_ReservoirMass          += out->ADP_ReservoirMass;
     BhP[i].ADP_DiscMass               += out->ADP_DiscMass;
     BhP[i].VelocityGasCircular[0]     += out->VelocityGasCircular[0];
     BhP[i].VelocityGasCircular[1]     += out->VelocityGasCircular[1];
     BhP[i].VelocityGasCircular[2]     += out->VelocityGasCircular[2];
     BhP[i].VelocityGas[0]             += out->VelocityGas[0];
     BhP[i].VelocityGas[1]             += out->VelocityGas[1];
     BhP[i].VelocityGas[2]             += out->VelocityGas[2];
#endif
    }
}


#include "../utils/generic_comm_helpers2.h"

/*! \brief Routine that defines what to do with local particles.
 *
 *  Calls the *_evaluate function in MODE_LOCAL_PARTICLES.
 *
 *  \return void
 */
static void kernel_local(void)
{
  int i, idx;

  {
    int j, threadid = get_thread_num();

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

        //i = NextParticle++;

        //if(i >= NumBhs)
        //break;
        
        idx = NextParticle++;

        if(idx >= TimeBinsBh.NActiveParticles)
          break;

        i = TimeBinsBh.ActiveParticleList[idx];

        if(bh_density_isactive(i))
          bh_density_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

/*! \brief Routine that defines what to do with imported particles.
 *
 *  Calls the *_evaluate function in MODE_IMPORTED_PARTICLES.
 *
 *  \return void
 */
static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
  {
    int threadid = get_thread_num();

    while(1)
      {
        i = cnt++;

        if(i >= Nimport)
          break;

        bh_density_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

/*! \brief Main function of SPH density calculation.
 *
 *  This function computes the local density for each active SPH particle and
 *  the number of weighted neighbors in the current smoothing radius. If a
 *  particle with its smoothing region is fully inside the local domain, it is
 *  not exported to the other processors. The function also detects particles
 *  that have a number of neighbors outside the allowed tolerance range. For
 *  these particles, the smoothing length is adjusted accordingly, and the
 *  computation is called again.
 *
 *  \return void
 */
void bh_density(void)
{
  MyFloat *Left, *Right;
  int idx, i, npleft, iter = 0;
  long long ntot;
  double t0, t1;

  CPU_Step[CPU_MISC] += measure_time();

  BhNumNgb  = (MyFloat *)mymalloc("BhNumNgb", NumBhs * sizeof(MyFloat));
  Left      = (MyFloat *)mymalloc("Left", NumBhs * sizeof(MyFloat));
  Right     = (MyFloat *)mymalloc("Right", NumBhs * sizeof(MyFloat));

  for(i = 0; i < NumBhs; i++)
    {
      Left[i] = Right[i] = 0;
      BhP[i].DensityFlag = 1;
    }

  mpi_printf("BH_DENSITY: Start density and neighbour search for %d black holes.\n", NumBhs);

  generic_set_MaxNexport();

    for(idx=0; idx<TimeBinsBh.NActiveParticles; idx++)
    {
      i = TimeBinsBh.ActiveParticleList[idx];
      if(BhP[i].Hsml <= 0)
      {
        mpi_printf("WARNING: BH %d has invalid Hsml=%g, reinitializing\n",
                    i, BhP[i].Hsml);
        // Use softening as fallback
        BhP[i].Hsml = All.SofteningTable[P[BhP[i].PID].SofteningType];
      }
    }
 
  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
      t0 = second();

      generic_comm_pattern(TimeBinsBh.NActiveParticles, kernel_local, kernel_imported);

      for(idx=0, npleft=0; idx<TimeBinsBh.NActiveParticles; idx++)
        {
          i = TimeBinsBh.ActiveParticleList[idx];

          if(BhNumNgb[i] < (All.BhDesNgb - All.BhDesDev) || BhNumNgb[i] > (All.BhDesNgb + All.BhDesDev))
          {
                  /* need to redo this particle */
            npleft++;

            if(Left[i] > 0 && Right[i] > 0)
              {
                if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
                  {
                        /* this one should be ok */
                    npleft--;
                    BhP[i].DensityFlag = -1; /* Mark as inactive */
                    continue;
                }
              } 

            if(BhNumNgb[i] < (All.BhDesNgb - All.BhDesDev))
              Left[i] = dmax(BhP[i].Hsml, Left[i]);
            else
              {
                if(Right[i] != 0)
                  {
                    if(BhP[i].Hsml < Right[i])
                        Right[i] = BhP[i].Hsml;
                  }
                    else
                        Right[i] = BhP[i].Hsml;
              }

            if(Right[i] > 0 && Left[i] > 0)
                BhP[i].Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
            else
              {
                if(Right[i] == 0 && Left[i] == 0)
                    terminate("should not occur");

                if(Right[i] == 0 && Left[i] > 0)
                  {
                    BhP[i].Hsml *= 1.26;
                  }

                if(Right[i] > 0 && Left[i] == 0)
                  {
                    BhP[i].Hsml /= 1.26;
                  }
              }
          }
        else
             BhP[i].DensityFlag = -1; /* Mark as inactive */ 
        }
        
      sumup_large_ints(1, &npleft, &ntot);

      t1 = second();

      if(ntot > 0)
        {
          iter++;

          if(iter > 0)
            mpi_printf("BH_DENSITY: ngb iteration %3d: need to repeat for %12lld particles. (took %g sec)\n", iter, ntot,
                       timediff(t0, t1));

          if(iter > MAXITER)
            terminate("failed to converge in neighbour iteration in bh_density()\n");
        }
    }
  while(ntot > 0);

  myfree(Right);
  myfree(Left);
  myfree(BhNumNgb);

  /* mark as active again */
  for(i = 0; i < NumBhs; i++)
    {
     BhP[i].DensityFlag = 1;
    }
  
  /* collect some timing information */
  CPU_Step[CPU_INIT] += measure_time();
}

/*! \brief Inner function of the SPH density calculation
 *
 *  This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 *
 *  \param[in] target Index of particle in local data/import buffer.
 *  \param[in] mode Mode in which function is called (local or impored data).
 *  \param[in] threadid ID of local thread.
 *
 *  \return 0
 */
static int bh_density_evaluate(int target, int mode, int threadid)
{
  int bin = TIMEBINS;
  double h, h2, hinv, hinv3, hinv4; 
  int j, n, numngb, numnodes, *firstnode;
  double dx, dy, dz, r, r2, u, wk, dwk;
  MyDouble *pos, mass_j, rho, mass; 
  integertime ngb_min_step;
  
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

  pos  = target_data->Pos;
  h    = target_data->Hsml;

#ifdef BONDI_ACCRETION
  MyDouble *vel;
  vel  = target_data->Vel;
  double dvx, dvy, dvz, rho_j;
  MyDouble internal_energy_gas = 0;
  MyDouble velocity_gas[3], velocity_gas_circular[3];
  velocity_gas[0] = velocity_gas[1] = velocity_gas[2] = 0;
  velocity_gas_circular[0] = velocity_gas_circular[1] = velocity_gas_circular[2] = 0;
#endif
#ifdef INFALL_ACCRETION
  MyDouble accretion = 0;
  double rbh  = h;
  double rbh2 = rbh * rbh;
#endif 
#ifdef TORQUE_ACCRETION
  /* Torque accretion model from Angles-Alcazar et al. 2016 */
  MyDouble torque_Mgas       = 0.0;  /* Total gas mass within R0 */
  MyDouble torque_Mstar      = 0.0;  /* Total stellar mass within R0 */
  MyDouble torque_Mgas_disk  = 0.0;  /* Disk component gas mass */
  MyDouble torque_Mstar_disk = 0.0;  /* Disk component stellar mass */
  MyDouble R0_torque = All.TorqueR0;   // 0.2-0.3 kpc in code units
  MyDouble R0_torque2 = R0_torque * R0_torque;
  /* Angular momentum for disk decomposition */
  MyDouble ang_mom[3]        = {0.0, 0.0, 0.0};
  MyDouble ang_mom_norm      = 0.0;
  MyDouble total_mass_for_angmom = 0.0;
  
  /* For angular momentum accretion tracking */
  MyDouble velocity_gas_circular[3] = {0.0, 0.0, 0.0};
  MyDouble *vel = target_data->Vel;
  double dvx, dvy, dvz;
#endif
#ifdef ADP_ACCRETION
  MyDouble adp_captured_mass = 0.0;
  MyDouble Racc  = h * 24 ;
  MyDouble Racc2 = Racc * Racc;
  MyDouble total_mass_for_angmom = 0.0;
  MyDouble ang_mom[3] = {0.0, 0.0, 0.0};
  MyDouble ang_mom_norm = 0.0;
  MyDouble *vel = target_data->Vel;
  MyDouble velocity_gas_circular[3] = {0.0, 0.0, 0.0};
  MyDouble velocity_gas[3] = {0.0, 0.0, 0.0};
  double dvx, dvy, dvz,dv2;
#endif
  h2   = h * h;
  hinv = 1.0 / h;
#ifndef TWODIMS
  hinv3 = hinv * hinv * hinv;
#else  /* #ifndef  TWODIMS */
  hinv3 = hinv * hinv / boxSize_Z;
#endif /* #ifndef  TWODIMS #else */
  hinv4 = hinv3 * hinv;

  numngb = rho = mass = 0;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

#ifdef ADP_ACCRETION
for(n = 0; n < nfound; n++)
  {
    j = Thread[threadid].Ngblist[n];

    dx = pos[0] - P[j].Pos[0];
    dy = pos[1] - P[j].Pos[1];
    dz = pos[2] - P[j].Pos[2];

    
    #ifndef REFLECTIVE_X
      if(dx >  boxHalf_X) dx -= boxSize_X;
      if(dx < -boxHalf_X) dx += boxSize_X;
    #endif
    #ifndef REFLECTIVE_Y
      if(dy >  boxHalf_Y) dy -= boxSize_Y;
      if(dy < -boxHalf_Y) dy += boxSize_Y;
    #endif
    #ifndef REFLECTIVE_Z
      if(dz >  boxHalf_Z) dz -= boxSize_Z;
      if(dz < -boxHalf_Z) dz += boxSize_Z;
    #endif

    r2 = dx*dx + dy*dy + dz*dz;

    if(r2 < Racc2 && P[j].Type == 0)   /* only gas within Racc */
      {
        mass_j = P[j].Mass;

        dvx = P[j].Vel[0] - vel[0];
        dvy = P[j].Vel[1] - vel[1];
        dvz = P[j].Vel[2] - vel[2];

        /* accumulate angular momentum: L = r × mv */
        ang_mom[0] += mass_j * (dy*dvz - dz*dvy);
        ang_mom[1] += mass_j * (dz*dvx - dx*dvz);
        ang_mom[2] += mass_j * (dx*dvy - dy*dvx);
        total_mass_for_angmom += mass_j;

        
        adp_captured_mass += mass_j;
      }
  }


ang_mom_norm = sqrt(ang_mom[0]*ang_mom[0] +
                    ang_mom[1]*ang_mom[1] +
                    ang_mom[2]*ang_mom[2]);
if(ang_mom_norm > 0)
  {
    ang_mom[0] /= ang_mom_norm;
    ang_mom[1] /= ang_mom_norm;
    ang_mom[2] /= ang_mom_norm;
  }
#endif
#ifdef ADP_ACCRETION
/* ── ADP PASS 2: Disk classification + circular velocity ── */
for(n = 0; n < nfound; n++)
  {
    j = Thread[threadid].Ngblist[n];

    dx = pos[0] - P[j].Pos[0];
    dy = pos[1] - P[j].Pos[1];
    dz = pos[2] - P[j].Pos[2];

    /* periodic wrapping */
    #ifndef REFLECTIVE_X
      if(dx >  boxHalf_X) dx -= boxSize_X;
      if(dx < -boxHalf_X) dx += boxSize_X;
    #endif
    #ifndef REFLECTIVE_Y
      if(dy >  boxHalf_Y) dy -= boxSize_Y;
      if(dy < -boxHalf_Y) dy += boxSize_Y;
    #endif
    #ifndef REFLECTIVE_Z
      if(dz >  boxHalf_Z) dz -= boxSize_Z;
      if(dz < -boxHalf_Z) dz += boxSize_Z;
    #endif

    r2 = dx*dx + dy*dy + dz*dz;

    if(r2 < Racc2 && P[j].Type == 0)   /* only gas within Racc */
      {
        mass_j = P[j].Mass;
        r = sqrt(r2);

        dvx = P[j].Vel[0] - vel[0];
        dvy = P[j].Vel[1] - vel[1];
        dvz = P[j].Vel[2] - vel[2];

        double v_cross[3];
        v_cross[0] = dy*dvz - dz*dvy;
        v_cross[1] = dz*dvx - dx*dvz;
        v_cross[2] = dx*dvy - dy*dvx;

        double v_phi = 0.0;
        if(ang_mom_norm > 0 && r > 0)
          v_phi = (v_cross[0]*ang_mom[0] +
                   v_cross[1]*ang_mom[1] +
                   v_cross[2]*ang_mom[2]) / r;

        int is_disk = (v_phi > 0.0);

        double rho_j = (SphP[j].Density > 0) ? SphP[j].Density : 1.0;
        velocity_gas_circular[0] += v_cross[0] * mass_j / rho_j;
        velocity_gas_circular[1] += v_cross[1] * mass_j / rho_j;
        velocity_gas_circular[2] += v_cross[2] * mass_j / rho_j;
      }
  }
#endif
/* ── PASS 1: Angular momentum axis ── */
#ifdef TORQUE_ACCRETION
  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      dx = pos[0] - P[j].Pos[0];
      dy = pos[1] - P[j].Pos[1];
      dz = pos[2] - P[j].Pos[2];

      /* periodic wrapping */
      #ifndef REFLECTIVE_X
        if(dx >  boxHalf_X) dx -= boxSize_X;
        if(dx < -boxHalf_X) dx += boxSize_X;
      #endif
      #ifndef REFLECTIVE_Y
        if(dy >  boxHalf_Y) dy -= boxSize_Y;
        if(dy < -boxHalf_Y) dy += boxSize_Y;
      #endif
      #ifndef REFLECTIVE_Z
        if(dz >  boxHalf_Z) dz -= boxSize_Z;
        if(dz < -boxHalf_Z) dz += boxSize_Z;
      #endif

      r2 = dx*dx + dy*dy + dz*dz;

      if(r2 < R0_torque2)
        {
          mass_j = P[j].Mass;
          dvx = P[j].Vel[0] - vel[0];
          dvy = P[j].Vel[1] - vel[1];
          dvz = P[j].Vel[2] - vel[2];

          ang_mom[0] += mass_j * (dy*dvz - dz*dvy);
          ang_mom[1] += mass_j * (dz*dvx - dx*dvz);
          ang_mom[2] += mass_j * (dx*dvy - dy*dvx);
          total_mass_for_angmom += mass_j;
        }
    }

  /* Normalise BEFORE Pass 2 */
  ang_mom_norm = sqrt(ang_mom[0]*ang_mom[0] +
                      ang_mom[1]*ang_mom[1] +
                      ang_mom[2]*ang_mom[2]);
  if(ang_mom_norm > 0)
    {
      ang_mom[0] /= ang_mom_norm;
      ang_mom[1] /= ang_mom_norm;
      ang_mom[2] /= ang_mom_norm;
    }
#endif


/* ── PASS 2: Disk decomposition ── */
#ifdef TORQUE_ACCRETION
  for(n = 0; n < nfound; n++)        
    {
      j = Thread[threadid].Ngblist[n];

      dx = pos[0] - P[j].Pos[0];
      dy = pos[1] - P[j].Pos[1];
      dz = pos[2] - P[j].Pos[2];

      /* periodic wrapping */
      #ifndef REFLECTIVE_X
        if(dx >  boxHalf_X) dx -= boxSize_X;
        if(dx < -boxHalf_X) dx += boxSize_X;
      #endif
      #ifndef REFLECTIVE_Y
        if(dy >  boxHalf_Y) dy -= boxSize_Y;
        if(dy < -boxHalf_Y) dy += boxSize_Y;
      #endif
      #ifndef REFLECTIVE_Z
        if(dz >  boxHalf_Z) dz -= boxSize_Z;
        if(dz < -boxHalf_Z) dz += boxSize_Z;
      #endif

      r2 = dx*dx + dy*dy + dz*dz;

      if(r2 < R0_torque2)            
        {
          mass_j = P[j].Mass;
          r = sqrt(r2);

          dvx = P[j].Vel[0] - vel[0];
          dvy = P[j].Vel[1] - vel[1];
          dvz = P[j].Vel[2] - vel[2];

          double v_cross[3];
          v_cross[0] = dy*dvz - dz*dvy;
          v_cross[1] = dz*dvx - dx*dvz;
          v_cross[2] = dx*dvy - dy*dvx;

          double v_phi = 0.0;
          if(ang_mom_norm > 0 && r > 0)   
            v_phi = (v_cross[0]*ang_mom[0] +
                     v_cross[1]*ang_mom[1] +
                     v_cross[2]*ang_mom[2]) / r;

          int is_disk = (v_phi > 0.0);

          if(P[j].Type == 0)   /* Gas */
            {
              torque_Mgas += mass_j;
              if(is_disk) torque_Mgas_disk += mass_j;

              double rho_j = (SphP[j].Density > 0) ? SphP[j].Density : 1.0;
              velocity_gas_circular[0] += v_cross[0] * mass_j / rho_j;
              velocity_gas_circular[1] += v_cross[1] * mass_j / rho_j;
              velocity_gas_circular[2] += v_cross[2] * mass_j / rho_j;
            }
          else if(P[j].Type == 4)   /* Stars */
            {
              torque_Mstar += mass_j;
              if(is_disk) torque_Mstar_disk += mass_j;
            }
        }
    }
#endif
#ifdef TORQUE_ACCRETION
  MyDouble M_disk  = torque_Mgas_disk + torque_Mstar_disk;
  MyDouble M_total = torque_Mgas + torque_Mstar;
  MyDouble f_d     = (M_total > 0) ? M_disk / M_total : 0.0;
#endif

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      dx = pos[0] - P[j].Pos[0];
      dy = pos[1] - P[j].Pos[1];
      dz = pos[2] - P[j].Pos[2];

#ifndef REFLECTIVE_X
      if(dx > boxHalf_X)
        dx -= boxSize_X;
      if(dx < -boxHalf_X)
        dx += boxSize_X;
#endif /* #ifndef REFLECTIVE_X */

#ifndef REFLECTIVE_Y
      if(dy > boxHalf_Y)
        dy -= boxSize_Y;
      if(dy < -boxHalf_Y)
        dy += boxSize_Y;
#endif /* #ifndef REFLECTIVE_Y */

#ifndef REFLECTIVE_Z
      if(dz > boxHalf_Z)
        dz -= boxSize_Z;
      if(dz < -boxHalf_Z)
        dz += boxSize_Z;
#endif /* #ifndef REFLECTIVE_Z */
      r2 = dx * dx + dy * dy + dz * dz;

      if(r2 < h2)
        {
          numngb++;

          r = sqrt(r2);

          u = r * hinv;

          kernel(u, hinv3, hinv4, &wk, &dwk);

          mass_j = P[j].Mass;

          /* compute bh density */
          rho +=  mass_j * wk;

          /* compute the bh-ngb-mass (sphere) */
          mass += mass_j;

          /* compute the min hydro step for neighbors */     
          if(bin > P[j].TimeBinHydro)
            bin = P[j].TimeBinHydro;

#ifdef BONDI_ACCRETION
          /* comute relative velocities, 
               relative specific angular momenta and internal energy of gas */
          dvx = P[j].Vel[0] - vel[0]; 
          dvy = P[j].Vel[1] - vel[1]; 
          dvz = P[j].Vel[2] - vel[2]; 

          if(SphP[j].Density > 0)
            rho_j  = SphP[j].Density;
          else
            rho_j = 1;

          velocity_gas[0] += dvx*mass_j/rho_j*wk;
          velocity_gas[1] += dvy*mass_j/rho_j*wk;
          velocity_gas[2] += dvz*mass_j/rho_j*wk;

          velocity_gas_circular[0] -= (dy * dvz - dz * dvy)*mass_j/rho_j*wk;
          velocity_gas_circular[1] -= (dz * dvx - dx * dvz)*mass_j/rho_j*wk;
          velocity_gas_circular[2] -= (dx * dvy - dy * dvx)*mass_j/rho_j*wk;

          internal_energy_gas += SphP[j].Utherm*mass_j/rho_j*wk;
#endif
#ifdef INFALL_ACCRETION
          /* cell nibbled */
          if(r < 2*rbh) 
            {
              accretion += P[j].Mass * exp(-r2/(2*rbh2));
              P[j].Mass -= P[j].Mass * exp(-r2/(2*rbh2));  
            }
#endif
        } // if(r2 < h2)
    } // for(n = 0; n < nfound; n++)

  
  /* compute bh timestep based on min ngb timestep */
  if(bin == 0)
    ngb_min_step = 0;
  else
    ngb_min_step   = (((integertime)1) << bin);
  
  out.Ngb                     = numngb;
  out.Rho                     = rho;
  out.Mass                    = mass;
  out.NgbMinStep              = ngb_min_step;
#ifdef BONDI_ACCRETION
  out.VelocityGas[0]          = velocity_gas[0];
  out.VelocityGas[1]          = velocity_gas[1];
  out.VelocityGas[2]          = velocity_gas[2];
  out.VelocityGasCircular[0]  = velocity_gas_circular[0];
  out.VelocityGasCircular[1]  = velocity_gas_circular[1];
  out.VelocityGasCircular[2]  = velocity_gas_circular[2];
  out.InternalEnergyGas       = internal_energy_gas;
#endif
#ifdef INFALL_ACCRETION
  out.Accretion               = accretion;
#endif
#ifdef TORQUE_ACCRETION
  out.TorqueMgas              = torque_Mgas;
  out.TorqueMstar             = torque_Mstar;
  out.TorqueMgasDisk          = torque_Mgas_disk;
  out.TorqueMstarDisk         = torque_Mstar_disk;
  out.TorqueR0                = h;  
  out.TorqueFd                = f_d;
  out.VelocityGasCircular[0]  = velocity_gas_circular[0];
  out.VelocityGasCircular[1]  = velocity_gas_circular[1];
  out.VelocityGasCircular[2]  = velocity_gas_circular[2];
#endif
#ifdef ADP_ACCRETION
  out.ADP_Racc                = Racc;
  out.ADP_CapturedMass        = adp_captured_mass;
  out.VelocityGasCircular[0]  = velocity_gas_circular[0];
  out.VelocityGasCircular[1]  = velocity_gas_circular[1];
  out.VelocityGasCircular[2]  = velocity_gas_circular[2];
  out.VelocityGas[0]          = velocity_gas[0];
  out.VelocityGas[1]          = velocity_gas[1];
  out.VelocityGas[2]          = velocity_gas[2];
#endif

  /* now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

/* \brief Determines if a BhP is active in current timestep.
 *
 *  \param[in] n Index of BhP in Particle array
 *
 *  \return 1: BhP active; 0: BhP not active.
 */
int bh_density_isactive(int n)
{
  if(BhP[n].DensityFlag < 0)
    return 0;

  return 1;
}

#endif /* #ifdef BLACKHOLES */
