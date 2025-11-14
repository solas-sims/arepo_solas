#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"

#if defined(BLACKHOLES) && defined(BLACKHOLES_FEEDBACK)

#define DEG_TO_RAD(deg) ((deg) * M_PI / 180.0)

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
  MyDouble Mass;
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
      BhP[i].NgbMassFeed               = out->Mass;
    }
  else /* combine */
    {
      BhNumNgb[i]                      += out->Ngb;
      BhP[i].NgbMassFeed               += out->Mass;
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
void bh_jet_density(void)
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

  generic_set_MaxNexport();

  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
      t0 = second();

      generic_comm_pattern(TimeBinsBh.NActiveParticles, kernel_local, kernel_imported);

      for(idx=0, npleft=0; idx<TimeBinsBh.NActiveParticles; idx++)
        {
          i = TimeBinsBh.ActiveParticleList[idx];

          if(BhP[i].NgbMassFeed < (All.BhDesNgb - All.BhDesDev) || BhP[i].NgbMassFeed > (All.BhDesNgb + All.BhDesDev))
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

            if(BhP[i].NgbMassFeed < (All.BhDesNgb - All.BhDesDev))
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
  double h, h2, hinv, hinv3, hinv4; 
  int j, n, numngb, numnodes, *firstnode;
  double dx, dy, dz, r, r2, u, wk, dwk;
  MyDouble *pos, mass; 
  
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
  
  h2   = h * h;
  hinv = 1.0 / h;
#ifndef TWODIMS
  hinv3 = hinv * hinv * hinv;
#else  /* #ifndef  TWODIMS */
  hinv3 = hinv * hinv / boxSize_Z;
#endif /* #ifndef  TWODIMS #else */
  hinv4 = hinv3 * hinv;

  numngb = mass = 0;

  /* jet axis and opening angle */    

  /* positive and negative jet axes */
  double pos_z_axis[3] = {0, 0, 1};
  double neg_z_axis[3] = {0, 0, -1};      
  /* jet angle */
  double theta = DEG_TO_RAD(10);
  double vx, vy, vz, pos_z_angle, neg_z_angle;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

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

          /* double cone jet setup */    
  
          /* calculate vector to cone vertex */
          vx = -dx; // x-component of the vector from the vertex to the point
          vy = -dy; // y-component of the vector from the vertex to the point
          vz = -dz; // z-component of the vector from the vertex to the point
          /* calculate angles */ 
          pos_z_angle = acos((vx*pos_z_axis[0] + vy*pos_z_axis[1] + vz*pos_z_axis[2]) / 
            (sqrt(pow(vx, 2) + pow(vy, 2) + pow(vz, 2)) * sqrt(pow(pos_z_axis[0], 2) + pow(pos_z_axis[1], 2) + pow(pos_z_axis[2], 2))));
          neg_z_angle = acos((vx*neg_z_axis[0] + vy*neg_z_axis[1] + vz*neg_z_axis[2]) / 
            (sqrt(pow(vx, 2) + pow(vy, 2) + pow(vz, 2)) * sqrt(pow(neg_z_axis[0], 2) + pow(neg_z_axis[1], 2) + pow(neg_z_axis[2], 2))));   
          /* check if particle is inside the cone */ 
          if((pos_z_angle <= theta) || (neg_z_angle <= theta))              
            {  
              /* compute ngb mass */
              mass += P[j].Mass;
            }

        }
    }
  
  out.Ngb                     = numngb;
  out.Mass                    = mass;

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

#endif /* #if defined(BLACKHOLES) && defined(BLACKHOLES_FEEDBACK) */