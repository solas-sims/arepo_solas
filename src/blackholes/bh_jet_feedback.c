#include <stdlib.h>       
#include <math.h>
#include <gsl/gsl_math.h>              
#include <mpi.h>            
  
#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"

static int bh_feedback_evaluate(int target, int mode, int threadid);

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  int Bin;
  MyDouble Pos[3];
  MyFloat Hsml;
  MyDouble Ngbsmass;
  MyDouble NgbsmassFeed;
  MyDouble Accretion;
#ifdef BH_JET_PRESSURE_TRIGGER
  int JetLaunchFlag;
#endif
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
  for(int j = 0; j < 3; j++)
    in->Pos[j] = PPB(i).Pos[j];
  in->Ngbsmass = BhP[i].Ngbsmass;
  in->NgbsmassFeed = BhP[i].NgbsmassFeed;
  in->Accretion = BhP[i].Accretion;
#ifdef BH_JET_PRESSURE_TRIGGER
  in->JetLaunchFlag = BhP[i].JetLaunchFlag;
#endif
  in->Hsml = BhP[i].Hsml;
  in->Firstnode = firstnode;
}

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.*/

typedef struct
{
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
 *  \return void*/
 
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
    }
  else /* combine */
    {
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
        
      bh_feedback_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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
  int threadid = get_thread_num();

  while(1)
    {
      i = cnt++;

      if(i >= Nimport)
        break;

      bh_feedback_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
    }
}

void bh_feedback(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  generic_set_MaxNexport();

  generic_comm_pattern(TimeBinsBh.NActiveParticles, kernel_local, kernel_imported);

  CPU_Step[CPU_MISC] += measure_time();
}

static int bh_feedback_evaluate(int target, int mode, int threadid)
{
  int i, n, numnodes, *firstnode; 
  double h, h2, r, r2, wk;
  double dx, dy, dz, dvx, dvy, dvz; 
  MyDouble *pos, ngbsmass, accretion, massloading, energyfeed;

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

#ifdef BH_JET_PRESSURE_TRIGGER
  /* jet ram pressure does not exceed the ambient pressure: no launch this step */
  if(!target_data->JetLaunchFlag)
    {
      if(mode == MODE_LOCAL_PARTICLES)
        out2particle(&out, target, MODE_LOCAL_PARTICLES);
      else
        DataResult[target] = out;
      return 0;
    }
#endif

  pos = target_data->Pos;
  h = target_data->Hsml;
  ngbsmass = target_data->Ngbsmass;
  ngbsmassfeed = target_data->NgbsmassFeed;
  accretion = target_data->Accretion;

  h2   = h * h;
  hinv = 1.0 / h;
#ifndef TWODIMS
  hinv3 = hinv * hinv * hinv;
#else  /* #ifndef  TWODIMS */
  hinv3 = hinv * hinv / boxSize_Z;
#endif /* #ifndef  TWODIMS #else */
  hinv4 = hinv3 * hinv;

  massloading = All.Mload * accretion; 
  energyfeed = All.Epsilon_f * All.Epsilon_r * accretion * (CLIGHT * CLIGHT / (All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s));

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

      if(P[j].Type != 0 || P[j].Mass == 0 || P[j].ID == 0)
        continue;

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
          r = sqrt(r2);

          u = r * hinv;

          bh_kernel(u, hinv3, hinv4, &wk, &dwk);

              /* double cone jet setup */

              /* calculate vector to cone vertex */
              vx = -dx; // x-component of the vector from the vertex to the point
              vy = -dy; // y-component of the vector from the vertex to the point
              vz = -dz; // z-component of the vector from the vertex to the point
              /* calculate angles */    
              pos_z_angle = acos((vx*pos_z_axis[0] + 
              vy*pos_z_axis[1] + vz*pos_z_axis[2]) / 
                (sqrt(pow(vx, 2) + pow(vy, 2) + pow(vz, 2)) * sqrt(pow(pos_z_axis[0], 2) + pow(pos_z_axis[1], 2) + pow(pos_z_axis[2], 2))));
              neg_z_angle = acos((vx*neg_z_axis[0] + vy*neg_z_axis[1] + vz*neg_z_axis[2]) / 
                (sqrt(pow(vx, 2) + pow(vy, 2) + pow(vz, 2)) * sqrt(pow(neg_z_axis[0], 2) + pow(neg_z_axis[1], 2) + pow(neg_z_axis[2], 2))));    
              /* check if particle is inside the cone */ 
              if((pos_z_angle <= theta) || (neg_z_angle <= theta))
                {
                  /* add mass */
                  SphP[j].MassLoading += massloading/ngbsmass_feed*P[j].Mass;
                  
                  /* split kinetic and thermal energy feed */ 
                  /* add kinetic energy in cone */
                  SphP[j].KineticFeed   += energyfeed/ngbsmass_feed*P[j].Mass;
                  All.EnergyExchange[0] += energyfeed/ngbsmass_feed*P[j].Mass;

                  /* set radial kick direction */      
                  SphP[j].BhKickVector[0] = vx/r;
                  SphP[j].BhKickVector[1] = vy/r;
                  SphP[j].BhKickVector[2] = vz/r;                
            }
        }
    }

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}