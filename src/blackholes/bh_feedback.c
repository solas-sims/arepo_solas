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
  MyDouble Pos[3];
  MyDouble NgbsVolume;
  MyDouble Accretion;
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
  for(int j = 0; j < 3; j++)
    in->Pos[j] = PPB(i).Pos[j];
  in->NgbsVolume = BhP[i].NgbsVolume;
  in->Accretion = BhP[i].Accretion;
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
  TIMER_START(CPU_BLACKHOLES_FEEDBACK);

  generic_set_MaxNexport();

  generic_comm_pattern(TimeBinsBh.NActiveParticles, kernel_local, kernel_imported);

  TIMER_STOP(CPU_BLACKHOLES_FEEDBACK);
}

static int bh_feedback_evaluate(int target, int mode, int threadid)
{
  int i, n, numnodes, *firstnode; 
  MyDouble xtmp, ytmp, ztmp;   
  MyDouble h, h2, dx, dy, dz, r2; 
  MyDouble *pos, ngbsvolume, accretion, mass_feed, energy_feed, factor;

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

  pos = target_data->Pos;
  h = target_data->Hsml;
  ngbsvolume = target_data->NgbsVolume;
  accretion = target_data->Accretion;

  mass_feed = All.Mload * All.Epsilon_r * accretion; 
  energy_feed = All.Epsilon_f * (1.0 - All.Mload) * All.Epsilon_r  * accretion * (CLIGHT*CLIGHT / (All.cf_UnitVelocity_in_cm_per_s*All.cf_UnitVelocity_in_cm_per_s));

  //MyDouble hinv, hinv3, hinv4, u, dwk;

  h2 = h * h;
  //hinv = 1.0 / h;
//#ifndef TWODIMS
//  hinv3 = hinv * hinv * hinv;
//#else  /* #ifndef  TWODIMS */
//  hinv3 = hinv * hinv / boxSize_Z;
//#endif /* #ifndef  TWODIMS #else */
//  hinv4 = hinv3 * hinv;

#ifdef BH_CONSTANT_RADIUS
  int nfound = ngb_treefind_variable_threads(pos, All.BhRadius, target, mode, threadid, numnodes, firstnode);
#else
  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);
#endif

  for(n = 0; n < nfound; n++)
    {
      i = Thread[threadid].Ngblist[n];

      if(P[i].Type != 0 || P[i].Mass == 0 || P[i].ID == 0)
        continue;

      /* Compute bh->cell position vector */
      dx = NEAREST_X(P[i].Pos[0] - pos[0]);
      dy = NEAREST_Y(P[i].Pos[1] - pos[1]);
      dz = NEAREST_Z(P[i].Pos[2] - pos[2]);

      r2 = dx * dx + dy * dy + dz * dz;

#ifdef BH_CONSTANT_RADIUS
      if(r2 < All.BhRadius*All.BhRadius)
#else
      if(r2 < h2)
#endif
        {
          //r = sqrt(r2);
          //u = r * hinv;

          //bh_kernel(u, hinv3, hinv4, &wk, &dwk);

          factor = SphP[i].Volume / ngbsvolume;

          /* Add thermal energy isotropically */
          SphP[i].BhMassFeed += mass_feed * factor;

          SphP[i].BhThermalFeed += energy_feed * factor;
          All.BhFeedbackLocal[0] += energy_feed * factor;
        }
    }

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}