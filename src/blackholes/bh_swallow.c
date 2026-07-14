#include <stdlib.h>       
#include <math.h>
#include <gsl/gsl_math.h>              
#include <mpi.h>            
  
#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"


static int bh_swallow_evaluate(int target, int mode, int threadid);

static MyFloat *AccretionLimited;

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Pos[3];
  MyDouble NgbsMass;
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
  in->NgbsMass = BhP[i].NgbsMass;
  in->Accretion = BhP[i].Accretion; 
  in->Hsml = BhP[i].Hsml;
  in->Firstnode = firstnode;
}  

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble AccretionLimited;
  MyDouble MassToDrain; 
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
      AccretionLimited[i] = out->AccretionLimited;
      BhP[i].MassToDrain = out->MassToDrain;
    }
  else /* combine */
    {
      AccretionLimited[i] += out->AccretionLimited;
      BhP[i].MassToDrain += out->MassToDrain;
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
        
      bh_swallow_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

      bh_swallow_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
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
void bh_swallow(void)
{
  TIMER_START(CPU_BLACKHOLES_SWALLOW);

  int i, idx;

  AccretionLimited = (MyFloat *)mymalloc("AccretionLimited", NumBhs * sizeof(MyFloat));

  generic_set_MaxNexport();

  generic_comm_pattern(TimeBinsBh.NActiveParticles, kernel_local, kernel_imported);

  for(idx = 0; idx < TimeBinsBh.NActiveParticles; idx++)
    {
      i = TimeBinsBh.ActiveParticleList[idx];
      
      if(AccretionLimited[i])
        BhP[i].Accretion = AccretionLimited[i];
    }

  myfree(AccretionLimited);

  TIMER_STOP(CPU_BLACKHOLES_SWALLOW);
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
static int bh_swallow_evaluate(int target, int mode, int threadid)
{
  int i, n, numnodes, *firstnode; 
  MyDouble xtmp, ytmp, ztmp;   
  MyDouble h, h2, dx, dy, dz, r, r2, wk; 
  MyDouble *pos, ngbsmass, accretion, factor, accretion_limited = 0, mass_to_drain = 0;
  
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
  ngbsmass = target_data->NgbsMass;
  accretion = target_data->Accretion;
  h = target_data->Hsml;

  //MyDouble hinv, hinv3, hinv4, u, dwk;

  //h2   = h * h;
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
          r = sqrt(r2);
          u = r * hinv;

          //bh_kernel(u, hinv3, hinv4, &wk, &dwk);

          factor = P[i].Mass / ngbsmass;

          if(accretion * factor > 0.9 * P[i].Mass)
            {
              SphP[i].BhMassDrain += 0.9 * P[i].Mass;
              
              accretion_limited += 0.9 * P[i].Mass;
              
              mass_to_drain += accretion * factor - 0.9 * P[i].Mass;
            }
          else
            {
              SphP[i].BhMassDrain += accretion * factor;
              
              accretion_limited += accretion * factor;
            }
        } // if(r2 < h2)
    } // for(n = 0; n < nfound; n++)

  out.AccretionLimited = accretion_limited;
  out.MassToDrain = mass_to_drain;

  /* now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}