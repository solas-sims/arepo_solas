#include <stdlib.h>       
#include <math.h>     
#include <mpi.h>            
  
#include "../../main/allvars.h"
#include "../../main/proto.h"

#include "../../domain/domain.h"

static int sf_massdrain_evaluate(int target, int mode, int threadid);

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Pos[3];
  MyDouble MassOfStar;
  MyDouble NgbsMass;
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
    in->Pos[j] = PPS(i).Pos[j];
  in->MassOfStar = SP[i].MassOfStar;
  in->NgbsMass = SP[i].NgbsMass;
  in->Hsml = SP[i].Hsml;
  in->Firstnode = firstnode;
}  

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.
 */
typedef struct
{ 
  MyDouble CM[3];
  MyDouble VM[3];
  MyDouble Metals;
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
      for(int j = 0; j < 3; j++)
        {
          PPS(i).Pos[j] = out->CM[j];
          PPS(i).Vel[j] = out->VM[j];
        }
#ifdef METALS
      SP[i].Metallicity = out->Metals;
#endif
    }
  else /* combine */
    {
      for(int j = 0; j < 3; j++)
        {
          PPS(i).Pos[j] += out->CM[j];
          PPS(i).Vel[j] += out->VM[j];
        }
#ifdef METALS
      SP[i].Metallicity += out->Metals;
#endif
    }
}


#include "../../utils/generic_comm_helpers2.h"

/*! \brief Routine that defines what to do with local particles.
 *
 *  Calls the *_evaluate function in MODE_LOCAL_PARTICLES.
 *
 *  \return void
 */
static void kernel_local(void)
{
  int i, idx, j;

  int threadid = get_thread_num();

  for(j = 0; j < NTask; j++)
    Thread[threadid].Exportflag[j] = -1;

  while(1)
    {
      if(Thread[threadid].ExportSpace < MinSpace)
        break;

      i = NextParticle++;

      if(i >= NumStars)
        break;

      if(PPS(i).Mass == 0)
        sf_massdrain_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

      sf_massdrain_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
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
void sf_massdrain()
{
  generic_set_MaxNexport();
  generic_comm_pattern(NumStars, kernel_local, kernel_imported);
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
static int sf_massdrain_evaluate(int target, int mode, int threadid)
{
  int i, n, numnodes, *firstnode; 
  double h, h2, dx, dy, dz, r, r2, wk; 
  MyDouble *pos, massofstar, ngbsmass, factor;
  MyDouble cm[3], vm[3];

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

  pos = target_data->Pos;
  h = target_data->Hsml;
  h2 = h * h;
  
  massofstar = target_data->MassOfStar;
  ngbsmass = target_data->NgbsMass;

  for(int k = 0; k < 3; k++) 
    cm[k] = vm[k] = 0;

#ifdef METALS
  MyDouble metals = 0;
#endif

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      i = Thread[threadid].Ngblist[n];

      if(P[i].Type != 0 || P[i].Mass == 0 || P[i].ID == 0)
        continue;

      /* compute cell->star position vectors */
      dx = P[i].Pos[0] - pos[0];
      dy = P[i].Pos[1] - pos[1]; 
      dz = P[i].Pos[2] - pos[2]; 

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

          wk = gaussian_weight(r, h);

          factor = P[i].Mass * wk / ngbsmass;

          // compute the mass drain
          SphP[i].StarMassDrain += massofstar * factor;
          // compute center of mass and velocity
          cm[0] += P[i].Pos[0] * factor;
          cm[1] += P[i].Pos[1] * factor;
          cm[2] += P[i].Pos[2] * factor;

          vm[0] += SphP[i].Momentum[0] / P[i].Mass * factor;
          vm[1] += SphP[i].Momentum[1] / P[i].Mass * factor;
          vm[2] += SphP[i].Momentum[2] / P[i].Mass * factor;
#ifdef METALS
          metals += SphP[i].GasMetals * factor;
#endif
        }
    }

  for(int k = 0; k < 3; k++) 
    {
      out.CM[k] = cm[k];
      out.VM[k] = vm[k];
    }
#ifdef METALS
    out.Metals = metals;
#endif

  /* now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}