#include <stdlib.h>       
#include <math.h>         
#include <mpi.h>            
  
#include "../../main/allvars.h"
#include "../../main/proto.h"

#include "../../domain/domain.h"

static int sf_starbystar_evaluate(int target, int mode, int threadid);
static int sf_starbystar_isactive(int n);

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
 typedef struct
{
  MyDouble Pos[3];
  MyDouble MassOfStar;
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
  in->Hsml = SP[i].Hsml;
  in->Firstnode = firstnode;
}  

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.
 */
typedef struct
{ 
  MyDouble NgbsMass;
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
      SP[i].NgbsMass = out->NgbsMass;
    }
  else /* combine */
    {
      SP[i].NgbsMass += out->NgbsMass;
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

      if((PPS(i).Mass == 0) && (sf_starbystar_isactive(i)))
        sf_starbystar_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

      sf_starbystar_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
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
void sf_starbystar()
{
  MyFloat *Left, *Right;
  int idx, i, npleft, iter = 0;
  long long ntot;
  double t0, t1;

  CPU_Step[CPU_MISC] += measure_time();

  Left = (MyFloat *)mymalloc("Left", NumStars * sizeof(MyFloat));
  Right = (MyFloat *)mymalloc("Right", NumStars * sizeof(MyFloat));

  for(i = 0; i < NumStars; i++)
    {
      Left[i] = Right[i] = 0;
      SP[i].DensityFlag = 1;
      
      if(SP[i].Hsml <= 0)
        SP[i].Hsml = All.SofteningTable[PPS(i).SofteningType];
    }

  generic_set_MaxNexport();

  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
      t0 = second();

      generic_comm_pattern(NumStars, kernel_local, kernel_imported);

      for(i = 0, npleft = 0; i < NumStars; i++)
        {
          if(PPS(i).Mass != 0)
            continue;
          
          if(SP[i].Hsml > 10 * All.SofteningTable[PPS(i).SofteningType])
            terminate("Star formation radius too large!");

          if(SP[i].NgbsMass < (5 * SP[i].MassOfStar) || SP[i].NgbsMass > (10 * SP[i].MassOfStar))
            {
              /* need to redo this particle */
              npleft++;

              if(Left[i] > 0 && Right[i] > 0)
                {
                  if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
                    {
                      /* this one should be ok */
                      npleft--;
                      SP[i].DensityFlag = -1; /* Mark as inactive */
                      continue;
                    }
                } 

              if(SP[i].NgbsMass < (5 * SP[i].MassOfStar))
                Left[i] = dmax(SP[i].Hsml, Left[i]);
              else
                {
                  if(Right[i] != 0)
                    {
                      if(SP[i].Hsml < Right[i])
                        Right[i] = SP[i].Hsml;
                    }
                  else
                    Right[i] = SP[i].Hsml;
                }

              if(Right[i] > 0 && Left[i] > 0)
                SP[i].Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
              else
                {
                  if(Right[i] == 0 && Left[i] == 0)
                    terminate("should not occur");

                  if(Right[i] == 0 && Left[i] > 0)
                    {
                      SP[i].Hsml *= 1.26;
                    }

                  if(Right[i] > 0 && Left[i] == 0)
                    {
                      SP[i].Hsml /= 1.26;
                    }
                }
            }
          else
            SP[i].DensityFlag = -1; /* Mark as inactive */ 
        }

      sumup_large_ints(1, &npleft, &ntot);

      t1 = second();

      if(ntot > 0)
        {
          iter++;

          if(iter > 0)
            mpi_printf("SF_STARBYSTAR: ngb iteration %3d: need to repeat for %12lld particles. (took %g sec)\n", iter, ntot,
                       timediff(t0, t1));

          if(iter > MAXITER)
            terminate("failed to converge in neighbour iteration in sf_starbystar()\n");
        }
    }
  while(ntot > 0);

  myfree(Right);
  myfree(Left);

  /* mark as active again */
  for(i = 0; i < NumStars; i++)
     SP[i].DensityFlag = 1;
  
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
static int sf_starbystar_evaluate(int target, int mode, int threadid)
{
  int i, n, numnodes, *firstnode; 
  double h, h2, dx, dy, dz, r, r2, wk; 
  MyDouble *pos, ngbsmass;

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
  h2 = h * h;

  ngbsmass = 0;

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
          double mu = compute_mu(i); 

          double number_dens = (SphP[i].Density * All.cf_UnitDensity_in_cgs) / mu / PROTONMASS;
          double temp = (SphP[i].Utherm * All.cf_UnitVelocity_in_cm_per_s * All.cf_UnitVelocity_in_cm_per_s) 
          * mu * PROTONMASS * GAMMA_MINUS1 / BOLTZMANN;

          if(number_dens < All.NumberDensThreshold / 10.0 || temp > All.TemperatureThreshold * 10)
            continue;

          //convergence criterion?

          r = sqrt(r2);

          wk = gaussian_weight(r, h);

          // compute the star-ngb-mass 
          ngbsmass += P[i].Mass * wk;
        }
    }

  out.NgbsMass = ngbsmass;

  /* now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

/* \brief Determines if a SP is active in current timestep.
 *
 *  \param[in] n Index of SP in Particle array
 *
 *  \return 1: SP active; 0: SP not active.
 */
int sf_starbystar_isactive(int n)
{
  if(SP[n].DensityFlag < 0)
    return 0;

  return 1;
}