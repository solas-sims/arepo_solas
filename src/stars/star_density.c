#include <stdlib.h>       
#include <math.h>
#include <gsl/gsl_math.h>              
#include <mpi.h>            
  
#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"


/* Pass counter: 1 = find host cell, 2 = gather feedback properties */
static int pass;

static int star_density_evaluate1(int target, int mode, int threadid);
static int star_density_evaluate2(int target, int mode, int threadid);

/* Thin dispatcher - keeps kernel_local/kernel_imported untouched */
static int star_density_evaluate(int target, int mode, int threadid)
{
  if(pass == 1)
    return star_density_evaluate1(target, mode, threadid);
  if(pass == 2)
    return star_density_evaluate2(target, mode, threadid);

  terminate("Star_density_evaluate: invalid pass value %d\n", pass);
  return -1;
}

static int star_density_isactive(int n);

static int feedback_compare(const void *a, const void *b)
{
  const Mechanical_Feedback_Data *da = a;
  const Mechanical_Feedback_Data *db = b;

  if(da->HostIndex < db->HostIndex) return -1;
  if(da->HostIndex > db->HostIndex) return 1;

  if(da->StarTask < db->StarTask) return -1;
  if(da->StarTask > db->StarTask) return 1;

  if(da->StarIndex < db->StarIndex) return -1;
  if(da->StarIndex > db->StarIndex) return 1;

  return 0;
}

static int *StarNgbs;
static int *StarHostIndex;
static int *StarHostTask;
static MyFloat *StarHostDistance;

struct Data 
{
  MyIDType StarID;

  int StarIndex; 
  int StarTask; 
  int HostIndex; 
  int HostTask; 
};

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Pos[3];
  
  struct Data Data;

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  MyDouble TimeToSN;
  MyDouble NextSNEnergy;
#endif  

  Mechanical_Feedback MechanicalFeedback;

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
  
  if(pass == 1)
    {
      in->Data.StarID = -1;

      in->Data.StarIndex = -1;
      in->Data.StarTask = -1;
      in->Data.HostIndex = -1;
      in->Data.HostTask = -1;
    }
  if(pass == 2)
    {
      in->Data.StarID = PPS(i).ID;

      in->Data.StarIndex = i;
      in->Data.StarTask = ThisTask;
      in->Data.HostIndex = StarHostIndex[i];
      in->Data.HostTask = StarHostTask[i];
    }

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  in->TimeToSN = SP[i].TimeToSN;
  in->NextSNEnergy = SP[i].NextSNEnergy;
#endif  

  in->MechanicalFeedback = SP[i].MechanicalFeedback;

  in->Hsml = SP[i].Hsml;
  in->Firstnode = firstnode;
}  

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.
 */
typedef struct
{ 
  /* Pass 1 outputs */
  int Ngbs;
  int HostIndex;
  int HostTask;
  MyFloat HostDistance;

  /* Pass 2 outputs */
  int HostHydroBin;
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
  /* Initial store */
  if(mode == MODE_LOCAL_PARTICLES) 
    {
      /* Pass 1 outputs */
      if(pass == 1)
        {
          StarNgbs[i] = out->Ngbs;
          StarHostIndex[i] = out->HostIndex;
          StarHostTask[i] = out->HostTask;
          StarHostDistance[i] = out->HostDistance;
        }

      /* Pass 2 outputs */
      if(pass == 2)
        {
          SP[i].HostHydroBin = out->HostHydroBin;
        }
    }
  /* Combine */
  else 
    {
      /* Pass 1 outputs */
      if(pass == 1)
        {
          StarNgbs[i] += out->Ngbs;
          if(out->HostDistance < StarHostDistance[i])
            {
              StarHostIndex[i] = out->HostIndex;
              StarHostTask[i] = out->HostTask;
              StarHostDistance[i] = out->HostDistance;
            }
        }

      /* Pass 2 outputs */
      if(pass == 2)
        {
          if(!SP[i].HostHydroBin)
            SP[i].HostHydroBin = out->HostHydroBin;
        }
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
  int idx, i, j;

  int threadid = get_thread_num();

  for(j = 0; j < NTask; j++)
    Thread[threadid].Exportflag[j] = -1;

  while(1)
    {
      if(Thread[threadid].ExportSpace < MinSpace)
        break;

      idx = NextParticle++;

      if(idx >= TimeBinsStar.NActiveParticles)
        break;

      i = TimeBinsStar.ActiveParticleList[idx];

      if(SP[i].WithFeedback == 0)
        continue;
      
      if(star_density_isactive(i))
        star_density_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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
  /* Now do the particles that were sent to us */
  int i, cnt = 0;

  int threadid = get_thread_num();

  while(1)
    {
      i = cnt++;

      if(i >= Nimport)
        break;

      star_density_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
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
void star_density(void)
{
  TIMER_START(CPU_STARS_DENSITY);
  
  int idx, i, npleft, iter = 0;
  long long ntot;
  double t0, t1;

  pass = 0;

  StarNgbs = (int *)mymalloc("StarNgbs", NumStars * sizeof(int));
  StarHostIndex = (int *)mymalloc("StarHostIndex", NumStars * sizeof(int));
  StarHostTask = (int *)mymalloc("StarHostTask", NumStars * sizeof(int));
  StarHostDistance = (MyFloat *)mymalloc("StarHostDistance", NumStars * sizeof(MyFloat));

  memset(StarNgbs, 0, NumStars * sizeof(int));
  memset(StarHostIndex, -1, NumStars * sizeof(int));
  memset(StarHostTask, -1, NumStars * sizeof(int));

  for(i = 0; i < NumStars; i++)
    StarHostDistance[i] = MAX_REAL_NUMBER;

  for(i = 0; i < NumStars; i++)
    {
      if(SP[i].WithFeedback == 0)
        continue;

      SP[i].DensityFlag = 1;
      
      if(SP[i].Hsml <= 0)
        SP[i].Hsml = All.SofteningTable[PPS(i).SofteningType]; 
    }
  
  /* Zero all hosts first */
  for(i = 0; i < NumGas; i++)
    SphP[i].Host = 0;

  generic_set_MaxNexport();
  
  /* Pass 1 - Expand Hsml until we enclose at least one gas cell, then record the closest cell as the host */

  pass++;

  do 
    {
      t0 = second();

      generic_comm_pattern(TimeBinsStar.NActiveParticles, kernel_local, kernel_imported);

      for(idx = 0, npleft = 0; idx < TimeBinsStar.NActiveParticles; idx++)
        {
          i = TimeBinsStar.ActiveParticleList[idx];

          if(SP[i].WithFeedback == 0)
            continue;

          if(StarNgbs[i] < 1)
            {
              npleft++;
              
              SP[i].Hsml *= 2;
            }
          else
            /* Mark as inactive */
            SP[i].DensityFlag = -1; 
        }
     
      sumup_large_ints(1, &npleft, &ntot);

      t1 = second();

      if(ntot > 0)
        {
          iter++;

          if(iter > 0)
            mpi_printf("STAR_DENSITY: ngb iteration %3d: need to repeat for %12lld particles. (took %g sec)\n", iter, ntot,
                       timediff(t0, t1));

          if(iter > MAXITER)
            terminate("failed to converge in neighbour iteration in star_density()\n");
        }
    }
  while(ntot > 0);

  /* Pass 2 - Add star feedback properties to the host cell */

  /* Re-activate all stars */
  for(i = 0; i < NumStars; i++)
    {
      if(SP[i].WithFeedback == 0)
        continue;

      SP[i].DensityFlag = 1;
    }

  pass++;

  generic_comm_pattern(TimeBinsStar.NActiveParticles, kernel_local, kernel_imported);

  /* Sort the hosts list */
  mysort(MechanicalFeedbackEvents.MechanicalFeedbackData, MechanicalFeedbackEvents.NumEvents, 
  sizeof(Mechanical_Feedback_Data), feedback_compare);

  /* Find the total number of events */
  sumup_large_ints(1, &MechanicalFeedbackEvents.NumEvents, &MechanicalFeedbackEvents.TotEvents);

  /* Free arrays */
  myfree(StarHostDistance); 
  myfree(StarHostTask); 
  myfree(StarHostIndex); 
  myfree(StarNgbs);
  
  /* Collect timing information */
  TIMER_STOP(CPU_STARS_DENSITY);
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
static int star_density_evaluate1(int target, int mode, int threadid)
{
  int i, n, numnodes, *firstnode;
  MyDouble xtmp, ytmp, ztmp;  
  MyDouble h, h2, dx, dy, dz, r, r2, wk;
  MyDouble *pos;

  int ngbs = 0, host_index = -1, host_task = -1;
  MyFloat host_distance = MAX_REAL_NUMBER;
  
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
  h2 = h * h;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      i = Thread[threadid].Ngblist[n];

      if(P[i].Type != 0 || P[i].Mass == 0 || P[i].ID == 0)
        continue;

      /* Compute star->cell position vector */
      dx = NEAREST_X(P[i].Pos[0] - pos[0]);
      dy = NEAREST_Y(P[i].Pos[1] - pos[1]);
      dz = NEAREST_Z(P[i].Pos[2] - pos[2]);

      r2 = dx * dx + dy * dy + dz * dz;

      if(r2 < h2)
        {
          ngbs++;

          r = sqrt(r2);
              
          if(r < host_distance)
            {
              host_distance = r;
              host_index = i;
              host_task = ThisTask;
            }
        }
    }

  out.Ngbs = ngbs;
  out.HostIndex = host_index;
  out.HostTask = host_task;
  out.HostDistance = host_distance;

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

static int star_density_evaluate2(int target, int mode, int threadid)
{
  int i, n, numnodes, *firstnode; 
  MyDouble xtmp, ytmp, ztmp;  
  MyDouble h, h2, dx, dy, dz, r, r2, wk;
  MyDouble *pos;

  int hosthydrobin = 0; 
  int star_id, star_index, star_task, host_index, host_task;

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
  h2 = h * h;

  star_id = target_data->Data.StarID;
  star_index = target_data->Data.StarIndex;
  star_task = target_data->Data.StarTask;
  host_index = target_data->Data.HostIndex;
  host_task = target_data->Data.HostTask;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      i = Thread[threadid].Ngblist[n];

      if(P[i].Type != 0 || P[i].Mass == 0 || P[i].ID == 0)
        continue;

      /* Compute star->cell position vector */
      dx = NEAREST_X(P[i].Pos[0] - pos[0]);
      dy = NEAREST_Y(P[i].Pos[1] - pos[1]);
      dz = NEAREST_Z(P[i].Pos[2] - pos[2]);

      r2 = dx * dx + dy * dy + dz * dz;

      if(r2 < h2)
        {
          if(i == host_index && ThisTask == host_task)
            {                
              hosthydrobin = P[i].TimeBinHydro;

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
              /* Limit timestep to resolve the SN */
              MyDouble time_to_sn = target_data->TimeToSN;
              MyDouble next_sn_energy = target_data->NextSNEnergy;

              MyDouble sn_lead_time = All.SN_LeadTime / All.UnitTime_in_Megayears;
          
              if(time_to_sn < sn_lead_time)
                {
                  /* Boost signal speed leading up to an event */
                  double E_inject_code = next_sn_energy;
                 
                  double f = 1.0 - time_to_sn / sn_lead_time;
                  f = fmin(fmax(f, 0.0), 1.0);

                  double Csn = sqrt(GAMMA * GAMMA_MINUS1 * E_inject_code / P[i].Mass) * f;
          
                  if(Csn > SphP[i].Csn)
                    SphP[i].Csn = Csn;
                }
#endif

              /* Setup feedback */
              SphP[i].Host++;

              /* Reallocate events if needed */
              if(MechanicalFeedbackEvents.NumEvents >= MechanicalFeedbackEvents.MaxEvents)
                feedback_reallocate(&MechanicalFeedbackEvents, 2 * MechanicalFeedbackEvents.MaxEvents);

              Mechanical_Feedback_Data *data = &MechanicalFeedbackEvents.MechanicalFeedbackData[MechanicalFeedbackEvents.NumEvents++];

              data->StarID = star_id;
              data->StarIndex = star_index;
              data->StarTask = star_task; 
              data->HostIndex = host_index;
              data->HostTask = host_task;
              data->MechanicalFeedback = target_data->MechanicalFeedback;

              /* Each star has exactly one host - no need to continue */
              break;
            }
        }
    }

  out.HostHydroBin = hosthydrobin;

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
int star_density_isactive(int n)
{
  if(SP[n].DensityFlag < 0)
    return 0;

  return 1;
}