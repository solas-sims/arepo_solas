#include "../main/allvars.h"
#include "../main/proto.h"

#include "../extern/chealpix.h"


double HealpixDirs[MAX_NUM_RAYS][3];

int NRays; 

/* Effective attenuation kappa_ext*(1 - a*<g>) [cm^2/g gas, solar Z]
   Band-averaged over Draine 2003 (renorm. WD01) MW R_V=3.1 model,
   kext_albedo_WD_MW_3.1_60_D03.all, energy and photon-weighted 4e4 K BB.
   Gas mass per H = 2.311e-24 g (M_dust/H = 1.398e-26, M_gas/M_dust = 165.3) */
double Kappa_E[WAVEBANDS] =
{
  [INFRARED] = 34.9, 
  [OPTICAL] = 278.3, 
  [ULTRAVIOLET] = 417.7,  
  [LYMAN_WERNER] = 736.6, /* Dust component */
  [IONIZING_HI] = 1.0,
  [IONIZING_HeI] = 1.0,
  [IONIZING_HeII] = 1.0,
};

double Kappa_N[WAVEBANDS] =
{
  [INFRARED] = 30.0, 
  [OPTICAL] = 242.3, 
  [ULTRAVIOLET] = 406.9,  
  [LYMAN_WERNER] = 731.4, /* Dust component */
  [IONIZING_HI] = 1.0,
  [IONIZING_HeI] = 1.0,
  [IONIZING_HeII] = 1.0,
};

/* Fraction of kappa_eff-attenuated energy that is truly absorbed (heats grains):
   f_abs = kappa_abs/kappa_eff = (1-a)/(1-a<g>), D03 MW dust, band-averaged.
   Remainder is non-forward-scattered light: removed from the ray and it does
   deliver momentum (kappa_eff is exactly the momentum-transfer opacity), but
   it must NOT contribute to heating. */
double TrueAbsorbedFraction[WAVEBANDS] =
{
  [INFRARED] = 0.54,
  [OPTICAL] = 0.62,
  [ULTRAVIOLET] = 0.81,
  [LYMAN_WERNER] = 0.88, /* Dust share only; H2 line share is pure absorption */
  [IONIZING_HI] = 1.0,
  [IONIZING_HeI] = 1.0,
  [IONIZING_HeII] = 1.0,
};

double ReradiatedFraction[WAVEBANDS] = 
{
  [INFRARED] = 0.54,
  [OPTICAL] = 0.62,
  [ULTRAVIOLET] = 0.77, /* 5% goes to pe heating */
  [LYMAN_WERNER] = 0.84, /* 5% goes to pe heating */
  [IONIZING_HI] = 0.0, /* No reradiation */
  [IONIZING_HeI] = 0.0, /* No reradiation */
  [IONIZING_HeII] = 0.0, /* No reradiation */
};

void update_dtau(void)
{
  for(int i = 0; i < NumGas; i++)
    {
      if(P[i].Type != 0 || P[i].Mass == 0 || P[i].ID == 0)
        continue;

      double Units = All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm / All.cf_UnitMass_in_g;
      
#ifdef METALS
      double Zsol = ((SphP[i].GasMetals + SphP[i].StarMetalsFeed) / (P[i].Mass + SphP[i].StarMassFeed)) / SOLAR_METALLICITY;
#else
      double Zsol = 0;
#endif

      double Density = (P[i].Mass + SphP[i].StarMassFeed) / SphP[i].Volume;

      double sigma_HI, sigma_HeI, sigma_HeII;

      SphP[i].DtauOverLength_E[INFRARED] = (Kappa_E[INFRARED] / Units) * Zsol * Density;
      SphP[i].DtauOverLength_E[OPTICAL] = (Kappa_E[OPTICAL] / Units) * Zsol * Density;
      SphP[i].DtauOverLength_E[ULTRAVIOLET] = (Kappa_E[ULTRAVIOLET] / Units) * Zsol * Density;      

      /* Dust only */
      SphP[i].DtauOverLength_E[LYMAN_WERNER] = (Kappa_E[LYMAN_WERNER] / Units) * Zsol * Density;
      
      /* Band averaged sigma0*(v/v0)^(-3), in cm^2 */
      sigma_HI = 3.25e-18;
      sigma_HeI = 5.04e-18;
      sigma_HeII = 1.30e-18; 
      SphP[i].DtauOverLength_E[IONIZING_HI] = (sigma_HI / (PROTONMASS) / Units) * SphP[i].GrackleSpeciesConserved(GRACKLE_HI) / SphP[i].Volume;
      SphP[i].DtauOverLength_E[IONIZING_HeI] = (sigma_HeI / (4.0 * PROTONMASS) / Units) * SphP[i].GrackleSpeciesConserved(GRACKLE_HeI) / SphP[i].Volume;
      SphP[i].DtauOverLength_E[IONIZING_HeII] = (sigma_HeII / (4.0 * PROTONMASS) / Units) * SphP[i].GrackleSpeciesConserved(GRACKLE_HeII) / SphP[i].Volume;  
      

      SphP[i].DtauOverLength_N[INFRARED] = (Kappa_N[INFRARED] / Units) * Zsol * Density;
      SphP[i].DtauOverLength_N[OPTICAL] = (Kappa_N[OPTICAL] / Units) * Zsol * Density;
      SphP[i].DtauOverLength_N[ULTRAVIOLET] = (Kappa_N[ULTRAVIOLET] / Units) * Zsol * Density;      

      /* Dust only */
      SphP[i].DtauOverLength_N[LYMAN_WERNER] = (Kappa_N[LYMAN_WERNER] / Units) * Zsol * Density;
      
      /* Band averaged sigma0*(v/v0)^(-3), in cm^2 */
      sigma_HI = 3.48e-18; 
      sigma_HeI = 5.27e-18;
      sigma_HeII = 1.31e-18; 
      SphP[i].DtauOverLength_N[IONIZING_HI] = (sigma_HI / (PROTONMASS) / Units) * SphP[i].GrackleSpeciesConserved(GRACKLE_HI) / SphP[i].Volume;
      SphP[i].DtauOverLength_N[IONIZING_HeI] = (sigma_HeI / (4.0 * PROTONMASS) / Units) * SphP[i].GrackleSpeciesConserved(GRACKLE_HeI) / SphP[i].Volume;
      SphP[i].DtauOverLength_N[IONIZING_HeII] = (sigma_HeII / (4.0 * PROTONMASS) / Units) * SphP[i].GrackleSpeciesConserved(GRACKLE_HeII) / SphP[i].Volume;  
    }
}

double dtau_IR(int i, double length)
{
  double Units = All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm / All.cf_UnitMass_in_g;
      
#ifdef METALS
  double Zsol = ((SphP[i].GasMetals + SphP[i].StarMetalsFeed) / (P[i].Mass + SphP[i].StarMassFeed)) / SOLAR_METALLICITY;
#else
  double Zsol = 0;
#endif

 double Density = (P[i].Mass + SphP[i].StarMassFeed) / SphP[i].Volume;

  double Dtau_IR = (1.0 / Units) * Zsol * Density * length;

  return Dtau_IR;
}

static double H2Tab_A[H2TAB_N]; /* A at table nodes */
static double H2Tab_dlogN; /* log10 spacing */
static double H2Tab_A_thinmin; /* A(NMIN) = SIGMA_PUMP * NMIN */
 
/* Wolcott-Green et al. (2011) self-shielding function */
static inline double f_selfshield_H2(double N_H2)
{
  double x  = N_H2 / 5.0e14;
  double sq = sqrt(1.0 + x);
  return 0.965 / pow(1.0 + x / H2_SHIELD_B5, 1.1)
       + 0.035 / sq * exp(-8.5e-4 * sq);
}
 
/* Build A(N) once at startup (trapezoid, log-spaced with linear-N areas,
   16 sub-steps per interval so table error << fit error) */
void init_h2shield(void)
{
  H2Tab_dlogN = (H2TAB_LOGNMAX - H2TAB_LOGNMIN) / (H2TAB_N - 1);
  H2Tab_A_thinmin = SIGMA_PUMP * pow(10.0, H2TAB_LOGNMIN);
 
  /* Thin part below NMIN: f_sh=1 */
  double A = H2Tab_A_thinmin;          
  H2Tab_A[0] = A;
 
  for(int i = 1; i < H2TAB_N; i++)
    {
      double N0 = pow(10.0, H2TAB_LOGNMIN + (i - 1) * H2Tab_dlogN);
      double N1 = pow(10.0, H2TAB_LOGNMIN + i * H2Tab_dlogN);
 
      const int nsub = 16;
      double dN = (N1 - N0) / nsub;
      for(int k = 0; k < nsub; k++)
        {
          double Na = N0 + k * dN;
          A += 0.5 * (f_selfshield_H2(Na) + f_selfshield_H2(Na + dN)) * dN * SIGMA_PUMP;
        }
      H2Tab_A[i] = A;
    }
}
 
/* A(N): thin analytic below NMIN, clamp above NMAX, linear-in-logN inside */
static inline double h2shield_A(double N_H2)
{
  if(N_H2 <= 0.0)
    return 0.0;
 
  double logN = log10(N_H2);
 
  /* f_sh = 1 exactly */
  if(logN <= H2TAB_LOGNMIN)
    return SIGMA_PUMP * N_H2; 
 
  /* Lines exhausted */
  if(logN >= H2TAB_LOGNMAX)
    return H2Tab_A[H2TAB_N - 1];     
 
  double u = (logN - H2TAB_LOGNMIN) / H2Tab_dlogN;
  int    j = (int)u;
  double f = u - j;
 
  return H2Tab_A[j] * (1.0 - f) + H2Tab_A[j + 1] * f;
}
 
/* Exact per-cell line optical depth for a cell adding dN_H2 to a ray
   that has already accumulated N_H2 */
double h2shield_dtau(double N_H2, double dN_H2)
{
  double N_H2_cgs = N_H2 * (All.cf_UnitMass_in_g / (All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm)) / (2.0 * PROTONMASS);
  double dN_H2_cgs = dN_H2 * (All.cf_UnitMass_in_g / (All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm)) / (2.0 * PROTONMASS);

  double dtau = h2shield_A(N_H2_cgs + dN_H2_cgs) - h2shield_A(N_H2_cgs);
  return dtau > 0.0 ? dtau : 0.0;
}

void start_healpix(void) 
{
  int nside = NSIDE_MIN;
  NRays = 12 * nside * nside;

  for(int ipix = 0; ipix < NRays; ipix++)
    {
      pix2vec_nest(nside, ipix, HealpixDirs[ipix]);
    }
}

static RayWorkStack *init_work_stack(long long capacity)
{
  RayWorkStack *w = malloc(sizeof(RayWorkStack));

  w->n = 0;
  w->capacity = capacity;
  w->rays = malloc(capacity * sizeof(RayPacket));
  return w;
}

void append_ray(RayWorkStack *w, const RayPacket *ray)
{
  if(w->n >= w->capacity)
    {
      w->capacity *= 2;
      w->rays = realloc(w->rays, w->capacity * sizeof(RayPacket));
    }
  w->rays[w->n++] = *ray;
}

static void free_work_stack(RayWorkStack *w)
{
  free(w->rays); free(w);
}

static void init_rays(RayWorkStack *work)
{
  double SQRT3 = sqrt(3);
    
  int ray_idx = 0;
     
  /* Act on host cells */
  for(int ev = 0; ev < MechanicalFeedbackEvents.NumEvents; ev++)
    {
      Mechanical_Feedback_Data *MechanicalFeedbackData = &MechanicalFeedbackEvents.MechanicalFeedbackData[ev];
      Mechanical_Feedback *MechanicalFeedback = &MechanicalFeedbackData->MechanicalFeedback;
        
      /* Loop over rays for this star */
      for(int iray = 0; iray < NRays; iray++)
        {  
          /* Initialize ray from star i */
          RayPacket ray = {0};

          ray.pos[0] = MechanicalFeedback->StarPosition[0];      
          ray.pos[1] = MechanicalFeedback->StarPosition[1]; 
          ray.pos[2] = MechanicalFeedback->StarPosition[2]; 
          ray.dir[0] = HealpixDirs[iray][0];        
          ray.dir[1] = HealpixDirs[iray][1];
          ray.dir[2] = HealpixDirs[iray][2];
          ray.t = 0.0;
          ray.t_exit = MAX_REAL_NUMBER;
          ray.t_maximum = SQRT3 * All.BoxSize;

          ray.active_bands = NO_IR_ACTIVE;

          for(int w = 0; w < WAVEBANDS; w++)
            { 
              ray.Radiated[w].Energy = MechanicalFeedback->Radiated[w].Energy / NRays;
              ray.Radiated[w].Photons = MechanicalFeedback->Radiated[w].Photons / NRays;

              ray.Radiated_Init[w].Energy = MechanicalFeedback->Radiated[w].Energy / NRays;
              ray.Radiated_Init[w].Photons = MechanicalFeedback->Radiated[w].Photons / NRays;

              if(ray.Radiated[w].Energy <= 0.0 && ray.Radiated[w].Photons <= 0.0)
                ray.active_bands &= (uint8_t)(~(1u << w));
            }
          
          ray.ray_id = ray_idx;
          ray.home_task = ThisTask;
          
          ray.n_pending = 0;
          ray.target_node = -1;

          ray.is_paused = 0;

          ray.nside = NSIDE_MIN;        
          ray.healpix_pixel = iray;            
            
          ray_idx++;

          append_ray(work, &ray);
        }
    }
}

/* Splits to 4 child rays */
void split_ray(const RayPacket *parent, RayPacket children[4])
{
  int new_nside = parent->nside * 2;

  for(int k = 0; k < 4; k++)
    {
      /* Copy all state including t, active_bands etc */
      children[k] = *parent;   

      children[k].nside = new_nside;
      children[k].healpix_pixel = 4 * parent->healpix_pixel + k;

      pix2vec_nest(new_nside, children[k].healpix_pixel, children[k].dir);

      for(int w = 0; w < WAVEBANDS; w++)
        {
          children[k].Radiated[w].Energy = parent->Radiated[w].Energy * 0.25;
          children[k].Radiated[w].Photons = parent->Radiated[w].Photons * 0.25;

          children[k].Radiated_Init[w].Energy = parent->Radiated_Init[w].Energy * 0.25;
          children[k].Radiated_Init[w].Photons = parent->Radiated_Init[w].Photons * 0.25;
        }
    }
}

static RayExportBuffer *init_export_buffer(long long capacity)
{
  RayExportBuffer *buf = malloc(sizeof(RayExportBuffer));
 
  buf->n = 0;
  buf->capacity = capacity;
  buf->task = malloc(capacity * sizeof(int));
  buf->rays = malloc(capacity * sizeof(RayPacket));

  return buf;
}

void append_export(RayExportBuffer *buf, const RayPacket *ray, int task)
{
  if(buf->n >= buf->capacity)
    {
      buf->capacity *= 2;
      buf->task = realloc(buf->task, buf->capacity * sizeof(int));
      buf->rays = realloc(buf->rays, buf->capacity * sizeof(RayPacket));
    }

  buf->task[buf->n] = task;
  buf->rays[buf->n] = *ray;
  buf->n++;
}

static void free_export_buffer(RayExportBuffer *buf)
{
  free(buf->rays); free(buf->task); free(buf);
}

static void sort_by_task(RayExportBuffer *buf)
{
  int *count = calloc(NTask, sizeof(int));
  int *cursor = malloc(NTask * sizeof(int));

  for(long long i = 0; i < buf->n; i++)
    count[buf->task[i]]++;

  cursor[0] = 0;
  for(int t = 1; t < NTask; t++)
    cursor[t] = cursor[t-1] + count[t-1];

  int *sorted_task = malloc(buf->n * sizeof(int));
  RayPacket *sorted_rays = malloc(buf->n * sizeof(RayPacket));
 
  for(long long i = 0; i < buf->n; i++)
    {
      int t = buf->task[i];
      long long pos = cursor[t]++;
      sorted_task[pos] = t;
      sorted_rays[pos] = buf->rays[i]; 
    }

  memcpy(buf->task, sorted_task, buf->n * sizeof(int));
  memcpy(buf->rays, sorted_rays, buf->n * sizeof(RayPacket));

  free(sorted_rays); 
  free(sorted_task); 
  free(cursor); 
  free(count); 
}

static void exchange_rays(RayExportBuffer *send, RayWorkStack *work)
{
  static MPI_Datatype MPI_RAYPACKET = MPI_DATATYPE_NULL;
  if(MPI_RAYPACKET == MPI_DATATYPE_NULL)
    {
      MPI_Type_contiguous(sizeof(RayPacket), MPI_BYTE, &MPI_RAYPACKET);
      MPI_Type_commit(&MPI_RAYPACKET);
    }

  /* Heap, not VLA: int send_count[NTask] blows the C stack for large NTask */
  int *send_count = malloc(NTask * sizeof(int));
  int *recv_count = malloc(NTask * sizeof(int));
  int *send_offset = malloc(NTask * sizeof(int));
  int *recv_offset = malloc(NTask * sizeof(int));

  memset(send_count, 0, NTask * sizeof(int));
  for(long long i = 0; i < send->n; i++)
    send_count[send->task[i]]++;

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  send_offset[0] = recv_offset[0] = 0;
  long long total_recv = recv_count[0];
  for(int i = 1; i < NTask; i++)
    {
      send_offset[i] = send_offset[i-1] + send_count[i-1];
      recv_offset[i] = recv_offset[i-1] + recv_count[i-1];
      total_recv += recv_count[i];
    }

  while(work->n + total_recv > work->capacity)
    {
      work->capacity *= 2;
      work->rays = realloc(work->rays, work->capacity * sizeof(RayPacket));
    }

  /* In ray Units */
  sort_by_task(send);   

  MPI_Alltoallv(send->rays, send_count, send_offset, MPI_RAYPACKET,
  work->rays + work->n, recv_count, recv_offset, MPI_RAYPACKET,
  MPI_COMM_WORLD);

  work->n += total_recv;

  free(recv_offset); 
  free(send_offset); 
  free(recv_count); 
  free(send_count);
}

#ifdef TREEPOINTS
/* Not operational yet -Ngb tree does not use Tree_Points */
struct rad_resultsactiveimported_data *Rad_ResultsActiveImported;

static void send_results_home(void)
{
  int i, j, n, k, ncount;
  int *Recv_count = malloc(NTask * sizeof(int));
  int *Send_count = malloc(NTask * sizeof(int));
  int *Recv_offset = malloc(NTask * sizeof(int));
  int *Send_offset = malloc(NTask * sizeof(int));

  /* Count gas cells among imported particles */
  for(i = 0, ncount = 0; i < Tree_NumPartImported; i++)
    if(Ngb_Tree_Points[i].Type == 0)
      ncount++;

  Rad_ResultsActiveImported = malloc(ncount * sizeof(struct rad_resultsactiveimported_data));

  /* Pack gas cell results */
  for(j = 0; j < NTask; j++)
    Recv_count[j] = 0;

  for(i = 0, n = 0, k = 0; i < NTask; i++)
    for(j = 0; j < Force_Recv_count[i]; j++, n++)
      if(Ngb_Tree_Points[n].Type == 0)
        {
          Rad_ResultsActiveImported[k].StarMomentumFeed[0] = Ngb_Tree_Points[n].StarMomentumFeed[0];
          Rad_ResultsActiveImported[k].StarMomentumFeed[1] = Ngb_Tree_Points[n].StarMomentumFeed[1];
          Rad_ResultsActiveImported[k].StarMomentumFeed[2] = Ngb_Tree_Points[n].StarMomentumFeed[2];

          for(int w = 0; w < WAVEBANDS; w++)
            {
              Rad_ResultsActiveImported[k].Absorbed[w].Energy = Ngb_Tree_Points[n].Absorbed[w].Energy;
              Rad_ResultsActiveImported[k].Absorbed[w].Photons = Ngb_Tree_Points[n].Absorbed[w].Photons;
            }

          Rad_ResultsActiveImported[k].index = Ngb_Tree_Points[n].index;
          Recv_count[i]++;
          k++;
        }

  MPI_Alltoall(Recv_count, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);

  int Nexport = 0, Nimport = 0;
  Send_offset[0] = Recv_offset[0] = 0;
  for(j = 0; j < NTask; j++)
    {
      Nexport += Send_count[j];
      Nimport += Recv_count[j];
      if(j > 0)
        {
          Send_offset[j] = Send_offset[j-1] + Send_count[j-1];
          Recv_offset[j] = Recv_offset[j-1] + Recv_count[j-1];
        }
    }

  struct rad_resultsactiveimported_data *tmp_results = malloc(Nexport * sizeof(struct rad_resultsactiveimported_data));

  /* Exchange results back to home ranks */
  for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&Rad_ResultsActiveImported[Recv_offset[recvTask]],
          Recv_count[recvTask] * sizeof(struct rad_resultsactiveimported_data), MPI_BYTE, recvTask, TAG_RAD,
          &tmp_results[Send_offset[recvTask]],
          Send_count[recvTask] * sizeof(struct rad_resultsactiveimported_data), MPI_BYTE, recvTask, TAG_RAD,
          MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  /* Apply results to local particles */
  for(i = 0; i < Nexport; i++)
    {  
      SphP[tmp_results[i].index].StarMomentumFeed[0] += tmp_results[i].StarMomentumFeed[0];
      SphP[tmp_results[i].index].StarMomentumFeed[1] += tmp_results[i].StarMomentumFeed[1];
      SphP[tmp_results[i].index].StarMomentumFeed[2] += tmp_results[i].StarMomentumFeed[2];
      
      for(int w = 0; w < WAVEBANDS; w++)
        {
          SphP[tmp_results[i].index].Absorbed[w].Energy += tmp_results[i].Absorbed[w].Energy;
          SphP[tmp_results[i].index].Absorbed[w].Photons += tmp_results[i].Absorbed[w].Photons;
        }
    }

  /* Free in reverse allocation order */
  free(tmp_results);
  free(Rad_ResultsActiveImported);
  free(Send_offset); free(Recv_offset);
  free(Send_count); free(Recv_count);
}
#endif

#ifdef RAD_OPENING_ANGLE
static void distribute_node_rad(int no)
{
  /* Quick check - skip empty nodes */
  int has_rad = 0;
        
  for(int w = 0; w < WAVEBANDS; w++)
    if(RtNgb_Nodes[no].Absorbed[w].Energy > 0 || RtNgb_Nodes[no].Absorbed[w].Photons > 0) 
      { 
        has_rad = 1; 
        break; 
      }
      
  if(!has_rad) 
    return;
  
  double node_dtau_E[WAVEBANDS], node_dtau_N[WAVEBANDS];

  for(int w = 0; w < WAVEBANDS; w++)
    {
      if(w == LYMAN_WERNER)
        {
          node_dtau_E[w] = RtNgb_Nodes[no].Volume * RtNgb_Nodes[no].DtauOverLength_E[w];
          node_dtau_N[w] = RtNgb_Nodes[no].Volume * RtNgb_Nodes[no].dN_H2_OverLength;
        }
      else
        {
          node_dtau_E[w] = RtNgb_Nodes[no].Volume * RtNgb_Nodes[no].DtauOverLength_E[w];
          node_dtau_N[w] = RtNgb_Nodes[no].Volume * RtNgb_Nodes[no].DtauOverLength_N[w];
        }
    }

  int child = Ngb_Nodes[no].u.d.nextnode;
  while(child != Ngb_Nodes[no].u.d.sibling && child >= 0)
    {
      /* Leaf particle - deposit directly */
      if(child < Ngb_MaxPart)
        {
          if(P[child].Type != 0 || P[child].Mass == 0 || P[child].ID == 0)
            {
              child = Ngb_Nextnode[child];
              continue;
            }
          
          for(int w = 0; w < WAVEBANDS; w++)
            {
              if(node_dtau_E[w] > 0)
                {
                  double child_dtau_E;
                  if(w == LYMAN_WERNER)
                    {
                       child_dtau_E = SphP[child].Volume * SphP[child].DtauOverLength_E[w];
                    }
                  else
                    {
                       child_dtau_E = SphP[child].Volume * SphP[child].DtauOverLength_E[w];
                    }
                  
                  double frac_E = child_dtau_E / node_dtau_E[w];
                  
                  SphP[child].Absorbed[w].Energy += frac_E * RtNgb_Nodes[no].Absorbed[w].Energy;
                }
              
              if(node_dtau_N[w] > 0)
                {
                  double child_dtau_N;
                  if(w == LYMAN_WERNER)
                    {
                      child_dtau_N = SphP[child].Volume * SphP[child].GrackleSpeciesConserved(GRACKLE_H2I) / SphP[child].Volume;
                    }
                  else
                    {
                      child_dtau_N = SphP[child].Volume * SphP[child].DtauOverLength_N[w];
                    }

                  double frac_N = child_dtau_N / node_dtau_N[w];
                  
                  SphP[child].Absorbed[w].Photons += frac_N * RtNgb_Nodes[no].Absorbed[w].Photons;
                }
            }

          child = Ngb_Nextnode[child];
        }
      /* Internal node - pass fraction down recursively */
      else if(child < Ngb_MaxPart + Ngb_MaxNodes)
        {
          for(int w = 0; w < WAVEBANDS; w++)
            {
              if(node_dtau_E[w] > 0)
                {
                  double child_dtau_E;
                  if(w == LYMAN_WERNER)
                    {
                      child_dtau_E = RtNgb_Nodes[child].Volume * RtNgb_Nodes[child].DtauOverLength_E[w];
                    }
                  else
                    {
                      child_dtau_E = RtNgb_Nodes[child].Volume * RtNgb_Nodes[child].DtauOverLength_E[w];
                    }
                  
                  double frac_E = child_dtau_E / node_dtau_E[w];
                  
                  RtNgb_Nodes[child].Absorbed[w].Energy += frac_E * RtNgb_Nodes[no].Absorbed[w].Energy;
                }
              
              if(node_dtau_N[w] > 0)
                {
                  double child_dtau_N;
                  if(w == LYMAN_WERNER)
                    {
                      child_dtau_N = RtNgb_Nodes[child].Volume * RtNgb_Nodes[child].dN_H2_OverLength;
                    }
                  else
                    {
                      child_dtau_N = RtNgb_Nodes[child].Volume * RtNgb_Nodes[child].DtauOverLength_N[w];
                    }
                  
                  double frac_N = child_dtau_N / node_dtau_N[w];
                  
                  RtNgb_Nodes[child].Absorbed[w].Photons += frac_N * RtNgb_Nodes[no].Absorbed[w].Photons;
                }
            }

          distribute_node_rad(child);
          child = Ngb_Nodes[child].u.d.sibling;
        }
      /* Pseudo-particle - skip, handled by export */
      else
        {    
          child = Ngb_Nextnode[child - Ngb_MaxNodes];
        }
    }

  for(int w = 0; w < WAVEBANDS; w++)
    RtNgb_Nodes[no].Absorbed[w].Energy = RtNgb_Nodes[no].Absorbed[w].Photons = 0.0;
}
#endif

static void radiation_feedback(void)
{
  /* Photoionization and photoelectric heating here -> we do rad pressure inside the tree walk */
  int i, w;
  
  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Type != 0 || P[i].Mass == 0 || P[i].ID == 0)
        continue;
      
      /* Volume, timestep */
      double V = SphP[i].Volume;
      double dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval; 

      /* In cgs */
      double V_cgs = V * (All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm);
      double dt_cgs = dt * All.cf_UnitTime_in_s;

#ifdef PHOTOELECTRIC_HEATING
      /* Photoelectric heating */
      double epsilon_pe = 0.05;

      double E_pe = (SphP[i].Absorbed[ULTRAVIOLET].Energy  * TrueAbsorbedFraction[ULTRAVIOLET]
                     + SphP[i].Absorbed[LYMAN_WERNER].Energy * TrueAbsorbedFraction[LYMAN_WERNER])
                     * epsilon_pe * All.cf_UnitEnergy_in_cgs;
      
      /* Volumetric_heating_rate: grackle docs say erg/(s cm^3), straight CGS, no conversion */
      SphP[i].PE_VolHeatingRate +=  E_pe / dt_cgs / V_cgs;
#endif

#ifdef DISSOCIATION
      /* H2 Dissociation */
      /* Number density */
      double n_H2 = SphP[i].GrackleSpeciesConserved(GRACKLE_H2I) / V / (2 * PROTONMASS / All.cf_UnitMass_in_g);

      double N_abs_H2 = SphP[i].Absorbed[LYMAN_WERNER].Photons;

      SphP[i].H2_DissociationRate += n_H2 > 0 ? (F_DISS * N_abs_H2 / (dt/All.cf_hubble_a/All.HubbleParam) / V) / n_H2 : 0.0;
#endif

#ifdef PHOTOIONIZATION
      /* Photoionization */     
      /* Number densities */
      double n_HI = SphP[i].GrackleSpeciesConserved(GRACKLE_HI) / V / (PROTONMASS / All.cf_UnitMass_in_g);
      double n_HeI = SphP[i].GrackleSpeciesConserved(GRACKLE_HeI) / V / (4 * PROTONMASS / All.cf_UnitMass_in_g);
      double n_HeII = SphP[i].GrackleSpeciesConserved(GRACKLE_HeII) / V / (4 * PROTONMASS / All.cf_UnitMass_in_g);

      /* In cgs */
      double n_HI_cgs = n_HI / (All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm);
      double n_HeI_cgs = n_HeI / (All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm);
      double n_HeII_cgs = n_HeII / (All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm);
     
      /* Threshold energies */ 
      double energy_thresh_HI = 13.6 * ELECTRONVOLT_IN_ERGS;
      double energy_thresh_HeI = 24.6 * ELECTRONVOLT_IN_ERGS;
      double energy_thresh_HeII = 54.4 * ELECTRONVOLT_IN_ERGS;

      double E_abs_HI = SphP[i].Absorbed[IONIZING_HI].Energy * All.cf_UnitEnergy_in_cgs;
      double E_abs_HeI = SphP[i].Absorbed[IONIZING_HeI].Energy * All.cf_UnitEnergy_in_cgs;
      double E_abs_HeII = SphP[i].Absorbed[IONIZING_HeII].Energy * All.cf_UnitEnergy_in_cgs;      

      double N_abs_HI = SphP[i].Absorbed[IONIZING_HI].Photons;
      double N_abs_HeI = SphP[i].Absorbed[IONIZING_HeI].Photons;
      double N_abs_HeII = SphP[i].Absorbed[IONIZING_HeII].Photons;

      /* RT_heating_rate: grackle docs say erg/(s cm^3) / n, straight CGS, no conversion */
      double E_threshold_HI = N_abs_HI * energy_thresh_HI; 
      
      if(n_HI_cgs)
        SphP[i].HI_HeatingRate += (E_abs_HI - E_threshold_HI) > 0 ? ((E_abs_HI - E_threshold_HI) / dt_cgs / V_cgs) / n_HI_cgs : 0.0;
      
      double E_threshold_HeI = N_abs_HeI * energy_thresh_HeI;
      
      if(n_HeI_cgs) 
        SphP[i].HeI_HeatingRate += (E_abs_HeI - E_threshold_HeI) > 0 ? ((E_abs_HeI - E_threshold_HeI) / dt_cgs / V_cgs) / n_HeI_cgs : 0.0;
      
      double E_threshold_HeII = N_abs_HeII * energy_thresh_HeII;
      
      if(n_HeII_cgs)
        SphP[i].HeII_HeatingRate += (E_abs_HeII - E_threshold_HeII) > 0 ? ((E_abs_HeII - E_threshold_HeII) / dt_cgs / V_cgs) / n_HeII_cgs : 0.0;

      /* RT_ionization_rate: 1 / (time Units) */
      SphP[i].HI_IonizationRate += n_HI > 0 ? (N_abs_HI / (dt/All.cf_hubble_a/All.HubbleParam) / V) / n_HI : 0.0;     
      
      SphP[i].HeI_IonizationRate += n_HeI > 0 ? (N_abs_HeI / (dt/All.cf_hubble_a/All.HubbleParam) / V) / n_HeI : 0.0;

      SphP[i].HeII_IonizationRate += n_HeII > 0 ? (N_abs_HeII / (dt/All.cf_hubble_a/All.HubbleParam) / V) / n_HeII : 0.0;
#endif

      for(w = 0; w < WAVEBANDS; w++)
        SphP[i].Absorbed[w].Energy = SphP[i].Absorbed[w].Photons = 0.0;
    }
}

#ifdef PHOTOIONIZATION
static void rt_timestep(void)
{
  int idx, i;
  double eps_ion = All.RTIonizationTimestepFraction;
  
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      /* Normalize ionization rate by total hydrogen (was normalized by neutral for grackle) */
      double m_HI = SphP[i].GrackleSpeciesConserved(GRACKLE_HI);

      double m_H = SphP[i].GrackleSpeciesConserved(GRACKLE_HI) + SphP[i].GrackleSpeciesConserved(GRACKLE_HII);

#if GRACKLE_CHEMISTRY >= 2
      m_H += SphP[i].GrackleSpeciesConserved(GRACKLE_H2I) + SphP[i].GrackleSpeciesConserved(GRACKLE_H2II) 
          + SphP[i].GrackleSpeciesConserved(GRACKLE_HM);
#endif

#if GRACKLE_CHEMISTRY >= 3
      m_H += SphP[i].GrackleSpeciesConserved(GRACKLE_DI) + SphP[i].GrackleSpeciesConserved(GRACKLE_DII) 
          + SphP[i].GrackleSpeciesConserved(GRACKLE_HDI);
#endif
    
      if(m_H > 0)
        {
          double x_HI = m_HI / m_H;
          rate = fmax(rate, SphP[i].HI_IonizationRate * x_HI);
        }

      SphP[i].RT_Timestep = (rate > 0.0) ? eps_ion / rate : All.MaxSizeTimestep / All.cf_hubble_a;
    }
}
#endif

void star_radiation(void)
{
  TIMER_START(CPU_STARS_RADIATION);

  double t0, t1;

  update_dtau();

#ifdef RAD_OPENING_ANGLE
  /* Zero accumulator on all nodes before treewalk - important for top nodes! */
  for(int no = Ngb_MaxPart; no < Ngb_MaxPart + Ngb_NumNodes; no++)
    {
      for(int w = 0; w < WAVEBANDS; w++)
        {
          RtNgb_Nodes[no].Absorbed[w].Energy = RtNgb_Nodes[no].Absorbed[w].Photons = 0.0;
        }
    }
#endif

  /* Zero accumulator on all leaves before treewalk */
  for(int i = 0; i < NumGas; i++)
    {
      for(int w = 0; w < WAVEBANDS; w++)
        {
          SphP[i].Absorbed[w].Energy = SphP[i].Absorbed[w].Photons = 0.0;
        }
    }
 
  int n_stars = MechanicalFeedbackEvents.NumEvents;
  long long n_rays_local = (long long)n_stars * NRays;

  long long n_rays_global;
  sumup_longs(1, &n_rays_local, &n_rays_global);

  mpi_printf("STAR_RADIATION: Initializing radiation with %12lld rays\n", n_rays_global);

  /* Floor so ranks with no local stars still have a buffer to receive imports */
  long long work_capacity = n_rays_local > 0 ? 4 * n_rays_local : 1024;
  long long export_capacity = n_rays_local > 0 ? n_rays_local : 1024;

  RayWorkStack *work = init_work_stack(work_capacity);
  RayExportBuffer *export_buf = init_export_buffer(export_capacity);
  
  init_rays(work);

  long long n_global;
  int iter = 0;
  do
    {
      t0 = second();

      while(work->n > 0)
        {
          RayPacket ray = work->rays[--work->n];
          raytrace_treewalk(&ray, work, export_buf);
        }

      /* Send rays to remote ranks, receive rays from remote ranks */
      exchange_rays(export_buf, work);

      /* Reset export buffer for this round */
      export_buf->n = 0;

      /* Check if anyone still has rays in flight */
      sumup_longs(1, &work->n, &n_global);

      t1 = second();

      iter++;

      if(n_global > 0 && iter > 0)
        mpi_printf("STAR_RADIATION: Rad iteration %3d: need to repeat for %12lld rays. (took %g sec)\n", iter, n_global,
        timediff(t0, t1));
      
    } while(n_global > 0);
    
  //send_results_home();

#ifdef RAD_OPENING_ANGLE
  distribute_node_rad(Ngb_MaxPart);
#endif

  radiation_feedback();

#ifdef PHOTOIONIZATION
  rt_timestep();
#endif

  free_export_buffer(export_buf);
  free_work_stack(work);

  TIMER_STOP(CPU_STARS_RADIATION);
}