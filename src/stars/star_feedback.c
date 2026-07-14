#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

/* SN host-injection modes, returned by SN_feedback_radius() */
#define SN_MESH          0 /* Couple across Voronoi faces */
#define SN_HOST          1 /* Thermal dump into host  */
#define SN_HOST_MOMENTUM 2 /* Momentum kick to host */

/* Kick packet sent to remote face-neighbor cells */
struct Feedback_Kick
{ 
  /* Cell index on the receiving task */
  int CellIndex; 

#ifdef WINDS
  MyDouble DeltaMass;
#ifdef METALS
  MyDouble DeltaMetals;
#endif
  MyDouble DeltaP[3];
  MyDouble DeltaE;
#endif

#ifdef SUPERNOVAE
  MyDouble SN_DeltaMass;
#ifdef METALS
  MyDouble SN_DeltaMetals;
#endif
  MyDouble SN_DeltaP[3];
  MyDouble SN_DeltaE;
#endif
};

/* Apply a kick packet to a local cell */
static void apply_kick(int j, const struct Feedback_Kick *Kick)
{
#ifdef WINDS
  SphP[j].StarMassFeed += Kick->DeltaMass;
  All.StarFeedbackLocal[0] += Kick->DeltaMass;
#ifdef METALS
  SphP[j].StarMetalsFeed += Kick->DeltaMetals;
  All.StarFeedbackLocal[1] += Kick->DeltaMetals;
#endif
  for(int k = 0; k < 3; k++)
    SphP[j].StarMomentumFeed[k] += Kick->DeltaP[k];

  SphP[j].StarEnergyFeed += Kick->DeltaE;
  All.StarFeedbackLocal[2] += Kick->DeltaE;
#endif 

#ifdef SUPERNOVAE
  SphP[j].StarMassFeed += Kick->SN_DeltaMass;
  All.StarFeedbackLocal[0] += Kick->SN_DeltaMass;
#ifdef METALS
  SphP[j].StarMetalsFeed += Kick->SN_DeltaMetals;
  All.StarFeedbackLocal[1] += Kick->SN_DeltaMetals;
#endif
  for(int k = 0; k < 3; k++)
    SphP[j].StarMomentumFeed[k] += Kick->SN_DeltaP[k];

  SphP[j].StarEnergyFeed += Kick->SN_DeltaE;
  All.StarFeedbackLocal[2] += Kick->SN_DeltaE;
#endif 
}

/* Mirror of apply_kick(), but for the local copy of a remote cell (ghost) */
static void apply_kick_primexch(int particle, const struct Feedback_Kick *Kick)
{
#ifdef WINDS
  PrimExch[particle].MassFeed += Kick->DeltaMass;
#ifdef METALS
  PrimExch[particle].MetalsFeed += Kick->DeltaMetals;
#endif
  for(int k = 0; k < 3; k++)
    PrimExch[particle].MomentumFeed[k] += Kick->DeltaP[k];
#endif

#ifdef SUPERNOVAE
  PrimExch[particle].MassFeed += Kick->SN_DeltaMass;
#ifdef METALS
  PrimExch[particle].MetalsFeed += Kick->SN_DeltaMetals;
#endif
  for(int k = 0; k < 3; k++)
    PrimExch[particle].MomentumFeed[k] += Kick->SN_DeltaP[k];
#endif
}

#ifdef SUPERNOVAE
/* Compute corrent p and E scaling of SN explosion */
static void SN_compute(int ev, int h, double e, double a, double b, double NgbsDensity, double NgbsMetallicity, double *p, double *E)
{
  Mechanical_Feedback_Data *MechanicalFeedbackData = &MechanicalFeedbackEvents.MechanicalFeedbackData[ev + h];
  Mechanical_Feedback *MechanicalFeedback = &MechanicalFeedbackData->MechanicalFeedback;
  
  double E_SN = MechanicalFeedback->SN_EnergyInject;
  double m_ej = MechanicalFeedback->SN_MassLoss;

  /* Pure thermal injection: no ejecta mass, skip momentum prescription */
  if(m_ej <= 0.0)
    {
      *p = 0.0;
      *E = E_SN; 
      
      return;
    }
  
  

  double E_SNR = E_SN + e;
  double E51 = E_SN * All.cf_UnitEnergy_in_cgs / 1.0e51;
  
  double n_H  = 0.76 * NgbsDensity * All.cf_UnitDensity_in_cgs / PROTONMASS; 

  double Zsol = fmax(NgbsMetallicity / 0.0127, 0.01);

  /* Terminal momentum: Kim & Ostriker (2015) */
  double p_term = 3.0e5 /* Msun km/s */
  * pow(E51, 16.0 / 17.0)
  * pow(n_H, -2.0 / 17.0)
  * pow(Zsol, -0.14);
  
  p_term /= (All.cf_UnitMass_in_Msun * All.cf_UnitVelocity_in_cm_per_s / 1.0e5);
  
  /* Boost momentum */
  double fkin = 0.28;

  a *= m_ej;
  b *= sqrt(m_ej / (2 * fkin * E_SNR));
  
  double fboost = fmin((sqrt(b*b + a) - b) / a, p_term / (sqrt(2 * fkin * E_SNR * m_ej)));

  double p_SNR = fboost * sqrt(2 * fkin * E_SNR * m_ej);

  double sq_vstar = MechanicalFeedback->StarVelocity[0]*MechanicalFeedback->StarVelocity[0] 
  + MechanicalFeedback->StarVelocity[1]*MechanicalFeedback->StarVelocity[1] 
  + MechanicalFeedback->StarVelocity[2]*MechanicalFeedback->StarVelocity[2];

  double E_tot = 0.5 * m_ej * sq_vstar + E_SN;

  *p = p_SNR;
  *E = E_tot;
}
#endif 

#ifdef SUPERNOVAE
/* 
 * Sedov/cooling radius check: is the host cell able to resolve the
 * pressure-driven expansion of this SN?
 * Kim & Ostriker (2015)-> r_SN = 30 pc * (E_SN / 1e51 erg)^0.29 * (rho_h / (1.4*m_p) / cm^-3)^(-0.46) 
 * Unphysical radius check: is the host cell too large for us to 
 * couple the SN to its neighbours?
 * Mathew Smith (2026)-> r_unphysical = 1 kpc
 */
static int SN_feedback_radius(int i, int ev, int h)
{
  double r_host = get_cell_radius(i); 

  Mechanical_Feedback_Data *MechanicalFeedbackData = &MechanicalFeedbackEvents.MechanicalFeedbackData[ev + h];
  Mechanical_Feedback *MechanicalFeedback = &MechanicalFeedbackData->MechanicalFeedback;
  
  double E_SN = MechanicalFeedback->SN_EnergyInject;
 
  double E51 = E_SN * All.cf_UnitEnergy_in_cgs / 1.0e51;
 
  double rho_cgs = SphP[i].Density * All.cf_UnitDensity_in_cgs;
  double n_cgs = rho_cgs / (1.4 * PROTONMASS); 
 
  double r_SN_pc = 30.0 * pow(E51, 0.29) * pow(n_cgs, -0.46); /* parsec */
  double r_SN = (r_SN_pc * PARSEC) / All.cf_UnitLength_in_cm; /* code units */

  /* Thermal dump into host  */
  if(r_host < r_SN / 10.0)
    return SN_HOST;

  double r_unphysical = 1.0 * KILOPARSEC / All.cf_UnitLength_in_cm;
  
  /* Momentum kick to host */
  if(r_host > r_unphysical)
    return SN_HOST_MOMENTUM;
  
  /* Couple across Voronoi faces */ 
  return SN_MESH;
}
 
/* Host-only injection path: deposit this star's SN mass, momentum, and
 * energy budget directly into its host cell, with no
 * mesh-neighbour loop. Density/metallicity feeding into SN_compute() are
 * the host cell's own (unweighted) values, since there is nothing to
 * average over.
 */
static void SN_feedback_host(int i, int ev, int h, int mode)
{
  Mechanical_Feedback_Data *MechanicalFeedbackData = &MechanicalFeedbackEvents.MechanicalFeedbackData[ev + h];
  Mechanical_Feedback *MechanicalFeedback = &MechanicalFeedbackData->MechanicalFeedback;
 
  struct Feedback_Kick Kick = {0};
  Kick.CellIndex = i;

  if(mode == SN_HOST)
    {
      double m_ej = MechanicalFeedback->SN_MassLoss;
      double E_SN = MechanicalFeedback->SN_EnergyInject;

      double sq_vstar = MechanicalFeedback->StarVelocity[0] * MechanicalFeedback->StarVelocity[0]
      + MechanicalFeedback->StarVelocity[1] * MechanicalFeedback->StarVelocity[1]
      + MechanicalFeedback->StarVelocity[2] * MechanicalFeedback->StarVelocity[2];

      double E = 0.5 * m_ej * sq_vstar + E_SN;

      Kick.SN_DeltaMass = m_ej;
#ifdef METALS
      Kick.SN_DeltaMetals = MechanicalFeedback->SN_MetalsLoss;
#endif
      for(int k = 0; k < 3; k++)
        Kick.SN_DeltaP[k] = m_ej * MechanicalFeedback->StarVelocity[k];
 
      Kick.SN_DeltaE = E;
    }
  else if(mode == SN_HOST_MOMENTUM)
    {
  /* Star -> host cell-centre direction */
  double d[3], dd; 
                  
  d[0] = NEAREST_X(P[i].Pos[0] - MechanicalFeedback->StarPosition[0]);
  d[1] = NEAREST_Y(P[i].Pos[1] - MechanicalFeedback->StarPosition[1]);
  d[2] = NEAREST_Z(P[i].Pos[2] - MechanicalFeedback->StarPosition[2]);
              
  dd = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);

  double wbar[3] = {0};

  if(dd > 0.0)
    {
      for(int k = 0; k < 3; k++)
        {
          wbar[k] = d[k] / dd;
        }    
    }        
  
  /* Helpers for supernovae injection */
  double num, den, e = 0.0, a = 0.0, b = 0.0;
  double m_ej = MechanicalFeedback->SN_MassLoss;
              
  /* Ngbs properties */
  double NgbsDensity = 0.0, NgbsMetallicity = 0.0;

  /* Single target: 100% of the deposit */
  const double sq_wbar = 1.0; 
  const double sqrtsq_wbar = 1.0;

  double mj, vj[3];
  
  mj = P[i].Mass + SphP[i].StarMassFeed;
  
  for(int k = 0; k < 3; k++)
    vj[k] = (SphP[i].Momentum[k] + SphP[i].StarMomentumFeed[k]) / mj;

  double sq_vj = vj[0]*vj[0] + vj[1]*vj[1] + vj[2]*vj[2];
      
  double sq_vstar = MechanicalFeedback->StarVelocity[0]*MechanicalFeedback->StarVelocity[0]
  + MechanicalFeedback->StarVelocity[1]*MechanicalFeedback->StarVelocity[1]
  + MechanicalFeedback->StarVelocity[2]*MechanicalFeedback->StarVelocity[2];
      
  double cross = 2.0 * (vj[0]*MechanicalFeedback->StarVelocity[0]
  + vj[1]*MechanicalFeedback->StarVelocity[1]
  + vj[2]*MechanicalFeedback->StarVelocity[2]);
  
  num = 0.5 * mj * m_ej * sqrtsq_wbar * (sq_vj + sq_vstar - cross);
  den = mj + m_ej * sqrtsq_wbar;
 
  e = num / den;

  num = sq_wbar;

  a = num / den;

  num =  mj * ((vj[0] - MechanicalFeedback->StarVelocity[0]) * wbar[0]
  + (vj[1] - MechanicalFeedback->StarVelocity[1]) * wbar[1]
  + (vj[2] - MechanicalFeedback->StarVelocity[2]) * wbar[2]);

  b = num / den;
 
  /* Host-only density/metallicity: the host cell's own values, unweighted */
  NgbsDensity = SphP[i].Density;
#ifdef METALS
  NgbsMetallicity = (SphP[i].GasMetals + SphP[i].StarMetalsFeed) / mj;
#else
  NgbsMetallicity = 0.0;
#endif
 
  double p, E;
      
  SN_compute(ev, h, e, a, b, NgbsDensity, NgbsMetallicity, &p, &E);
 
  Kick.SN_DeltaMass = m_ej;
#ifdef METALS
  Kick.SN_DeltaMetals = MechanicalFeedback->SN_MetalsLoss;
#endif
  for(int k = 0; k < 3; k++)
    Kick.SN_DeltaP[k] = m_ej * MechanicalFeedback->StarVelocity[k] + p * wbar[k];
 
  Kick.SN_DeltaE = E;
  }
  else
    terminate("Should not have another mode!");
 
  apply_kick(i, &Kick);
}
#endif

/* star_feedback() -> main entry point */
void star_feedback(void)
{
  TIMER_START(CPU_STARS_FEEDBACK);

  #define MAX_FACES 128

  int ev, h, i, k, q, f;

  int n_export = 0;
  int max_export = 20 * MechanicalFeedbackEvents.NumEvents;
  
  int *ExportTask = (int *) mymalloc_movable(&ExportTask, "ExportTask", max_export * sizeof(int));
  struct Feedback_Kick *ExportBuf =
  (struct Feedback_Kick *) mymalloc_movable(&ExportBuf, "ExportBuf", max_export * sizeof(struct Feedback_Kick));

  /* Act on host cells */
  for(ev = 0; ev < MechanicalFeedbackEvents.NumEvents;)
    {
      i = MechanicalFeedbackEvents.MechanicalFeedbackData[ev].HostIndex;

      for(h = 0; h < SphP[i].Host; h++)
        {                    
          int flag_wind = 0, flag_sn = 0;
          int flag_wind_host = 0, flag_sn_host = 0;

          Mechanical_Feedback_Data *MechanicalFeedbackData = &MechanicalFeedbackEvents.MechanicalFeedbackData[ev + h];
          Mechanical_Feedback *MechanicalFeedback = &MechanicalFeedbackData->MechanicalFeedback;

#ifdef WINDS
          if(MechanicalFeedback->MassLoss)
            flag_wind = 1;
#endif

#ifdef SUPERNOVAE
          if(MechanicalFeedback->SN_MassLoss || MechanicalFeedback->SN_EnergyInject)
            flag_sn = 1;
#endif
        
#ifdef WINDS
          /* Wind radius */
          /* Skip the mesh-neighbour geometry entirely and dump Wind into the host cell */
          if(flag_wind)
            {
              //if(Wind_feedback_radius(i, ev, h))
              //  {
              //    Wind_feedback_host(i, ev, h);
              //    flag_wind_host = 1;
              //  }
            }
#endif

#ifdef SUPERNOVAE
          /* SN radius */
          /* Skip the mesh-neighbour geometry entirely and dump SN into the host cell */
          if(flag_sn)
            {
              int sn_mode = SN_feedback_radius(i, ev, h);

              if(sn_mode != SN_MESH)
                {
                  SN_feedback_host(i, ev, h, sn_mode);
                  flag_sn_host = 1;
                }
            }
#endif
          
          /* No feedback star */
          if(!flag_wind && !flag_sn)
            continue;
          
          /* Host deposition */
          if((!flag_wind || flag_wind_host) && (!flag_sn || flag_sn_host))
            continue;

          /* Mesh deposition */
          /* We will loop over the cell faces 4 times */
          int n_faces = 0;
          int dc_list[MAX_FACES];     

          /* GEOMETRY: passes 1 & 2 (Voronoi mesh) */
          /* Compute weights */
          double nplus[3], nminus[3], fplus[3], fminus[3];
          
          for(k = 0; k < 3; k++) 
            nplus[k] = nminus[k] = fplus[k] = fminus[k] = 0.0;
      
          /* Accumulate helper */
          double Splus[3], Sminus[3];
          
          for(k = 0; k < 3; k++)
            Splus[k] = Sminus[k] = 0.0;

          /* First pass */            
          q = SphP[i].first_connection;

          while(q >= 0)
            {
              if(q < 0 || q >= MaxNvc)
                {
                  char buf[1000];
                  sprintf(buf, "Strange connectivity q=%d Nvc=%d", q, MaxNvc);
                  terminate(buf);
                }
              
              int dp = DC[q].dp_index;
              int vf = DC[q].vf_index;
              int particle = Mesh.DP[dp].index;
          
              /* Cell has been removed */
              if(particle < 0)
                {
                  if(q == SphP[i].last_connection)
                    break;
                  
                  q = DC[q].next;
                    continue;
                }

              /* Face normal - from cell generator to cell generator */
              double n[3], nn; 

              n[0] = Mesh.DP[dp].x - P[i].Pos[0];
              n[1] = Mesh.DP[dp].y - P[i].Pos[1];
              n[2] = Mesh.DP[dp].z - P[i].Pos[2];
          
              nn = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

              /* Star-to-face direction */
              double d[3], dd; 
                  
              d[0] = Mesh.VF[vf].cx - MechanicalFeedback->StarPosition[0];
              d[1] = Mesh.VF[vf].cy - MechanicalFeedback->StarPosition[1];
              d[2] = Mesh.VF[vf].cz - MechanicalFeedback->StarPosition[2];
              
              dd = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);

              if(nn > 0.0 && dd > 0.0 && Mesh.VF[vf].area > 0.0)
                {
                  n[0] /= nn;  n[1] /= nn;  n[2] /= nn;
                  
                  d[0] /= dd;  d[1] /= dd;  d[2] /= dd;

                  double costheta = n[0]*d[0] + n[1]*d[1] + n[2]*d[2];
                  if(costheta < 0.0) costheta = 0.0;
            
                  double omega = 0.5 * (1 - 1 / sqrt(1 + Mesh.VF[vf].area * costheta / (M_PI * dd*dd)));

                  /* Star-to-cell direction */
                  double r[3], rr; 
                  
                  r[0] = Mesh.DP[dp].x - MechanicalFeedback->StarPosition[0];
                  r[1] = Mesh.DP[dp].y - MechanicalFeedback->StarPosition[1];
                  r[2] = Mesh.DP[dp].z - MechanicalFeedback->StarPosition[2];

                  rr = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

                  if(rr > 0.0)
                    {
                      r[0] /= rr;  r[1] /= rr;  r[2] /= rr;

                      for(k = 0; k < 3; k++)
                        {
                          nplus[k] = r[k] >= 0 ? r[k] : 0;
                          nminus[k] = r[k] < 0 ? r[k] : 0; 
                        
                          Splus[k] += omega * fabs(nplus[k]);
                          Sminus[k] += omega * fabs(nminus[k]);
                        }
                    }
                }
          
              if(n_faces >= MAX_FACES)
                terminate("star_feedback: MAX_FACES exceeded for cell %d\n", i);

              dc_list[n_faces++] = q;

              if(q == SphP[i].last_connection) 
                break;
          
              q = DC[q].next;
            }

          double weights[MAX_FACES][3];
          double wtot = 0.0; 
  
          /* Second pass */
          for(f = 0; f < n_faces; f++)
            {
              q = dc_list[f];
          
              int dp = DC[q].dp_index;
              int vf = DC[q].vf_index;

              double w[3] = {0.0, 0.0, 0.0};

             /* Face normal - from cell generator to cell generator */
              double n[3], nn; 

              n[0] = Mesh.DP[dp].x - P[i].Pos[0];
              n[1] = Mesh.DP[dp].y - P[i].Pos[1];
              n[2] = Mesh.DP[dp].z - P[i].Pos[2];
          
              nn = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

              /* Star-to-face direction */
              double d[3], dd; 
                  
              d[0] = Mesh.VF[vf].cx - MechanicalFeedback->StarPosition[0];
              d[1] = Mesh.VF[vf].cy - MechanicalFeedback->StarPosition[1];
              d[2] = Mesh.VF[vf].cz - MechanicalFeedback->StarPosition[2];
              
              dd = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);

              if(nn > 0.0 && dd > 0.0 && Mesh.VF[vf].area > 0.0)
                {
                  n[0] /= nn;  n[1] /= nn;  n[2] /= nn;
                  
                  d[0] /= dd;  d[1] /= dd;  d[2] /= dd;

                  double costheta = n[0]*d[0] + n[1]*d[1] + n[2]*d[2];
                  if(costheta < 0.0) costheta = 0.0;
            
                  double omega = 0.5 * (1 - 1 / sqrt(1 + Mesh.VF[vf].area * costheta / (M_PI * dd*dd)));

                  /* Star-to-cell direction */
                  double r[3], rr; 
                  
                  r[0] = Mesh.DP[dp].x - MechanicalFeedback->StarPosition[0];
                  r[1] = Mesh.DP[dp].y - MechanicalFeedback->StarPosition[1];
                  r[2] = Mesh.DP[dp].z - MechanicalFeedback->StarPosition[2];

                  rr = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

                  if(rr > 0.0)
                    {
                      r[0] /= rr;  r[1] /= rr;  r[2] /= rr;
                    
                      for(k = 0; k < 3; k++)
                        {
                          nplus[k] = r[k] >= 0 ? r[k] : 0;
                          nminus[k] = r[k] < 0 ? r[k] : 0; 
                          fplus[k] = (Splus[k] > 0.0) ? sqrt(0.5 * (1 + Sminus[k]*Sminus[k] / (Splus[k]*Splus[k]))) : 0;
                          fminus[k] =(Sminus[k] > 0.0) ? sqrt(0.5 * (1 + Splus[k]*Splus[k] / (Sminus[k]*Sminus[k]))) : 0;

                          w[k] = omega * (nplus[k] * fplus[k] + nminus[k] * fminus[k]); 
                        }
                    }                                         
                }
                
              for(k = 0; k < 3; k++)
                weights[f][k] = w[k];
                                    
              wtot += sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
            }

          if(wtot <= 0.0)
            terminate("STAR_FEEDBACK: invalid weight for host cell %d\n", i);
    
#ifdef SUPERNOVAE     
          double p, E;

          if(flag_sn && !flag_sn_host)
            {
              /* Helpers for supernovae injection */
              double num, den, e = 0.0, a = 0.0, b = 0.0;
              double m_ej = MechanicalFeedback->SN_MassLoss;
              
              /* Ngbs properties */
              int Ngbs = 0;
              double NgbsMass = 0.0, NgbsDensity = 0.0, NgbsMetallicity = 0.0;
   
              /* Third pass */ 
              for(f = 0; f < n_faces; f++)
                {
                  q = dc_list[f];

                  int dp = DC[q].dp_index;
                  int particle = Mesh.DP[dp].index;

                  if(particle >= NumGas && Mesh.DP[dp].task == ThisTask)
                    particle -= NumGas;

                  double wbar[3]; 
          
                  for(k = 0; k < 3; k++)
                    wbar[k] = weights[f][k] / wtot;

                  double sq_wbar = (wbar[0]*wbar[0] + wbar[1]*wbar[1] + wbar[2]*wbar[2]);
                  double sqrtsq_wbar = sqrt(sq_wbar);  

                  double mj, vj[3];
                  if(Mesh.DP[dp].task == ThisTask)
                    {
                      mj = P[particle].Mass + SphP[particle].StarMassFeed;
                      
                      for(k = 0; k < 3; k++)
                        vj[k] = (SphP[particle].Momentum[k] + SphP[particle].StarMomentumFeed[k]) / mj;   
                    }
                  else
                    {
                      mj = PrimExch[particle].Density * PrimExch[particle].Volume + PrimExch[particle].MassFeed;
                      
                      for(k = 0; k < 3; k++)
                        vj[k] = (PrimExch[particle].Density * PrimExch[particle].Volume * PrimExch[particle].VelGas[k]
                        + PrimExch[particle].MomentumFeed[k]) / mj;

                    }

                  double sq_vj = vj[0]*vj[0] + vj[1]*vj[1] + vj[2]*vj[2];
                  
                  double sq_vstar = MechanicalFeedback->StarVelocity[0]*MechanicalFeedback->StarVelocity[0]
                  + MechanicalFeedback->StarVelocity[1]*MechanicalFeedback->StarVelocity[1]
                  + MechanicalFeedback->StarVelocity[2]*MechanicalFeedback->StarVelocity[2];
                  
                  double cross = 2.0 * (vj[0]*MechanicalFeedback->StarVelocity[0]
                  + vj[1]*MechanicalFeedback->StarVelocity[1]
                  + vj[2]*MechanicalFeedback->StarVelocity[2]);
         
                  num = 0.5 * mj * m_ej * sqrtsq_wbar * (sq_vj + sq_vstar - cross);
                  den = mj + m_ej * sqrtsq_wbar;
          
                  e += num / den;
               
                  num = sq_wbar;
          
                  a += num / den;

                  num = mj * ((vj[0] - MechanicalFeedback->StarVelocity[0]) * wbar[0] 
                  + (vj[1] - MechanicalFeedback->StarVelocity[1]) * wbar[1]
                  + (vj[2] - MechanicalFeedback->StarVelocity[2]) * wbar[2]);

                  b += num / den;

                  if(Mesh.DP[dp].task == ThisTask)
                    {
                      Ngbs++;
                      NgbsMass += mj * sqrtsq_wbar;
                      NgbsDensity += mj / SphP[particle].Volume * sqrtsq_wbar;
#ifdef METALS
                      NgbsMetallicity += (SphP[particle].GasMetals + SphP[particle].StarMetalsFeed) / mj * sqrtsq_wbar;
#endif
                    }
                  else
                    {
                      Ngbs++;
                      NgbsMass += mj * sqrtsq_wbar;
                      NgbsDensity +=  mj / PrimExch[particle].Volume * sqrtsq_wbar;
#ifdef METALS
                      NgbsMetallicity += (PrimExch[particle].Density * PrimExch[particle].Volume * PrimExch[particle].Scalars[METALS_INDEX] 
                      + PrimExch[particle].MetalsFeed) / mj * sqrtsq_wbar;
#endif
                    }  
                }

              SN_compute(ev, h, e, a, b, NgbsDensity, NgbsMetallicity, &p, &E);
            }
#endif
      
          /* Fourth pass */  
          for(f = 0; f < n_faces; f++)
            {
              q = dc_list[f];

              int dp = DC[q].dp_index;
              int particle = Mesh.DP[dp].index;

              if(particle >= NumGas && Mesh.DP[dp].task == ThisTask)
                particle -= NumGas;

              double wbar[3]; 
          
              for(k = 0; k < 3; k++)
                wbar[k] = weights[f][k] / wtot;

              double sq_wbar = (wbar[0]*wbar[0] + wbar[1]*wbar[1] + wbar[2]*wbar[2]);
              double sqrtsq_wbar = sqrt(sq_wbar);
      
              struct Feedback_Kick Kick = {0};
              Kick.CellIndex = DC[q].index;
           
              /* Mesh ngbs feedback */
#ifdef WINDS
              if(flag_wind && !flag_wind_host)
                {
                  Kick.DeltaMass = MechanicalFeedback->MassLoss * sqrtsq_wbar;
#ifdef METALS
                  Kick.DeltaMetals = MechanicalFeedback->MetalsLoss * sqrtsq_wbar;
#endif    
                  for(k = 0; k < 3; k++)
                    Kick.DeltaP[k] = MechanicalFeedback->MassLoss * sqrtsq_wbar * MechanicalFeedback->StarVelocity[k] + MechanicalFeedback->WindMomentum * wbar[k];
          
                  double sq_vstar = MechanicalFeedback->StarVelocity[0]*MechanicalFeedback->StarVelocity[0] 
                  + MechanicalFeedback->StarVelocity[1]*MechanicalFeedback->StarVelocity[1] 
                  + MechanicalFeedback->StarVelocity[2]*MechanicalFeedback->StarVelocity[2];

                  double sq_vwind = MechanicalFeedback->WindMomentum / MechanicalFeedback->MassLoss
                  * MechanicalFeedback->WindMomentum / MechanicalFeedback->MassLoss; 

                  //double cross = 2.0 * (MechanicalFeedback->StarVelocity[0] * MechanicalFeedback->WindMomentum / MechanicalFeedback->MassLoss * wbar[0] 
                  //+ MechanicalFeedback->StarVelocity[1] * MechanicalFeedback->WindMomentum / MechanicalFeedback->MassLoss * wbar[1] 
                  //+ MechanicalFeedback->StarVelocity[2] * MechanicalFeedback->WindMomentum / MechanicalFeedback->MassLoss * wbar[2]);

                  Kick.DeltaE = 0.5 * MechanicalFeedback->MassLoss * (sq_vstar + sq_vwind) * sqrtsq_wbar;
                }   
#endif
 
#ifdef SUPERNOVAE
              if(flag_sn && !flag_sn_host)
                {
                  Kick.SN_DeltaMass = MechanicalFeedback->SN_MassLoss * sqrtsq_wbar;
#ifdef METALS
                  Kick.SN_DeltaMetals = MechanicalFeedback->SN_MetalsLoss * sqrtsq_wbar;
#endif 
                  for(k = 0; k < 3; k++)
                    Kick.SN_DeltaP[k] = MechanicalFeedback->SN_MassLoss * sqrtsq_wbar * MechanicalFeedback->StarVelocity[k] + p * wbar[k];

                  Kick.SN_DeltaE = E * sqrtsq_wbar;
                }
#endif

              if(DC[q].task == ThisTask)
                apply_kick(DC[q].index, &Kick);
              else
                {
                  /* Keep local copy current */
                  apply_kick_primexch(particle, &Kick);   

                  if(n_export >= max_export)
                    {
                      max_export *= 2;
 
                      ExportTask = (int *) myrealloc_movable(ExportTask, max_export * sizeof(int));
                      ExportBuf = (struct Feedback_Kick *) myrealloc_movable(ExportBuf, max_export * sizeof(struct Feedback_Kick));
                    }
                  
                  ExportTask[n_export] = DC[q].task;
                  ExportBuf[n_export] = Kick;
                  n_export++;
                }
            }
        } //for(int h = 0; h < SphP[i].Host; h++)
      
      /* Go to next host */
      ev += SphP[i].Host;
      /* All stars processed: release host slot */        
      //SphP[i].Host = 0;
    } //for(int ev = 0; ev < MechanicalFeedbackEvents.NumEvents;)

  /* MPI exchange of remote kick packets via MPI_Alltoallv */
  int *SendCount = mymalloc("FBSendCount", NTask * sizeof(int));
  int *RecvCount = mymalloc("FBRecvCount", NTask * sizeof(int));
  int *SendDisp = mymalloc("FBSendDisp",  NTask * sizeof(int));
  int *RecvDisp = mymalloc("FBRecvDisp",  NTask * sizeof(int));
 
  memset(SendCount, 0, NTask * sizeof(int));
  for(k = 0; k < n_export; k++)
    SendCount[ExportTask[k]]++;
 
  MPI_Alltoall(SendCount, 1, MPI_INT, RecvCount, 1, MPI_INT, MPI_COMM_WORLD);
 
  SendDisp[0] = RecvDisp[0] = 0;
  for(int t = 1; t < NTask; t++)
    {
      SendDisp[t] = SendDisp[t-1] + SendCount[t-1];
      RecvDisp[t] = RecvDisp[t-1] + RecvCount[t-1];
    }
  int n_recv = RecvDisp[NTask-1] + RecvCount[NTask-1];
 
  /* Sort ExportBuf into task-contiguous order for Alltoallv */
  struct Feedback_Kick *SortedExport = mymalloc("FBSortedExport", (n_export > 0 ? n_export : 1) * sizeof(struct Feedback_Kick));
  
  int *tmp_offset = mymalloc("FBTmpOffset", NTask * sizeof(int));
  memcpy(tmp_offset, SendDisp, NTask * sizeof(int));
  for(k = 0; k < n_export; k++)
    {
      int t = ExportTask[k];
      SortedExport[tmp_offset[t]++] = ExportBuf[k];
    }
  myfree(tmp_offset);
 
  struct Feedback_Kick *RecvBuf = mymalloc("RecvBuf",
  (n_recv > 0 ? n_recv : 1) * sizeof(struct Feedback_Kick));
 
  const int sz = (int)sizeof(struct Feedback_Kick);
 
  int *SendCountB = mymalloc("SendCountB", NTask * sizeof(int));
  int *RecvCountB = mymalloc("RecvCountB", NTask * sizeof(int));
  int *SendDispB  = mymalloc("SendDispB",  NTask * sizeof(int));
  int *RecvDispB  = mymalloc("RecvDispB",  NTask * sizeof(int));
 
  for(int t = 0; t < NTask; t++)
    {
      SendCountB[t] = SendCount[t] * sz;
      RecvCountB[t] = RecvCount[t] * sz;
      SendDispB[t]  = SendDisp[t]  * sz;
      RecvDispB[t]  = RecvDisp[t]  * sz;
    }
 
  MPI_Alltoallv(SortedExport, SendCountB, SendDispB, MPI_BYTE,
  RecvBuf, RecvCountB, RecvDispB, MPI_BYTE,
  MPI_COMM_WORLD);
 
  /* Free byte-count arrays in reverse allocation order (LIFO stack) */
  myfree(RecvDispB); myfree(SendDispB); myfree(RecvCountB); myfree(SendCountB);
 
  /* Apply received kicks to local cells */
  for(k = 0; k < n_recv; k++)
    apply_kick(RecvBuf[k].CellIndex, &RecvBuf[k]);
 
  /* Cleanup in reverse allocation order */
  myfree(RecvBuf); myfree(SortedExport);

  myfree(RecvDisp); myfree(SendDisp); myfree(RecvCount); myfree(SendCount);

  myfree(ExportBuf); myfree(ExportTask);
 
  TIMER_STOP(CPU_STARS_FEEDBACK);
}