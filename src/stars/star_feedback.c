#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

/* SN host-injection modes, returned by SN_feedback_radius() */
#define MESH 0 /* Couple across Voronoi faces */
#define HOST 1 /* Thermal dump into host  */

/* Percent thermal energy in the energy conserving phase (Sedov-Taylor) of an SN */
#define SN_F_THERMAL 0.72 

/* Kick packet sent to remote face-neighbor cells */
struct Feedback_Kick
{ 
  /* Cell index on the receiving task */
  int CellIndex; 

#ifdef WINDS
  MyDouble DeltaMass;
#if GRACKLE_CHEMISTRY >= 1
  MyDouble DeltaChem[GRACKLE_SPECIES_NUMBER];
#endif
#ifdef METALS
  MyDouble DeltaMetals;
#endif
  MyDouble DeltaP[3];
  MyDouble DeltaE;
#endif

#ifdef SUPERNOVAE
  MyDouble SN_DeltaMass;
#if GRACKLE_CHEMISTRY >= 1
  MyDouble SN_DeltaChem[GRACKLE_SPECIES_NUMBER];
#endif
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
#if GRACKLE_CHEMISTRY >= 1
  for(int s = 0; s < GRACKLE_SPECIES_NUMBER; s++)
    SphP[j].StarChemFeed[s] += Kick->DeltaChem[s];
#endif
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
#if GRACKLE_CHEMISTRY >= 1
  for(int s = 0; s < GRACKLE_SPECIES_NUMBER; s++)
    SphP[j].StarChemFeed[s] += Kick->SN_DeltaChem[s];
#endif
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
/*
 * Solve the energy-conserving SN momentum quadratic
 *
 * alpha p^2 + beta p + gamma = 0 , gamma = dU_SN,th - E_SN <= 0
 *
 * and return the physical (positive) root
 *
 * p_SN = (-beta + sqrt(beta^2 - 4 alpha gamma)) / (2 alpha)
 */
static void SN_compute(int ev, int h, double alpha, double beta, double gamma, double NgbsDensity, double NgbsMetallicity, double *p)
{
  if(alpha <= 0.0)
    {
      *p = 0.0;
      return;
    }

  double disc = beta * beta - 4.0 * alpha * gamma;
  
  if(disc < 0.0)
    disc = 0.0;

  double p_SN = (-beta + sqrt(disc)) / (2.0 * alpha);
  
  if(p_SN < 0.0)
    p_SN = 0.0;

  /* Terminal momentum: T. Kimm et al. (2015) and Matthew Smith (2026) */
  Mechanical_Feedback_Data *MechanicalFeedbackData = &MechanicalFeedbackEvents.MechanicalFeedbackData[ev + h];
  Mechanical_Feedback *MechanicalFeedback = &MechanicalFeedbackData->MechanicalFeedback;

  double E_SN  = MechanicalFeedback->SN_EnergyInject;
  double E51 = E_SN * All.cf_UnitEnergy_in_cgs / SN_ENERGY;
  
  double n_H = HYDROGEN_MASSFRAC * NgbsDensity * All.cf_UnitDensity_in_cgs / PROTONMASS;
  
  double Zsol = fmax(NgbsMetallicity / 0.0127, 0.01);

 
  double p_term = 3.0e5 /* Msun km/s */
                * pow(E51, 16.0 / 17.0)
                * pow(n_H, -2.0 / 17.0)
                * pow(Zsol, -0.14);
  
  p_term /= (All.cf_UnitMass_in_Msun * All.cf_UnitVelocity_in_cm_per_s / 1.0e5);

  if(p_SN > p_term)
    p_SN = p_term;

  *p = p_SN;
}
#endif

#ifdef WINDS
/*
 * TODO
 */
static int Wind_feedback_radius(int i, int ev, int h)
{
  return MESH;
}

/*
 * Host-only injection path: deposit this star's wind mass, metals,
 * momentum, and energy budget directly into its host cell, with no
 * mesh-neighbour loop.
 */
static void Wind_feedback_host(int i, int ev, int h, int mode)
{
  terminate("Wind_feedback_host: host-mode winds not yet implemented (cell %d, mode %d)\n", i, mode);
}
#endif

#ifdef SUPERNOVAE
/* 
 * Sedov/cooling radius check: is the host cell able to resolve the
 * pressure-driven expansion of this SN?
 * Kim & Ostriker (2015)-> r_SN = 22.6 pc * (E_SN / SN_ENERGY erg)^0.29 * (rho_h / (1.4*m_p) / cm^-3)^(-0.42) 
 */
static int SN_feedback_radius(int i, int ev, int h)
{
  double r_host = get_cell_radius(i); 

  Mechanical_Feedback_Data *MechanicalFeedbackData = &MechanicalFeedbackEvents.MechanicalFeedbackData[ev + h];
  Mechanical_Feedback *MechanicalFeedback = &MechanicalFeedbackData->MechanicalFeedback;
  
  double E_SN = MechanicalFeedback->SN_EnergyInject;
 
  double E51 = E_SN * All.cf_UnitEnergy_in_cgs / SN_ENERGY;
  
  double rho = (P[i].Mass + SphP[i].StarMassFeed) / SphP[i].Volume;
  double rho_cgs = rho * All.cf_UnitDensity_in_cgs;

  /* Based on Kim & Ostriker */
  double n_cgs = rho_cgs / (1.4 * PROTONMASS); 
 
  double r_SN_pc = 22.6 * pow(E51, 0.29) * pow(n_cgs, -0.42); /* parsec */
  double r_SN = (r_SN_pc * PARSEC) / All.cf_UnitLength_in_cm; /* code units */

  /* Thermal dump into host  */
  if(r_host < r_SN / 10.0)
    return HOST;
  
  /* Couple across Voronoi faces */ 
  return MESH;
}
 
/* 
 * Host-only injection path: deposit this star's SN mass, momentum, and
 * energy budget directly into its host cell, with no
 * mesh-neighbour loop. Density/metallicity feeding into SN_compute() are
 * the host cell's own (unweighted) values, since there is nothing to average over.
 */
static void SN_feedback_host(int i, int ev, int h, int mode)
{
  Mechanical_Feedback_Data *MechanicalFeedbackData = &MechanicalFeedbackEvents.MechanicalFeedbackData[ev + h];
  Mechanical_Feedback *MechanicalFeedback = &MechanicalFeedbackData->MechanicalFeedback;

  int k;

  double m_ej = MechanicalFeedback->SN_MassLoss;

  double sq_vstar = MechanicalFeedback->StarVelocity[0]*MechanicalFeedback->StarVelocity[0]
  + MechanicalFeedback->StarVelocity[1]*MechanicalFeedback->StarVelocity[1]
  + MechanicalFeedback->StarVelocity[2]*MechanicalFeedback->StarVelocity[2];

  double E;

  if(mode == HOST)
    {
      /* Thermal + advected-mass dump: no directed momentum into the host */
      E = 0.5 * m_ej * sq_vstar + MechanicalFeedback->SN_EnergyInject;
    }
  else
    terminate("SN_feedback_host: unknown mode %d for host cell %d\n", mode, i);

  struct Feedback_Kick Kick = {0};
  Kick.CellIndex = i;

  Kick.SN_DeltaMass = m_ej;
#if GRACKLE_CHEMISTRY >= 1
  Kick.SN_DeltaChem[CHEM_HII] = MechanicalFeedback->SN_HLoss;
  Kick.SN_DeltaChem[CHEM_HeIII] = MechanicalFeedback->SN_HeLoss;
#endif
#ifdef METALS
  Kick.SN_DeltaMetals = MechanicalFeedback->SN_MetalsLoss;
#endif
  /* Advected mass only: ejecta carries the star's own velocity */
  for(k = 0; k < 3; k++)
    Kick.SN_DeltaP[k] = m_ej * MechanicalFeedback->StarVelocity[k];

  Kick.SN_DeltaE = E;

  apply_kick(i, &Kick);
}
#endif

/* star_feedback() -> main entry point */
void star_feedback(void)
{
  TIMER_START(CPU_STARS_FEEDBACK);

  #define MAX_FACES 128

  int ev, h, i, k, q, f;
  double xtmp, ytmp, ztmp;

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
              int wind_mode = Wind_feedback_radius(i, ev, h);
              
              if(wind_mode != MESH)
                {
                  Wind_feedback_host(i, ev, h, wind_mode);
                  flag_wind_host = 1;
                }
            }
#endif

#ifdef SUPERNOVAE
          /* SN radius */
          /* Skip the mesh-neighbour geometry entirely and dump SN into the host cell */
          if(flag_sn)
            {
              int sn_mode = SN_feedback_radius(i, ev, h);

              if(sn_mode != MESH)
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

          /* Treat periodic boundaries */
          double xstar[3];

          xstar[0] = P[i].Pos[0] - NEAREST_X(P[i].Pos[0] - MechanicalFeedback->StarPosition[0]);
          xstar[1] = P[i].Pos[1] - NEAREST_Y(P[i].Pos[1] - MechanicalFeedback->StarPosition[1]);
          xstar[2] = P[i].Pos[2] - NEAREST_Z(P[i].Pos[2] - MechanicalFeedback->StarPosition[2]);

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
                  
              d[0] = Mesh.VF[vf].cx - xstar[0];
              d[1] = Mesh.VF[vf].cy - xstar[1];
              d[2] = Mesh.VF[vf].cz - xstar[2];
              
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
                  
                  r[0] = Mesh.DP[dp].x - xstar[0];
                  r[1] = Mesh.DP[dp].y - xstar[1];
                  r[2] = Mesh.DP[dp].z - xstar[2];

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
                  
              d[0] = Mesh.VF[vf].cx - xstar[0];
              d[1] = Mesh.VF[vf].cy - xstar[1];
              d[2] = Mesh.VF[vf].cz - xstar[2];
              
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
                  
                  r[0] = Mesh.DP[dp].x - xstar[0];
                  r[1] = Mesh.DP[dp].y - xstar[1];
                  r[2] = Mesh.DP[dp].z - xstar[2];

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
          /* Directed momentum magnitude and thermal budget for this event */
          double p_SN = 0.0;
          double dU_SN_th = 0.0;
          
          /* Ejecta mass */
          double m_ej = MechanicalFeedback->SN_MassLoss;

          /* Ejecta internal energy; cold ejecta default */
          double U_ej = 0.0;

          /* Approximate SN shell mass */
          double E_SN = MechanicalFeedback->SN_EnergyInject;
 
          double E51 = E_SN * All.cf_UnitEnergy_in_cgs / SN_ENERGY;

          double rho = (P[i].Mass + SphP[i].StarMassFeed) / SphP[i].Volume;
          double rho_cgs = rho * All.cf_UnitDensity_in_cgs;

          double n_cgs = rho_cgs / (1.4 * PROTONMASS);
          
          /* Shell mass from Kim & Ostriker (2015) */
          double Msh = (1680.0 / All.cf_UnitMass_in_Msun) * pow(E51, 0.87) * pow(n_cgs, -0.26);   

          /* Host swept mass and internal energy */
          double m_h = P[i].Mass + SphP[i].StarMassFeed;
          
          double shell_sweep_frac = fmin(All.SN_HostShellSweepFrac, 0.9);
          double dm_h = shell_sweep_frac * fmin(Msh, m_h);
          dm_h = fmin(dm_h, m_h - 0.1 * P[i].Mass); 
          dm_h = fmax(dm_h, 0.0);

#if GRACKLE_CHEMISTRY >= 1
          double dmChem_h[GRACKLE_SPECIES_NUMBER];
            for(int s = 0; s < GRACKLE_SPECIES_NUMBER; s++)
              dmChem_h[s] = dm_h * (SphP[i].GrackleSpeciesConserved(GRACKLE_SPECIES_INDEX + s) + SphP[i].StarChemFeed[s]) / m_h;
#endif

#ifdef METALS
          double dmZ_h = dm_h * (SphP[i].GasMetals + SphP[i].StarMetalsFeed) / m_h;
#endif

          double p_h[3], v_h[3];
          for(k = 0; k < 3; k++)
            {
              p_h[k] = (SphP[i].Momentum[k] + SphP[i].StarMomentumFeed[k]);
              v_h[k] = (SphP[i].Momentum[k] + SphP[i].StarMomentumFeed[k]) / m_h;
            } 

          double K_h = (p_h[0]*p_h[0] + p_h[1]*p_h[1] + p_h[2]*p_h[2]) / (2 * m_h);
          double U_h = SphP[i].Energy + SphP[i].StarEnergyFeed - K_h;
          
          if(U_h < 0.0)
            U_h = 0.0;    
          
          double dU_h = (dm_h / m_h) * U_h; 

          /* Total advected (ejecta + swept host) mass */
          double m_feed = m_ej + dm_h; 

          /* Per-face state carried from the coefficient pass to deposition */
          double SN_ptilde[MAX_FACES][3]; /* pre-kick momentum m_j v_j + advected */
          double SN_mfj[MAX_FACES];       /* final neighbour mass m_j + |w|m_dep   */
          double SN_KEj[MAX_FACES];       /* initial neighbour kinetic energy      */
          double SN_Dj[MAX_FACES];        /* kinetic energy dissipated in-place    */

          if(flag_sn && !flag_sn_host)
            {
              /* Quadratic coefficients alpha p^2 + beta p + gamma = 0 */
              double alpha = 0.0, beta = 0.0;

              /* Ngbs properties (for the optional terminal-momentum cap) */
              double NgbsDensity = 0.0, NgbsMetallicity = 0.0;

              double sq_vstar = MechanicalFeedback->StarVelocity[0]*MechanicalFeedback->StarVelocity[0]
              + MechanicalFeedback->StarVelocity[1]*MechanicalFeedback->StarVelocity[1]
              + MechanicalFeedback->StarVelocity[2]*MechanicalFeedback->StarVelocity[2];

              double sq_vh = v_h[0]*v_h[0] + v_h[1]*v_h[1] + v_h[2]*v_h[2];

              /* Coefficient pass */
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
                      mj = PrimExch[particle].Mass + PrimExch[particle].MassFeed;

                      for(k = 0; k < 3; k++)
                        vj[k] = (PrimExch[particle].Momentum[k] + PrimExch[particle].MomentumFeed[k]) / mj;
                    }

                  double sq_vj = vj[0]*vj[0] + vj[1]*vj[1] + vj[2]*vj[2];

                  /* Pre-kick momentum: neighbour + advected ejecta (at v_star) + advected swept host (at v_host) */
                  double ptilde[3];
                  for(k = 0; k < 3; k++)
                    ptilde[k] = mj * vj[k] + (m_ej * MechanicalFeedback->StarVelocity[k] + dm_h * v_h[k]) * sqrtsq_wbar;

                  double mfj = mj + m_feed * sqrtsq_wbar;
                  double sq_ptilde = ptilde[0]*ptilde[0] + ptilde[1]*ptilde[1] + ptilde[2]*ptilde[2];
                  double bw_dot_pt = wbar[0]*ptilde[0] + wbar[1]*ptilde[1] + wbar[2]*ptilde[2];

                  /* alpha = sum |wbar|^2 / (2 m_f) ; beta = sum wbar.ptilde / m_f */
                  alpha += sq_wbar / (2.0 * mfj);
                  beta += bw_dot_pt / mfj;

                  /* Kinetic energy dissipated in the inelastic merge of the
                   * (neighbour, ejecta, host) streams; thermalised in place. */
                  double ke_streams = 0.5 * mj * sq_vj + 0.5 * m_ej * sqrtsq_wbar * sq_vstar + 0.5 * dm_h * sqrtsq_wbar * sq_vh;

                  double D = ke_streams - 0.5 * sq_ptilde / mfj;
                  
                  if(D < 0.0)
                    D = 0.0;

                  for(k = 0; k < 3; k++)
                    SN_ptilde[f][k] = ptilde[k];
                  
                  SN_mfj[f] = mfj;
                  SN_KEj[f] = 0.5 * mj * sq_vj;
                  SN_Dj[f] = D;

                  if(Mesh.DP[dp].task == ThisTask)
                    {
                      NgbsDensity += mj / SphP[particle].Volume * sqrtsq_wbar;
#ifdef METALS
                      NgbsMetallicity += (SphP[particle].GasMetals + SphP[particle].StarMetalsFeed) / mj * sqrtsq_wbar;
#endif
                    }
                  else
                    {
                      NgbsDensity += mj / PrimExch[particle].Volume * sqrtsq_wbar;
#ifdef METALS
                      NgbsMetallicity += (PrimExch[particle].Metals + PrimExch[particle].MetalsFeed) / mj * sqrtsq_wbar;
#endif
                    }
                }

              /* gamma = dU_SN,th - E_SN, then solve the quadratic for p_SN */
              dU_SN_th = SN_F_THERMAL * E_SN;
              double gamma = dU_SN_th - E_SN;

              SN_compute(ev, h, alpha, beta, gamma, NgbsDensity, NgbsMetallicity, &p_SN);
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
#if GRACKLE_CHEMISTRY >= 1
                  Kick.DeltaChem[CHEM_HI] = MechanicalFeedback->HLoss * sqrtsq_wbar;
                  Kick.DeltaChem[CHEM_HeI] = MechanicalFeedback->HeLoss * sqrtsq_wbar;
#endif
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

                  Kick.DeltaE = 0.5 * MechanicalFeedback->MassLoss * (sq_vstar + sq_vwind) * sqrtsq_wbar;
                }   
#endif
 
#ifdef SUPERNOVAE
              if(flag_sn && !flag_sn_host)
                {
                  /* Advected mass = ejecta + swept host, shared by |wbar| */
                  Kick.SN_DeltaMass = m_feed * sqrtsq_wbar;
#if GRACKLE_CHEMISTRY >= 1
                  for(int s = 0; s < GRACKLE_SPECIES_NUMBER; s++)
                    Kick.SN_DeltaChem[s] = dmChem_h[s] * sqrtsq_wbar;
                  Kick.SN_DeltaChem[CHEM_HII] += MechanicalFeedback->SN_HLoss * sqrtsq_wbar;
                  Kick.SN_DeltaChem[CHEM_HeIII] += MechanicalFeedback->SN_HeLoss * sqrtsq_wbar;
#endif 
#ifdef METALS
                  Kick.SN_DeltaMetals = (MechanicalFeedback->SN_MetalsLoss + dmZ_h) * sqrtsq_wbar ;
#endif
                  /* Momentum = advected (ejecta @ v_star + host @ v_host)
                   * + directed energy-conserving kick p_SN * wbar.          */
                  for(k = 0; k < 3; k++)
                    Kick.SN_DeltaP[k] = (m_ej * MechanicalFeedback->StarVelocity[k] + dm_h * v_h[k]) * sqrtsq_wbar
                    + p_SN * wbar[k];

                  /* Total energy = kinetic change + internal-energy gain.
                   * StarEnergyFeed tracks *total* energy, so we form the
                   * exact final-minus-initial cell energy here:
                   *   DKE = |ptilde + p_SN wbar|^2 / (2 m_f) - KE_j
                   *   dU  = |wbar|(U_ej + dU_h + dU_SN,th) + D_j            */
                  double pf[3];
                  for(k = 0; k < 3; k++)
                    pf[k] = SN_ptilde[f][k] + p_SN * wbar[k];

                  double sq_pf = pf[0]*pf[0] + pf[1]*pf[1] + pf[2]*pf[2];

                  double DKE = sq_pf / (2.0 * SN_mfj[f]) - SN_KEj[f];
                  double dU = (U_ej + dU_h + dU_SN_th) * sqrtsq_wbar + SN_Dj[f];

                  Kick.SN_DeltaE = DKE + dU;
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

#ifdef SUPERNOVAE
          if(flag_sn && !flag_sn_host && dm_h > 0.0)
            {
              double sq_vh = v_h[0]*v_h[0] + v_h[1]*v_h[1] + v_h[2]*v_h[2];

              struct Feedback_Kick Kick_h = {0};
              Kick_h.CellIndex = i;

              Kick_h.SN_DeltaMass = -dm_h;
#if GRACKLE_CHEMISTRY >= 1
              for(int s = 0; s < GRACKLE_SPECIES_NUMBER; s++)
                Kick_h.SN_DeltaChem[s] = -dmChem_h[s];
#endif
#ifdef METALS
              Kick_h.SN_DeltaMetals = -dmZ_h;
#endif
              for(k = 0; k < 3; k++)
                Kick_h.SN_DeltaP[k] = -dm_h * v_h[k];

              Kick_h.SN_DeltaE = -(0.5 * dm_h * sq_vh + dU_h);

              apply_kick(i, &Kick_h);
            }
#endif
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