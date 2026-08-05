#include <math.h>

#include "../main/allvars.h"
#include "../main/proto.h"


static inline int rt_node_inverted(int no)
{
  return RtNgb_Nodes[no].rt_range_min[0] >= RtNgb_Nodes[no].rt_range_max[0];
}

static inline int ray_box_intersect(const double *ray_pos, const double *ray_dir,
                                    const MyNgbTreeFloat *rmin, const MyNgbTreeFloat *rmax,
                                    double *t_enter, double *t_exit)
{
  double xtmp, ytmp, ztmp;

  double center[3] = {0.5 * (rmin[0] + rmax[0]), 0.5 * (rmin[1] + rmax[1]), 0.5 * (rmin[2] + rmax[2])};

  /* Minimum-image displacement of box centre relative to ray origin */
  double d[3];

  d[0] = NEAREST_X(center[0] - ray_pos[0]);
  d[1] = NEAREST_Y(center[1] - ray_pos[1]);
  d[2] = NEAREST_Z(center[2] - ray_pos[2]);

  double halfextent[3] = {0.5 * (rmax[0] - rmin[0]), 0.5 * (rmax[1] - rmin[1]), 0.5 * (rmax[2] - rmin[2])};

  double halfdomain[3] = {boxHalf_X, boxHalf_Y, boxHalf_Z};

  double tmin = -MAX_REAL_NUMBER, tmax = MAX_REAL_NUMBER;

  for(int i = 0; i < 3; i++)
    {
      /* Box too wide along this axis for a single minimum-image shift to be trustworthy */
      if(halfextent[i] >= halfdomain[i])
        { 
          /* Conservative hit */ 
          *t_enter = 0.0;
          *t_exit = MAX_REAL_NUMBER;
          return 1;
        }

      double shifted_min = ray_pos[i] + d[i] - halfextent[i];
      double shifted_max = ray_pos[i] + d[i] + halfextent[i];

      /* Ray parallel to this slab */
      if(fabs(ray_dir[i]) < 1e-12)
        {
          if(ray_pos[i] < shifted_min || ray_pos[i] > shifted_max)
            return 0;
        }
      else
        {
          double inv_dir = 1.0 / ray_dir[i];
          double t1 = (shifted_min - ray_pos[i]) * inv_dir;
          double t2 = (shifted_max - ray_pos[i]) * inv_dir;

          if(t1 > t2)
            {
              double tmp = t1;
              t1 = t2;
              t2 = tmp;
            }

          tmin = t1 > tmin ? t1 : tmin;
          tmax = t2 < tmax ? t2 : tmax;

          if(tmin > tmax)
            return 0;
        }
    }

  if(tmax < 0)
    return 0;

  *t_enter = fmax(tmin, 0.0);
  *t_exit  = tmax;

  return 1;
}

static inline int ray_sphere_intersect(const double *ray_pos, const double *ray_dir, 
                                       const double *center, const double r2,
                                       double *t_enter, double *t_exit)
{
  double xtmp, ytmp, ztmp;

  /* Minimum-image displacement of sphere centre relative to ray origin */
  double d[3];

  d[0] = NEAREST_X(center[0] - ray_pos[0]);
  d[1] = NEAREST_Y(center[1] - ray_pos[1]);
  d[2] = NEAREST_Z(center[2] - ray_pos[2]);

  double dist2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
  double t_closest = d[0] * ray_dir[0] + d[1] * ray_dir[1] + d[2] * ray_dir[2];
  int origin_inside = (dist2 < r2);

  /* Sphere is behind the ray, and ray starts outside it */
  if(!origin_inside && t_closest <= 0.0)
    return 0;

  /* Squared impact parameter */
  double b2 = dist2 - t_closest * t_closest;

  /* Ray misses the sphere entirely */
  if(!origin_inside && b2 >= r2)
    return 0;

  double dt = sqrt(r2 - b2);

  *t_enter = origin_inside ? 0.0 : t_closest - dt;
  *t_exit  = t_closest + dt;

  return 1;
}

static inline int ray_absorb(RayPacket *ray, double Dtau_E[WAVEBANDS], double Dtau_N[WAVEBANDS], WavebandData absorbed[WAVEBANDS],
                             double dN_H2, double *lw_line)
{
  for(int w = 0; w < WAVEBANDS; w++)
    {
      absorbed[w].Energy = absorbed[w].Photons = 0.0;

      /* Deactivate band if it has fallen below the dead-fraction threshold */
      if(ray->Radiated[w].Energy < RAD_TRUNC_FRAC * ray->Radiated_Init[w].Energy && 
        ray->Radiated[w].Photons < RAD_TRUNC_FRAC * ray->Radiated_Init[w].Photons)
        ray->active_bands &= (uint8_t)(~(1u << w));

      if(!(ray->active_bands & (1u << w)))
        continue;

      /* Separate treatment of LW*/
      if(w == LYMAN_WERNER)
        continue;
   
      double absorbed_energy = ray->Radiated[w].Energy * (1.0 - exp(-Dtau_E[w]));
      double absorbed_photons = ray->Radiated[w].Photons * (1.0 - exp(-Dtau_N[w]));

      absorbed[w].Energy += absorbed_energy;
      ray->Radiated[w].Energy -= absorbed_energy;
      
      absorbed[w].Photons += absorbed_photons;
      ray->Radiated[w].Photons -= absorbed_photons; 
    }

  /* LW band: H2 line self shielding + dust absorption */
  if(ray->active_bands & (1u << LYMAN_WERNER))
    {
      double Dtau_line = h2shield_dtau(ray->N_H2, dN_H2);
      double Dtau_dust_E = Dtau_E[LYMAN_WERNER];
      double Dtau_dust_N = Dtau_N[LYMAN_WERNER];

      double tot_E = Dtau_line + Dtau_dust_E; 
      double tot_N = Dtau_line + Dtau_dust_N;
  
      double absorbed_energy = ray->Radiated[LYMAN_WERNER].Energy * (1.0 - exp(-tot_E));
      double absorbed_photons = ray->Radiated[LYMAN_WERNER].Photons * (1.0 - exp(-tot_N));

      absorbed[LYMAN_WERNER].Energy += absorbed_energy;
      ray->Radiated[LYMAN_WERNER].Energy -= absorbed_energy;
      
      absorbed[LYMAN_WERNER].Photons += absorbed_photons;
      ray->Radiated[LYMAN_WERNER].Photons -= absorbed_photons;

      if(tot_E > 0) 
        lw_line[0] = Dtau_line / tot_E;
      if(tot_N > 0) 
        lw_line[1] = Dtau_line / tot_N;
    } 
  
  ray->N_H2 += dN_H2;
    
  return ray->active_bands != 0;
}

/* 
Tree indices are organized as follows:

[0 ... Ngb_MaxPart-1] -> real particles

[Ngb_MaxPart ... Ngb_MaxPart+Ngb_MaxNodes-1] -> internal nodes

    └── [Ngb_MaxPart ... Tree_FirstNonTopLevelNode-1] -> top-level nodes (replicated everywhere)
 
                             └──  [Tree_FirstNonTopLevelNode ... Ngb_MaxPart+Ngb_MaxNodes-1] -> local branch nodes

[Ngb_MaxPart+Ngb_MaxNodes ... Ngb_MaxPart+Ngb_MaxNodes+NTopleaves-1] -> pseudo-particles
*/

void raytrace_treewalk(RayPacket *ray, RayWorkStack *work, RayExportBuffer *export_buf)
{
  double xtmp, ytmp, ztmp;
  
  /* Local stack for ordering within this domain */
  StackEntry stack[RAY_STACK_SIZE];
  int stack_top = 0;

  /* Push entry point */
  if(ray->target_node < 0 )
    /* Root */
    stack[stack_top++] = (StackEntry){0.0, MAX_REAL_NUMBER, Ngb_MaxPart}; 
  else
    {
      memcpy(stack, ray->pending, ray->n_pending * sizeof(StackEntry));
      stack_top = ray->n_pending;
      ray->n_pending = 0;
      /* Push the target node on top - it goes first */
      stack[stack_top++] = (StackEntry){ray->t, ray->t_exit, ray->target_node};
    }
  
  while(stack_top > 0)
    {
      StackEntry cur = stack[--stack_top];
      int no = cur.node;

      /* ---- Cell ---- */
      if(no < Ngb_MaxPart)
        {     
          if(P[no].Type != 0 || P[no].Mass == 0 || P[no].ID == 0)
            continue;

          double length = cur.t_exit - cur.t_enter;
              
          double Dtau_E[WAVEBANDS];
          for(int w = 0; w < WAVEBANDS; w++)
            Dtau_E[w] = SphP[no].DtauOverLength_E[w] * length;

          double Dtau_N[WAVEBANDS];
          for(int w = 0; w < WAVEBANDS; w++)
            Dtau_N[w] = SphP[no].DtauOverLength_N[w] * length;
          
          WavebandData absorbed[WAVEBANDS];

          /* Line Dissociation */
          double dN_H2 = SphP[no].GrackleSpeciesConserved(GRACKLE_H2I) / SphP[no].Volume * length;
          /* Percent LW absorption that goes into H2 line dissociation */
          double lw_line[2] = {0.0, 0.0};

          /* Proccess ray */
          int still_alive = ray_absorb(ray, Dtau_E, Dtau_N, absorbed, dN_H2, lw_line);

          /* Reradiation in the IR (Boosts momentum) */
          double Dtau_IR = dtau_IR(no, length);
          
          /* Deposit absorbed energy into cells, one band at a time */
          double dK_total = 0.0;
          for(int w = 0; w < WAVEBANDS; w++)
            {
              double dp;
              
              /* No IR reradiation */
              if(w == IONIZING_HI || w == IONIZING_HeI || w == IONIZING_HeII)
                dp = absorbed[w].Energy / (CLIGHT / All.cf_UnitVelocity_in_cm_per_s) / All.cf_atime;
              /* IR reradiation (dust only) */
              else if(w == LYMAN_WERNER)
                dp = absorbed[w].Energy * (1.0 + (1.0 - lw_line[0]) * Dtau_IR * ReradiatedFraction[w]) / (CLIGHT / All.cf_UnitVelocity_in_cm_per_s) / All.cf_atime;
              /* IR reradiation (full) */
              else
                dp = absorbed[w].Energy * (1.0 + Dtau_IR * ReradiatedFraction[w]) / (CLIGHT / All.cf_UnitVelocity_in_cm_per_s) / All.cf_atime;        
            
              double dp_vec[3] = {dp * ray->dir[0], dp * ray->dir[1], dp * ray->dir[2]};

              /* Partially updated state */
              double mj, pj[3];

              mj = P[no].Mass + SphP[no].StarMassFeed;
              for(int k = 0; k < 3; k++)
                pj[k] = SphP[no].Momentum[k] + SphP[no].StarMomentumFeed[k];

              double sq_momentum = dp_vec[0]*dp_vec[0] + dp_vec[1]*dp_vec[1] + dp_vec[2]*dp_vec[2];

              double cross = 2 * (pj[0] * dp_vec[0] + pj[1] * dp_vec[1] + pj[2] * dp_vec[2]);

              double dK = (sq_momentum + cross) / (2 * mj);

              SphP[no].StarMomentumFeed[0] += dp_vec[0];
              SphP[no].StarMomentumFeed[1] += dp_vec[1];
              SphP[no].StarMomentumFeed[2] += dp_vec[2];

              if(w == LYMAN_WERNER)
                SphP[no].Absorbed[w].Energy += (1.0 - lw_line[0]) * (absorbed[w].Energy - dK);
              else
                SphP[no].Absorbed[w].Energy += absorbed[w].Energy - dK;

              dK_total += dK;
            }

          SphP[no].StarEnergyFeed += dK_total;
          All.StarFeedbackLocal[2] += dK_total;
          
          /* Deposit absorbed photons into cells, one band at a time */
          /* Dissociating Photons */
           SphP[no].Absorbed[LYMAN_WERNER].Photons += lw_line[1] * absorbed[LYMAN_WERNER].Photons; 

          /* Ionizing Photons */
          SphP[no].Absorbed[IONIZING_HI].Photons += absorbed[IONIZING_HI].Photons;
          SphP[no].Absorbed[IONIZING_HeI].Photons += absorbed[IONIZING_HeI].Photons;
          SphP[no].Absorbed[IONIZING_HeII].Photons += absorbed[IONIZING_HeII].Photons;
         
          ray->t = cur.t_exit;

          if(ray->t >= ray->t_maximum) 
            {
              ray->is_paused = 1; 
              return;
            }
          
          /* All bands are exhausted */
          if(!still_alive) 
            return;      
        }
      /* ---- Internal node ---- */
      else if(no < Ngb_MaxPart + Ngb_MaxNodes)
        {
          struct NgbNODE *nop = &Ngb_Nodes[no];
          struct RtNgbNODE *rt_nop = &RtNgb_Nodes[no];

          /* Node geometry */
          /* Node center */
          double cx = 0.5 * (rt_nop->rt_range_max[0] + rt_nop->rt_range_min[0]);
          double cy = 0.5 * (rt_nop->rt_range_max[1] + rt_nop->rt_range_min[1]);
          double cz = 0.5 * (rt_nop->rt_range_max[2] + rt_nop->rt_range_min[2]);
                            
          /* Node extent */ 
          double dx = rt_nop->rt_range_max[0] - rt_nop->rt_range_min[0];
          double dy = rt_nop->rt_range_max[1] - rt_nop->rt_range_min[1];
          double dz = rt_nop->rt_range_max[2] - rt_nop->rt_range_min[2];

          /* Node silhouette */
          double ax = fabs(ray->dir[0]), ay = fabs(ray->dir[1]), az = fabs(ray->dir[2]);
          double A_proj = dy*dz*ax + dx*dz*ay + dx*dy*az;   

          /* Node center to ray origin distance */
          double ddx = NEAREST_X(cx - ray->pos[0]);
          double ddy = NEAREST_Y(cy - ray->pos[1]);
          double ddz = NEAREST_Z(cz - ray->pos[2]);

          double dist2 = ddx*ddx + ddy*ddy + ddz*ddz;

#ifdef RAD_OPENING_ANGLE
          /* -- Barnes-Hut opening criterion -- */
          /* This should only trigger for non-top level nodes */ 
          if(no >= Ngb_FirstNonTopLevelNode)
            {     
              /* Node is far enough - treat as single slab */
              if(dist2 > 0 && A_proj / dist2 < All.RadOpeningAngle * All.RadOpeningAngle)
                {
                  /* Node aspect */
                  double lo = fmin(dx, fmin(dy, dz));
                  double hi = fmax(dx, fmax(dy, dz));
                  double aspect = (lo > 0.0) ? hi / lo : MAX_REAL_NUMBER;
                  
                  /* This should only trigger for not too elongated nodes */
                  if(aspect < All.NodeAspectRatio)
                    {
                      double length = cur.t_exit - cur.t_enter;
              
                      double Dtau_E[WAVEBANDS];
                      for(int w = 0; w < WAVEBANDS; w++)
                        /* Volume-weighted mean kappa * density */
                        Dtau_E[w] = RtNgb_Nodes[no].DtauOverLength_E[w] * length;
                        
                      double Dtau_N[WAVEBANDS];
                      for(int w = 0; w < WAVEBANDS; w++)
                        /* Volume-weighted mean kappa * density */
                        Dtau_N[w] = RtNgb_Nodes[no].DtauOverLength_N[w] * length;

                      WavebandData absorbed[WAVEBANDS];

                      double dN_H2 = RtNgb_Nodes[no].dN_H2_OverLength * length;
                      double lw_line[2] = {0.0, 0.0};

                      int still_alive = ray_absorb(ray, Dtau_E, Dtau_N, absorbed, dN_H2, lw_line);

                      /* Accumulate for later distribution to children */
                      for(int w = 0; w < WAVEBANDS; w++)
                        {
                          if(w == LYMAN_WERNER)
                            {
                              RtNgb_Nodes[no].Absorbed[w].Energy += (1.0 - lw_line[0]) * absorbed[w].Energy;
                              RtNgb_Nodes[no].Absorbed[w].Photons += lw_line[1] * absorbed[w].Photons;
                            }
                          else
                            {
                              RtNgb_Nodes[no].Absorbed[w].Energy += absorbed[w].Energy;
                              RtNgb_Nodes[no].Absorbed[w].Photons += absorbed[w].Photons;
                            }
                        }

                      ray->t = cur.t_exit;

                      if(ray->t >= ray->t_maximum) 
                        {
                          ray->is_paused = 1; 
                          return;
                        }

                      if(!still_alive) 
                        return;
                  
                      /* Don't open this node */
                      continue;  
                    } /* Else: Elongated node */
                } /* Else: Node too big/close */
            } /* Else: Top level node */
#endif      
          /* -- Adaptive splitting criterion -- */
          /* This should not trigger for the root node */
          if(no > Ngb_MaxPart)
            {
              /* This should not trigger for maximum resolution rays */
              if(ray->nside < NSIDE_MAX)
                {
                  /* Use number of actual children for adaptive f */
                  /* At least 1 ray per child */
                  double f_eff = All.RaySplitFactor * (double)RtNgb_Nodes[no].Nchildren; 

                  /* Ray is too coarse - push split children to split_buf, consume parent */
                  /* Criterion: Omega_node = A_proj / dist2 < f * Omega_ray = f * 4pi/(12*nside^2) */
                  if(dist2 > 0.0 && A_proj / dist2 < f_eff * 4.0 * M_PI / (12 * (double)(ray->nside * ray->nside)))
                    {
                      /* Pack pending */
                      if(stack_top > RAY_PENDING_SIZE)
                        terminate("Too many pending entries to split!");
                      
                      ray->n_pending = stack_top;
                      memcpy(ray->pending, stack, stack_top * sizeof(StackEntry));

                      ray->t = cur.t_enter;
                      ray->t_exit = cur.t_exit;
                  
                      /* Store for re-entry */
                      ray->target_node = no;
                      
                      RayPacket children[4];
                      split_ray(ray, children);
                    
                      for(int k = 0; k < 4; k++)                      
                        append_ray(work, &children[k]);
                   
                      /* Parent ray is consumed */
                      return;   
                                  
                    } /* Else: ray fine enough for node */
                } /* Else: ray at max resolution, open node */
            } /* Else: Root node */
         
          /* Open node and enumerate children -> sort by t_enter, push */
          StackEntry children[8];
          int nchildren = 0;

          int child = nop->u.d.nextnode;
          while(child != nop->u.d.sibling && child >= 0)
            {
              double t_enter, t_exit;
              int hit = 0;

              /* Cell */
              if(child < Ngb_MaxPart) 
                {             
                  double r = get_cell_radius(child);
                  double r2 = r * r;
                      
                  double cpos[3] = {P[child].Pos[0], P[child].Pos[1], P[child].Pos[2]};
                  
                  hit = ray_sphere_intersect(ray->pos, ray->dir, cpos, r2, &t_enter, &t_exit);            
                }
              /* Internal node */  
              else if(child < Ngb_MaxPart + Ngb_MaxNodes) 
                {   
                  if(!rt_node_inverted(child))
                    hit = ray_box_intersect(ray->pos, ray->dir, RtNgb_Nodes[child].rt_range_min, RtNgb_Nodes[child].rt_range_max, &t_enter, &t_exit);
                }
              /* Pseudo-particle: remote domain */
              else 
                {
                  int pseudo_idx = child - (Ngb_MaxPart + Ngb_MaxNodes);
                  int top_node = Ngb_DomainNodeIndex[pseudo_idx];

                  if(!rt_node_inverted(top_node))
                    hit = ray_box_intersect(ray->pos, ray->dir, RtNgb_Nodes[top_node].rt_range_min, RtNgb_Nodes[top_node].rt_range_max, &t_enter, &t_exit);
                }

              if(hit)
                {
                  if(t_enter < ray->t_maximum)
                    {
                      /* Limit traversal distance */
                      t_exit = fmin(t_exit, ray->t_maximum);  
                  
                      if(nchildren >= 8)
                        terminate("Too many children!");

                      children[nchildren++] = (StackEntry){t_enter, t_exit, child};
                    }
                }

              /* Advance to next direct child via sibling */
              if(child < Ngb_MaxPart)
                child = Ngb_Nextnode[child];
              else if(child < Ngb_MaxPart + Ngb_MaxNodes)
                child = Ngb_Nodes[child].u.d.sibling;
              else
                child = Ngb_Nextnode[child - Ngb_MaxNodes];
            }

          /* Sort ascending by t_enter */
          for(int i = 1; i < nchildren; i++)
            {
              StackEntry key = children[i];
              int j = i - 1;
              while(j >= 0 && children[j].t_enter > key.t_enter)
                {
                  children[j+1] = children[j];
                  j--;
                }
              children[j+1] = key;
            }

          /* Push in reverse so smallest t_enter is popped first */
          for(int i = nchildren - 1; i >= 0; i--)
            {
              if(stack_top >= RAY_STACK_SIZE)
              terminate("Ray stack overflow!");

              stack[stack_top++] = children[i];
            }
        }
    
      /* ---- Pseudo-particle: remote domain ---- */  
      else
        {
          int task = DomainTask[no - (Ngb_MaxPart + Ngb_MaxNodes)];
          int remote_node = Ngb_DomainNodeIndex[no - (Ngb_MaxPart + Ngb_MaxNodes)];

          /* Pack pending */
          if(stack_top > RAY_PENDING_SIZE) 
            terminate("Too many pending entries to export!");

          ray->n_pending = stack_top;
          memcpy(ray->pending, stack, stack_top * sizeof(StackEntry));

          ray->t = cur.t_enter;
          ray->t_exit = cur.t_exit;
          /* Store for re-entry */
          ray->target_node = remote_node;  

          /* Add to export buffer */
          append_export(export_buf, ray, task);

        /* This rank is done with this ray */
        return;
        }
       
      if(stack_top == 0)
        if(ray->t < ray->t_maximum) 
          {
            ray->t = ray->t_maximum;
            ray->is_paused = 1; 
            return;
          }
    }
}