#ifndef BH_PROTO_H
#define BH_PROTO_H


/* black hole functions */

/* Memory allocation */
void reallocate_memory_maxpartbhs(void);
void domain_resize_storage_bhs(int count_get_bh);

#ifdef BH_ACTIVE
/* Timesteps */
integertime bh_timestep(int p);
void bh_update_timesteps(void);
void bh_reconstruct_timebins(void);
void bh_update_list_of_active_particles(void);

/* Density loop */
void bh_density(void);
void bh_perform_end_of_step_physics(void);

void bh_kernel(double u, double hinv3, double hinv4, double *wk, double *dwk);
#endif

#ifdef BH_ACCRETION_ACTIVE
/* Accretion loops */
void bh_accretion(void);
void bh_swallow(void);
#endif

#ifdef BH_MERGER
/* BH-BH merging */
void bh_merger(void);

/* Gravitational force-softening length for the BH at local BhP[]/PPB() index idx.
 * Mirrors how the gravity routines look up per-particle softening (All.ForceSoftening[],
 * indexed by P[].SofteningType); exposed here (not static in bh_merger.c) so a future
 * BH-star capture routine can reuse the same lookup for a "close enough" radius test. */
double bh_softening_for_index(int idx);
#endif

/* Feedback loops */
#ifdef BH_THERMAL_FEEDBACK
void bh_feedback(void);
#endif

#ifdef BH_JET_PRESSURE_TRIGGER
void bh_jet_pressure_trigger(void);
#endif

#if defined(BH_JET_FEEDBACK) || defined(BH_JET_PRESSURE_TRIGGER)
/* kinetic energy [erg] the jet of BH i produces over its current step, and the length of that step [s] */
double bh_jet_step_energy_cgs(int i);
double bh_jet_step_time_cgs(int i);
#endif

#ifdef BH_JET_FEEDBACK
void bh_jet_prepare(void);   /* jet axis, energy reservoir, release decision (needs JetLaunchFlag) */
void bh_jet_density(void);   /* mass of the cone cells in each active BH kernel; marks cone cells for the timestep limiter */
void bh_jet_feedback(void);  /* distributes the released energy over the cone cells (buffers in SphP) */

/*! \brief Is the BH->cell vector d inside either lobe of the double cone about axis?
 *
 *  \param[in] d BH->cell separation (nearest image)
 *  \param[in] axis unit vector along the jet axis
 *  \param[in] cos_theta cosine of the half-opening angle
 *  \param[out] r distance |d|
 *
 *  \return 1 if inside the cone and r > 0, else 0
 */
static inline int bh_jet_in_cone(const double d[3], const MyDouble axis[3], double cos_theta, double *r)
{
  double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];

  *r = sqrt(r2);

  if(r2 <= 0)
    return 0;

  return fabs(d[0] * axis[0] + d[1] * axis[1] + d[2] * axis[2]) >= (*r) * cos_theta;
}
#endif

void blackhole_mark_cells_for_refinement(void);

#ifdef BLACKHOLES_FEEDBACK
void bh_jet_density(void);
void bh_ngb_feedback(void);
#endif

#endif /* #ifdef BLACKHOLES */

/* black hole seeding from on-the-fly FOF: see src/fof/fof_seeding.h */
