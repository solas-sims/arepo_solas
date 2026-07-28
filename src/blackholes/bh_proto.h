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

/* Feedback loops */
#ifdef BH_THERMAL_FEEDBACK
void bh_feedback(void);
#endif

#ifdef BH_JET_FEEDBACK
void bh_jet_density(void);
void bh_jet_feedback(void);
#endif

void blackhole_mark_cells_for_refinement(void);

#ifdef BLACKHOLES_FEEDBACK
void bh_jet_density(void);
void bh_ngb_feedback(void);
#endif

#endif /* #ifdef BLACKHOLES */

/* black hole seeding from on-the-fly FOF: see src/fof/fof_seeding.h */
