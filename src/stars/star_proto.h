#ifndef STAR_PROTO_H
#define STAR_PROTO_H


/* Star functions */

/* Memory allocation */
void reallocate_memory_maxpartstars(void);
void domain_resize_storage_stars(int count_get_star);

#ifdef STAR_PARTICLES
/* IMF */
double imf_kroupa(double m); 
double imf_chabrier(double m); 
double imf_salpeter(double m);
double imf(double m); 
double m_times_imf(double m); 
void build_imf_cdf(void);
double sample_imf(double u);
#endif

#if defined(STAR_PARTICLES) && STAR_PARTICLES < 2
void setup_mass_bins(void);
void sample_star_particle(double m, int *bins);
#ifdef STAR_FEEDBACK_ACTIVE
Star_Feedback star_particle_feedback(int index, double dt, double z, double a);
#endif
#endif

#if STAR_PARTICLES == 0 
void setup_imf_integrals(void);
#endif

#ifdef INDIVIDUAL_STAR_BY_STAR_FORMATION
void individual_starbystar_formation(void);
void sf_starbystar(void);
void sf_massdrain(void);
double gaussian_weight(double r, double h);
#endif

#ifdef STAR_FEEDBACK_ACTIVE
/* Feedback tables interpolation */
void load_star_tables(const char *filename);
void free_stellar_tables(void);
Star_Feedback star_feedback_compute(double dt, double z_val, double m_val, double a);
Star_Feedback units_for_feedback(Star_Feedback star);

double IntegralTrapezoidal(double a, double b, int N, double (*f)(double));

/* Wrap stars module */
void star_in(void);
void star_exit(void);

/* Feedback setup */
void feedback_init(struct Mechanical_Feedback_Events *MFEvents);
void feedback_allocate(struct Mechanical_Feedback_Events *MFEvents, int MaxEvents);
void feedback_reallocate(struct Mechanical_Feedback_Events *MFEvents, int NewMaxEvents);
void feedback_free(struct Mechanical_Feedback_Events *MFEvents);

/* Timesteps */
integertime star_timestep(int p);
void star_update_timesteps(void);
void star_reconstruct_timebins(void);
void star_update_list_of_active_particles(void);

/* Prep and density */
void star_prep(void);
void star_density(void);
#endif

#ifdef STAR_RADIATION_ACTIVE
/* Radiation */
#include "../stars/star_radiation.h"

void update_dtau(void);
double dtau_IR(int i, double length);

void init_h2shield(void);
double h2shield_dtau(double N_H2, double dN_H2);

void append_ray(RayWorkStack *w, const RayPacket *ray);
void split_ray(const RayPacket *parent, RayPacket children[4]);

void append_export(RayExportBuffer *buf, const RayPacket *ray, int task);

void star_radiation(void);
void raytrace_treewalk(RayPacket *ray, RayWorkStack *work, RayExportBuffer *export_buf);
#endif

#if defined(WINDS) || defined(SUPERNOVAE)
void star_feedback(void);
#endif

#ifdef STAR_FEEDBACK_ACTIVE
void star_perform_end_of_step_physics(void);
#endif

#endif