#ifndef SIDM_H
#define SIDM_H

#include "../main/allvars.h"

#ifdef SIDM

void sidm_density(void);
void sidm_scatter(int timebin);

/*! Safety factor bounding per-step scattering probability (rate*dt) for
 *  the SIDM timestep criterion in get_timestep_gravity() (timestep.c).
 *  Hardcoded for v1 rather than a param.txt entry -- promote to a
 *  runtime parameter later if tuning turns out to be needed. */
#define SIDM_TIMESTEP_SAFETY_FACTOR 0.1

#endif /* #ifdef SIDM */

#endif /* #ifndef SIDM_H */
