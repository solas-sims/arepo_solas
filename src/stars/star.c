#include "../main/allvars.h"

#ifdef STARS

int NumStars;

#ifdef STAR_FEEDBACK_ACTIVE
struct TimeBinData TimeBinsStar;
#endif

#ifdef STAR_FEEDBACK_ACTIVE
Mechanical_Feedback_Events MechanicalFeedbackEvents;
#endif

Star_Particle_Data *SP;

#endif /* #ifdef STARS */
