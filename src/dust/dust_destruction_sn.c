#include "dust.h"
#include "dust_proto.h"

/*! \brief Dust mass destroyed in a gas cell by a single SN shock.
 *
 *  McKee (1989)-style swept-mass form: the destroyed fraction of a cell's
 *  dust is destroy_mass_threshold (the ISM mass processed by this SN's
 *  shock -- see the call sites in src/stars/star_feedback.c, which pass in
 *  this code's own Kim & Ostriker 2015 swept-mass estimate) as a fraction
 *  of the cell's gas mass, weighted by the same per-face/per-cell weight
 *  the SN event's metal deposit uses, capped at destroying all of the
 *  cell's dust.
 *
 *  Pure function -- M_dust, cell_gas_mass and destroy_mass_threshold must
 *  all be passed in the same consistent mass unit.
 *
 *  \param[in] species              dust species index (Phase 1: always DUST_PHASE1_SPECIES)
 *  \param[in] M_dust               current dust mass in the cell
 *  \param[in] cell_gas_mass        current gas mass of the cell
 *  \param[in] destroy_mass_threshold ISM mass processed by this SN's shock, in the caller's mass unit
 *  \param[in] weight               this SN event's fractional deposit onto this cell, in [0,1]
 *
 *  \return dust mass destroyed (>= 0, <= M_dust)
 */
double dust_destruction_sn(int species, double M_dust, double cell_gas_mass, double destroy_mass_threshold,
                            double weight)
{
  (void)species; /* Phase 1: single species; kept for the Phase 2b API */

  if(M_dust <= 0 || cell_gas_mass <= 0 || weight <= 0)
    return 0.0;

  double destroyed_fraction = weight * destroy_mass_threshold / cell_gas_mass;

  if(destroyed_fraction > 1.0)
    destroyed_fraction = 1.0;

  return destroyed_fraction * M_dust;
}
