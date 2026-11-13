#include "dust.h"
#include "dust_proto.h"

/*! \brief Dust mass condensed from a single SN event's ejected metal mass.
 *
 *  Phase 1: SN-channel only (see dust.h). Pure function -- SN_MetalsLoss and
 *  the return value share whatever consistent mass unit the caller uses.
 *
 *  \param[in] species     dust species index (Phase 1: always DUST_PHASE1_SPECIES)
 *  \param[in] SN_MetalsLoss metal mass ejected by this SN event
 *
 *  \return dust mass condensed from this event
 */
double dust_production_sn(int species, double SN_MetalsLoss)
{
  (void)species; /* Phase 1: single species; kept for the Phase 2b API */

  if(SN_MetalsLoss <= 0)
    return 0.0;

  return DUST_SN_CONDENSATION_EFFICIENCY * SN_MetalsLoss;
}
