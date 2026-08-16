#include <mpi.h>

#include "../main/allvars.h"
#include "../main/proto.h"
#include "dust.h"
#include "dust_proto.h"

#ifdef DUST

/*! \brief Register dust snapshot fields.
 *
 *  Called from init_io_fields() in src/io/io_fields.c under #ifdef DUST.
 *  No module currently has its own dedicated io file (all fields funnel
 *  through the shared io_fields.c) -- this is a deliberate small deviation
 *  from that convention, kept to a single hook call in the shared file so
 *  it stays low-risk.
 *
 *  Phase 1 keeps the snapshot footprint minimal: just the dust mass itself.
 *  D/Z and the gas-phase/dust-phase metal-budget split are left to
 *  post-processing (they're trivially derived from GasDustMass and the
 *  existing PassiveScalars metallicity output).
 */
void dust_init_io_fields(void)
{
  init_field(IO_DUST_MASS, "DUST", "DustMass", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].GasDustMass, 0,
             GAS_ONLY);
  init_units(IO_DUST_MASS, 0., -1., 0., 1., 0., All.UnitMass_in_g);
  init_snapshot_type(IO_DUST_MASS, SN_MINI);
}

#endif /* DUST */
