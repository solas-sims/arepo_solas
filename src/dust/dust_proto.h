#ifndef DUST_PROTO_H
#define DUST_PROTO_H

/* No #ifdef DUST guard here -- see the note at the top of dust.h. */

/* SN-channel production/destruction (event-driven, called from
 * src/stars/star_feedback.c per SN event). Pure functions: all masses must
 * be passed in the same (arbitrary but consistent) unit system -- unit
 * conversion is the caller's responsibility, so these stay testable
 * standalone. */
double dust_production_sn(int species, double SN_MetalsLoss);
double dust_destruction_sn(int species, double M_dust, double cell_gas_mass, double destroy_mass_threshold,
                            double weight);

/* Continuous rate terms (per active gas cell), used by dust_cell(). Return
 * dM/dt in [mass-unit of M_dust / Gyr]; density/temperature must be passed
 * in physical cgs units (n_H in cm^-3, T in K). */
double dust_growth_rate(int species, double M_dust, double M_metal, double n_H_cgs, double temp_K);
double dust_sputtering_rate(int species, double M_dust, double n_H_cgs, double temp_K);

/* Per-cell rate-ODE integration and driver, mirrors cool_cell()/cooling_only()
 * in src/cooling/cooling.c */
void dust_cell(int i);
void dust_processes(void);

/* Snapshot I/O registration, called from init_io_fields() in
 * src/io/io_fields.c */
void dust_init_io_fields(void);

#endif /* DUST_PROTO_H */
