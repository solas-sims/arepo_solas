#ifndef DUST_H
#define DUST_H

/* Deliberately no #include of ../main/allvars.h here: this header only
 * defines rate-law constants and is shared by the pure, Arepo-independent
 * rate-law source files (dust_production_sn.c, dust_destruction_sn.c,
 * dust_growth.c, dust_sputtering.c) so they can be compiled standalone by
 * the tests/dust/ unit-test harness (no MPI, no mesh, no SphP/P globals). */

/*
 * Dust module, Phase 1 (single scalar dust mass per gas cell, see the
 * dust-module phaseplan doc). Physical processes implemented here:
 *
 *   - production:   SN-channel dust condensation only (Phase 1 scope decision:
 *                    the WINDS channel in this code blends massive-star and
 *                    AGB-type mass loss into one table with no way to isolate
 *                    true AGB ejecta, so wind-channel dust production is
 *                    deferred until that split exists)
 *   - growth:        Dwek (1998) grain-growth accretion timescale form
 *   - sputtering:     thermal sputtering, Tsai & Mathews (1995) timescale form
 *   - destruction:    SN shock destruction, McKee (1989) swept-mass form
 *
 * Rate-law constants below are set to match Li, Narayanan & Dave (2019,
 * MNRAS, arXiv:1906.09277) directly -- verified against the paper's primary
 * text and Table 1 (eq. 3, 5-8), not an AI-summarized secondary source. See
 * the verification note in DUST_PHASE1_IMPLEMENTATION.md for the full
 * paper-vs-code comparison this was derived from.
 */

/* No #ifdef DUST guard on the content below: this header (and dust_proto.h)
 * are only ever included from files that are themselves only compiled when
 * DUST is set -- the source files under src/dust (gated in the Makefile)
 * and other modules' "#ifdef DUST #include ..." sites. Guarding the
 * content here too would silently compile the pure rate-law source files
 * (which don't include allvars.h/arepoconfig.h, by design, for standalone
 * testability) into empty translation units in the real Arepo build. */

/* Every dust_* function takes an explicit species index so that Phase 2b
 * (DUST_NUMBER > 1) is a matter of looping this index, not a struct/API
 * rewrite. Phase 1 only ever calls these with species == 0. */
#define DUST_PHASE1_SPECIES 0

/* --- SN-channel production --------------------------------------------- */
/* Fraction of freshly-synthesized SN metal mass that condenses into dust.
 * Li+2019 Table 1: delta_i,dust^SNII = 0.15, uniform across elements (their
 * eq. 3 is per-element on elemental SNII ejecta mass; Phase 1's single
 * total-metals scalar applies this efficiency to total SN_MetalsLoss
 * instead -- a deliberate Phase 1 scope simplification, see "Scope
 * decisions" in DUST_PHASE1_IMPLEMENTATION.md). */
#define DUST_SN_CONDENSATION_EFFICIENCY 0.15

/* --- SN shock destruction (McKee 1989) ----------------------------------
 * Li+2019 eq. 8: tau_de = M_g / (epsilon * gamma * M_s), epsilon = 0.3
 * (Table 1). Phase 1 reuses this code's own Kim & Ostriker (2015)
 * shell-mass estimate (Msh -> dm_h in src/stars/star_feedback.c) in place
 * of the paper's independent M_s, but does apply the paper's epsilon=0.3
 * destruction efficiency to it -- see dust_destruction_sn() in
 * dust_destruction_sn.c and its call sites in star_feedback.c. */
#define DUST_SN_DESTRUCTION_EFFICIENCY 0.3

/* --- Grain growth (Dwek 1998) -------------------------------------------
 * Li+2019 eq. 5: tau_accr = tau_ref * (rho_ref/rho_g) * (T_ref/T_g) *
 * (Z_sun/Z_g), values from their text/Table 1: tau_ref=10 Myr,
 * rho_ref=100 H atoms cm^-3, T_ref=20 K, Z_sun=0.0134 (Asplund et al. 2009).
 * Both density and temperature ratios are LINEAR (not sqrt), and the
 * metallicity term is required -- growth speeds up at higher gas-phase
 * metallicity (more metals available to accrete). */
#define DUST_GROWTH_TAU_REF_GYR 0.01
#define DUST_GROWTH_REF_NH_CGS 1.0e2   /* cm^-3 */
#define DUST_GROWTH_REF_TEMP_K 20.0    /* K */
#define DUST_SOLAR_METALLICITY_MASSFRAC 0.0134 /* Asplund et al. 2009, as used by Li+2019 */

/* --- Thermal sputtering (Tsai & Mathews 1995) ---------------------------
 * Li+2019 eq. 6-7: tau_sp ~ (0.17 Gyr) * (rho_ref/rho_g) *
 * [(T_0/T_g)^omega + 1], omega=2.5, T_0=2e6 K, rho_ref=1e-27 g cm^-3 (a
 * TOTAL gas mass density, not a hydrogen number density); then
 * dM/dt|sputter = -M_dust / (tau_sp / 3). This code passes hydrogen number
 * density (n_H_cgs, matching dust_get_nH_and_temp() in dust_update.c) into
 * this rate law, so rho_ref is converted to the equivalent n_H via this
 * code's own n_H = HYDROGEN_MASSFRAC * rho / PROTONMASS convention:
 * n_H_ref = 0.76 * 1e-27 g cm^-3 / 1.67262178e-24 g = 4.543e-4 cm^-3. */
#define DUST_SPUTTER_TAU_REF_GYR 0.17
#define DUST_SPUTTER_REF_NH_CGS 4.543e-4  /* cm^-3, see conversion note above */
#define DUST_SPUTTER_REF_TEMP_K 2.0e6  /* K */
#define DUST_SPUTTER_OMEGA 2.5

/* Rate-ODE bisection tolerance / iteration cap, mirroring DoCooling() in
 * src/cooling/cooling.c */
#define DUST_ODE_TOLERANCE 1.0e-6
#define DUST_ODE_MAXITER 300

#endif /* DUST_H */
