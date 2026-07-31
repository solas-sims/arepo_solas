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
 * All rate-law constants below are the standard literature forms as commonly
 * cited in the dust-evolution literature (e.g. McKinnon et al. 2016; Li,
 * Narayanan & Dave 2019), NOT copied directly from a specific paper's table.
 * They have not been checked against Li et al. (2019) itself and should be
 * verified against that paper (or with collaborators) before treating any
 * output as science-grade -- see the "Rate-law coefficient caveat" note in
 * the Phase 1 implementation plan.
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
 * Placeholder value: SN dust condensation is generally found to be
 * inefficient once reverse-shock processing is accounted for. VERIFY. */
#define DUST_SN_CONDENSATION_EFFICIENCY 0.01

/* --- SN shock destruction (McKee 1989) ----------------------------------
 * Rather than an independent literature swept-mass constant, Phase 1 reuses
 * this code's own Kim & Ostriker (2015) shell-mass estimate (Msh -> dm_h in
 * src/stars/star_feedback.c) as the McKee (1989)-style "ISM mass cleared of
 * dust per SN": the dust contained in the swept mass is treated as fully
 * destroyed (sputtered) rather than surviving intact like the metal tracer
 * mass does. See dust_destruction_sn() in dust_destruction_sn.c and its
 * call sites in star_feedback.c for the exact accounting. */

/* --- Grain growth (Dwek 1998) -------------------------------------------
 * tau_growth = DUST_GROWTH_TAU_REF_GYR * (n_ref / n_H) * sqrt(T_ref / T)
 * dM/dt|growth = (M_dust / tau_growth) * (1 - M_dust / M_metal)
 * Timescale shortens (rate rises) with increasing T, reflecting faster
 * grain-gas collisions at higher thermal velocity (Hirashita 2000); the
 * reference density is chosen for cold, dense (molecular-cloud-like) gas,
 * so growth is suppressed mainly by the 1/n_H scaling in diffuse gas, not
 * by temperature. VERIFY against Li+2019 -- some implementations instead
 * restrict growth to cold/dense gas by phase cut rather than relying on
 * this T-scaling alone. */
#define DUST_GROWTH_TAU_REF_GYR 0.03
#define DUST_GROWTH_REF_NH_CGS 1.0e3   /* cm^-3 */
#define DUST_GROWTH_REF_TEMP_K 20.0    /* K */

/* --- Thermal sputtering (Tsai & Mathews 1995) ---------------------------
 * tau_sp = DUST_SPUTTER_TAU_REF_GYR * (n_ref / n_H) * [(T_ref/T)^omega + 1]
 * dM/dt|sputter = -M_dust / tau_sp
 * VERIFY against Li+2019 / Tsai & Mathews (1995) directly. */
#define DUST_SPUTTER_TAU_REF_GYR 0.17
#define DUST_SPUTTER_REF_NH_CGS 1.0    /* cm^-3 */
#define DUST_SPUTTER_REF_TEMP_K 2.0e6  /* K */
#define DUST_SPUTTER_OMEGA 2.5

/* Rate-ODE bisection tolerance / iteration cap, mirroring DoCooling() in
 * src/cooling/cooling.c */
#define DUST_ODE_TOLERANCE 1.0e-6
#define DUST_ODE_MAXITER 300

#endif /* DUST_H */
