/*
 * Standalone unit test for the dust module's pure rate-law functions
 * (src/dust/dust_production_sn.c, dust_destruction_sn.c, dust_growth.c,
 * dust_sputtering.c).
 *
 * Deliberately decoupled from the rest of Arepo: no MPI, no mesh, no
 * SphP/P globals -- these rate-law files only depend on dust.h/dust_proto.h
 * and libm, so they can be compiled and exercised standalone. See
 * tests/dust/build.sh.
 *
 * This checks conservation/limit properties of the rate laws themselves,
 * not their numerical integration in dust_cell() (src/dust/dust_update.c),
 * which requires the full Arepo build (see the Phase 1 implementation
 * plan's "What I can verify myself vs. what needs you" section).
 */

#include <math.h>
#include <stdio.h>

#include "../../src/dust/dust.h"
#include "../../src/dust/dust_proto.h"

static int g_failures = 0;

#define CHECK(cond, msg)                                    \
  do                                                         \
    {                                                        \
      if(!(cond))                                            \
        {                                                    \
          printf("FAIL: %s (%s:%d)\n", msg, __FILE__, __LINE__); \
          g_failures++;                                      \
        }                                                    \
      else                                                   \
        {                                                    \
          printf("PASS: %s\n", msg);                         \
        }                                                    \
    }                                                        \
  while(0)

int main(void)
{
  /* --- dust_production_sn ------------------------------------------- */
  CHECK(dust_production_sn(DUST_PHASE1_SPECIES, 0.0) == 0.0, "production_sn: zero metal loss -> zero dust");
  CHECK(dust_production_sn(DUST_PHASE1_SPECIES, -1.0) == 0.0, "production_sn: negative input -> zero dust (guarded)");

  {
    double m = 100.0;
    double d = dust_production_sn(DUST_PHASE1_SPECIES, m);
    CHECK(d > 0.0 && d < m, "production_sn: positive input -> 0 < dust < metal loss");
    CHECK(fabs(d - DUST_SN_CONDENSATION_EFFICIENCY * m) < 1e-12, "production_sn: matches condensation efficiency exactly");
  }

  /* --- dust_destruction_sn -------------------------------------------- */
  CHECK(dust_destruction_sn(DUST_PHASE1_SPECIES, 0.0, 1e4, 1e3, 1.0) == 0.0, "destruction_sn: zero dust -> zero destroyed");
  CHECK(dust_destruction_sn(DUST_PHASE1_SPECIES, 10.0, 1e4, 1e3, 0.0) == 0.0, "destruction_sn: zero weight -> zero destroyed");

  {
    /* threshold == cell_gas_mass, weight=1 -> full destruction */
    double destroyed = dust_destruction_sn(DUST_PHASE1_SPECIES, 10.0, 1000.0, 1000.0, 1.0);
    CHECK(fabs(destroyed - 10.0) < 1e-9, "destruction_sn: threshold==cell mass, weight=1 -> destroys all dust");
  }
  {
    /* threshold > cell_gas_mass must still cap at destroying all dust, not overshoot */
    double destroyed = dust_destruction_sn(DUST_PHASE1_SPECIES, 10.0, 1000.0, 5000.0, 1.0);
    CHECK(fabs(destroyed - 10.0) < 1e-9, "destruction_sn: oversized threshold caps at full destruction, no overshoot");
  }
  {
    /* linear scaling in threshold, below the cap */
    double d1 = dust_destruction_sn(DUST_PHASE1_SPECIES, 10.0, 1000.0, 100.0, 1.0);
    double d2 = dust_destruction_sn(DUST_PHASE1_SPECIES, 10.0, 1000.0, 200.0, 1.0);
    CHECK(fabs(d2 - 2.0 * d1) < 1e-9, "destruction_sn: linear in destroy_mass_threshold below the cap");
  }

  /* --- dust_growth_rate ------------------------------------------------ */
  CHECK(dust_growth_rate(DUST_PHASE1_SPECIES, 0.0, 100.0, 1e3, 20.0) == 0.0, "growth_rate: zero dust mass -> zero rate");
  CHECK(dust_growth_rate(DUST_PHASE1_SPECIES, 50.0, 50.0, 1e3, 20.0) == 0.0,
        "growth_rate: M_dust==M_metal (fully condensed) -> zero rate (saturated)");

  {
    double rate = dust_growth_rate(DUST_PHASE1_SPECIES, 10.0, 100.0, 1e3, 20.0);
    CHECK(rate > 0.0, "growth_rate: 0 < M_dust < M_metal -> positive growth rate");
  }
  {
    /* Higher density -> shorter growth timescale -> higher rate */
    double rate_lo = dust_growth_rate(DUST_PHASE1_SPECIES, 10.0, 100.0, 1e2, 20.0);
    double rate_hi = dust_growth_rate(DUST_PHASE1_SPECIES, 10.0, 100.0, 1e4, 20.0);
    CHECK(rate_hi > rate_lo, "growth_rate: rate increases with density");
  }
  {
    /* Dwek (1998)/Hirashita (2000) accretion timescale scales as 1/sqrt(T)
     * (collision rate with thermal velocity), so tau shortens -- and the
     * rate rises -- as T increases. */
    double rate_cold = dust_growth_rate(DUST_PHASE1_SPECIES, 10.0, 100.0, 1e3, 20.0);
    double rate_hot   = dust_growth_rate(DUST_PHASE1_SPECIES, 10.0, 100.0, 1e3, 2000.0);
    CHECK(rate_hot > rate_cold, "growth_rate: rate increases with temperature (faster thermal collisions)");
  }

  /* --- dust_sputtering_rate --------------------------------------------- */
  CHECK(dust_sputtering_rate(DUST_PHASE1_SPECIES, 0.0, 1.0, 1e6) == 0.0, "sputtering_rate: zero dust mass -> zero rate");

  {
    double rate = dust_sputtering_rate(DUST_PHASE1_SPECIES, 10.0, 1.0, 1e6);
    CHECK(rate < 0.0, "sputtering_rate: positive dust mass -> negative (destructive) rate");
    CHECK(fabs(rate) < 10.0 * 1e6, "sputtering_rate: magnitude is finite and not absurd");
  }
  {
    /* Lower density -> longer sputtering timescale -> smaller |rate| */
    double rate_lo_n = dust_sputtering_rate(DUST_PHASE1_SPECIES, 10.0, 0.01, 1e6);
    double rate_hi_n = dust_sputtering_rate(DUST_PHASE1_SPECIES, 10.0, 100.0, 1e6);
    CHECK(fabs(rate_hi_n) > fabs(rate_lo_n), "sputtering_rate: |rate| increases with density");
  }

  /* --- combined: zero-rate no-op (Phase 0 exit-gate spirit) -------------- */
  {
    /* With M_dust already at the metal budget and negligible density/temp
     * driving growth (M_dust == M_metal keeps growth at exactly zero via
     * saturation), only sputtering can move M -- confirm it does not
     * spontaneously produce dust out of nothing at M_dust == 0. */
    double g = dust_growth_rate(DUST_PHASE1_SPECIES, 0.0, 100.0, 1e3, 20.0);
    double s = dust_sputtering_rate(DUST_PHASE1_SPECIES, 0.0, 1.0, 1e6);
    CHECK(g == 0.0 && s == 0.0, "combined: M_dust=0 is a fixed point of growth+sputtering (no spontaneous creation)");
  }

  printf("\n%s (%d failure%s)\n", g_failures == 0 ? "ALL TESTS PASSED" : "TESTS FAILED", g_failures, g_failures == 1 ? "" : "s");

  return g_failures == 0 ? 0 : 1;
}
