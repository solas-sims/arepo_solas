# Self-interacting dark matter (SIDM)

Elastic, single-state SIDM with a constant cross-section per unit mass, following the Monte Carlo scattering
scheme of Vogelsberger, Zavala & Loeb (2012) as described in Valdarnini (2023, arXiv:2309.10374). Enable it with
the `SIDM` Config flag; `Config_SIDM.sh` in the repository root is a worked example (isolated halo).

This page is the user-facing summary. For design decisions, debugging history and status, read
[`SIDM.md`](../../SIDM.md); for a function-by-function reference, [`IMPLEMENTATION.md`](../../IMPLEMENTATION.md).

**Status:** implemented and running; momentum and energy conservation are checked live, and a first isolated
NFW halo shows the expected cusp-to-core behaviour. Quantitative validation against published results is **not
done**, so treat results as unvalidated.

## Config flags

| Flag | Default | Effect |
|---|---|---|
| `SIDM` | off | Enables the module: a dedicated dark-matter-only neighbour tree, SIDM density estimate, a timestep criterion and the scattering step. |
| `SIDM_KICK_POSITION_CORRECTION` | off | Corrects positions for the velocity change at a scattering and wakes a passively kicked inactive neighbour. Unvalidated; compare against a run without it. |
| `SIDM_LOG_COLLISIONS` | off | Diagnostic only: writes detected collisions to `sidm_collisions_<task>.txt` instead of applying them. Never use in production runs. |

## Parameters

| Parameter | Type | Meaning |
|---|---|---|
| `SidmCrossSection` | REAL | Cross-section per unit mass, **in code units** (see below). |
| `SidmDesNumNgb` | REAL | Target number of DM neighbours for the SIDM smoothing length. |
| `SidmDesNumNgbDev` | REAL | Allowed deviation from `SidmDesNumNgb`; the smoothing-length iteration stops when the count is within this tolerance. |

### Unit conversion for `SidmCrossSection`

```
[sigma/m]_code = (sigma/m)_cgs * UnitMass_in_g / UnitLength_in_cm^2
```

The factor depends on the unit system in your parameter file, so recompute it whenever `UnitLength_in_cm` or
`UnitMass_in_g` change. For the common conventions (1 cm^2/g in cgs):

- Mpc/h and 1e10 Msun/h (cosmological boxes): multiply by 2.089e-6.
- kpc and 1e10 Msun (isolated halo): multiply by 2.089.

Getting this wrong is easy and costly: a value of 0.1 in the cosmological convention is about 47,900 cm^2/g.

## Behaviour to be aware of

- The SIDM density, and the dedicated tree it uses, are rebuilt together with domain decomposition, not every
  timestep.
- Only the lower-ID particle of a pair evaluates a scattering, so the double-counting factor of 1/2 in the source
  formula is not applied.
- Cross-task scattering partners are supported. Between candidate gathering and kick delivery a remote particle
  could in principle have been kicked by another local scatter; this is an accepted approximation.
- The scattering time interval goes through cosmological time, not a raw timebin multiply.
- Dark-matter side data lives in a `DMSP[]` array that follows `P[]` through domain exchange, restarts and
  snapshots. Four fields are written to snapshots: `SidmDensity`, `SidmHsml`, `SidmNumNgb`, `SidmVelDisp`.
- When changing a file that also has `#ifdef STARS` blocks, build with `STARS` both on and off.

## Diagnostics in the log

The `SIDM_SCATTER:` log line reports scatter counts (including `cross_task=N`), the maximum momentum and
energy conservation error split by local and cross-task partners, and how often the SIDM term in the timestep
criterion was the binding one (`SidmTimestepChecks`, `SidmTimestepBinding`).
