# Examples

Examples live in `examples/<name>/`. Two kinds exist:

- **Test-suite examples** have a `create.py` (generates initial conditions) and `check.py` (compares the result
  with a reference). They run through `test.sh`; see [simtypes.md](./simtypes.md).
- **Run-by-hand examples** download external ICs and tables with `get_ics.sh` and `get_tables.sh` and have their
  own `README.md`. They are not part of `test.sh`. Run the scripts from inside the example directory, since
  downloads go to `./ICs/`, `./grackle_data/` and `./star_tables/` relative to the current directory.

## Hydrodynamics and MHD tests (test suite)

`shocktube_1d`, `interacting_blastwaves_1d`, `wave_1d`, `alfven_wave_1d`, `mhd_shocktube_1d`,
`polytrope_1d_spherical`, `gresho_2d`, `noh_2d`, `noh_3d`, `yee_2d`, `current_sheet_2d`.

## Gravity and galaxy tests (test suite)

| Example | What it runs |
|---|---|
| `cosmo_box_gravity_only_3d` | Cosmological box, gravity only. |
| `cosmo_zoom_gravity_only_3d` | Cosmological zoom, gravity only. |
| `isolated_galaxy_collisionless_3d` | Isolated collisionless galaxy (GalIC ICs). |
| `cosmo_box_star_formation_3d` | Small cosmological box with the default `EEOS_SF` star formation. |
| `galaxy_merger_star_formation_3d` | Galaxy merger with star formation. |

## Star formation and feedback (run by hand)

| Example | What it runs |
|---|---|
| `agora_disc_star_formation_3d` | AGORA isolated disc: `AGORA_SF`, Grackle, winds and supernovae. |
| `cosmo_L10N128_star_formation_3d` | Small cosmological box (10 Mpc/h, 128^3) with cooling, star formation and feedback. |
| `cosmo_zoom_L50N128_star_formation_3d` | Cosmological zoom with star formation (MUSIC2 configuration included). |

Each of these has a `README.md` with the download, build and run commands.

## Examples on other branches

Not yet on `main`: `agora_disc_bh_star_formation_3d` (AGORA disc with a central black hole, branch
`agora_disc_bh_example`) and `bh_jet_stratified_3d` (branch `jet_launching_condition`). Add them here when they
merge.

## Adding an example

Follow the guidelines in [development.md](./development.md). For a run-by-hand example, include a `README.md`
that says what is downloaded, which parameters differ from a related example, and what has and has not been
checked.
