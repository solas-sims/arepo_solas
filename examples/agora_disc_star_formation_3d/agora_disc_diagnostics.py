#!/usr/bin/env python3
"""
AGORA disc diagnostics command-line script.

Computes radial surface density, scale height, and vertical velocity
dispersion profiles for a chosen particle type, for one or more Arepo HDF5
snapshots. Profiles from multiple snapshots are overlaid on a single set of
axes; the smoothed-grid projection plot is always produced, one image per
snapshot.

Dependencies
------------
If the `analysistools` package is installed, it is used directly
(SnapshotTools, ProfileTools, GriddingTools). Otherwise this script falls
back to:
  - h5py, for reading the Arepo HDF5 snapshot directly
  - local ports of ProfileTools and GriddingTools (ported from
    profile_tools.py / gridding_tools.py), so the script has no dependency
    on analysistools itself and can be shared/run standalone.

The fallback GriddingTools drops the `numba` dependency (vectorized NumPy
is used instead of JIT'd loops) and omits the Voronoi/nearest-neighbour
slice helpers, which this script doesn't use.

Examples:
    # single snapshot
    python agora_disc_diagnostics.py --snapshot snap_200.hdf5

    # multiple snapshots, overlaid profiles + per-snapshot grid images
    python agora_disc_diagnostics.py \
        --snapshot snap_100.hdf5 snap_150.hdf5 snap_200.hdf5 \
        --output disc_profiles.png \
        --grid-output "disc_grid_{name}.png" \
        --mask-radius 15 \
        --save-data disc_profiles.npz
"""

import argparse
import logging
import os
import sys

import matplotlib
matplotlib.use("Agg")  # safe for headless / CLI use
import matplotlib.pyplot as plt
import numpy as np

log = logging.getLogger("agora_disc_diagnostics")

# ----------------------------------------------------------------------
# Prefer analysistools if it's installed; otherwise fall back to local,
# dependency-light ports of the same classes plus a raw h5py reader.
# ----------------------------------------------------------------------
try:
    from analysistools import SnapshotTools, ProfileTools, GriddingTools
    HAVE_ANALYSISTOOLS = True
except ImportError:
    HAVE_ANALYSISTOOLS = False
    SnapshotTools = None

    # ------------------------------------------------------------------
    # Fallback: ProfileTools (ported from profile_tools.py)
    # ------------------------------------------------------------------
    class ProfileTools:
        """Local port of analysistools.ProfileTools (subset used here, plus
        volume_density/velocity_dispersion/profile for API parity)."""

        def __init__(self, **kwargs):
            self.comoving_units = kwargs.get("comoving_units", True)
            self.fr_cut = kwargs.get("fr_cut", 5)
            self.numbins = kwargs.get("numbins", 25)
            self.halo_id = kwargs.get("halo_id", 0)

        @staticmethod
        def relative_position(pos, centre=None):
            if centre is None:
                centre = np.zeros(3)
            return pos - centre

        def spherical_radius(self, pos, centre=None):
            dpos = self.relative_position(pos, centre)
            return np.sqrt(np.sum(dpos ** 2, axis=1))

        def cylindrical_radius(self, pos, centre=None):
            dpos = self.relative_position(pos, centre)
            return np.sqrt(dpos[:, 0] ** 2 + dpos[:, 1] ** 2)

        def radial_bins(self, rmin, rmax, nbins=None, logarithmic=True):
            if nbins is None:
                nbins = self.numbins
            if logarithmic:
                return np.logspace(np.log10(rmin), np.log10(rmax), nbins + 1)
            return np.linspace(rmin, rmax, nbins + 1)

        def bin_indices(self, r, bins):
            return np.digitize(r, bins) - 1

        def profile(self, r, quantity, bins, weights=None):
            index = self.bin_indices(r, bins)
            nbins = len(bins) - 1
            result = {
                "r": np.zeros(nbins),
                "mean": np.full(nbins, np.nan),
                "median": np.full(nbins, np.nan),
                "count": np.zeros(nbins),
            }
            for i in range(nbins):
                mask = index == i
                result["r"][i] = 0.5 * (bins[i] + bins[i + 1])
                if np.any(mask):
                    q = quantity[mask]
                    if weights is None:
                        result["mean"][i] = np.mean(q)
                    else:
                        result["mean"][i] = np.average(q, weights=weights[mask])
                    result["median"][i] = np.median(q)
                    result["count"][i] = np.sum(mask)
            return result

        def volume_density(self, pos, mass, centre, rmin, rmax, nbins=None):
            r = self.spherical_radius(pos, centre)
            bins = self.radial_bins(rmin, rmax, nbins)
            index = self.bin_indices(r, bins)

            rho, radius = [], []
            for i in range(len(bins) - 1):
                mask = index == i
                volume = (4 * np.pi / 3.0) * (bins[i + 1] ** 3 - bins[i] ** 3)
                rho.append(np.sum(mass[mask]) / volume)
                radius.append(0.5 * (bins[i] + bins[i + 1]))
            return {"r": np.array(radius), "density": np.array(rho)}

        def surface_density(self, pos, mass, centre, rmin, rmax, nbins=None):
            R = self.cylindrical_radius(pos, centre)
            bins = self.radial_bins(rmin, rmax, nbins)
            index = self.bin_indices(R, bins)

            sigma, radius = [], []
            for i in range(len(bins) - 1):
                mask = index == i
                area = np.pi * (bins[i + 1] ** 2 - bins[i] ** 2)
                sigma.append(np.sum(mass[mask]) / area)
                radius.append(0.5 * (bins[i] + bins[i + 1]))
            return {"r": np.array(radius), "density": np.array(sigma)}

        def velocity_dispersion(self, pos, vel, centre, rmin, rmax, nbins=None):
            r = self.spherical_radius(pos, centre)
            bins = self.radial_bins(rmin, rmax, nbins)
            index = self.bin_indices(r, bins)

            sigma, radius = [], []
            for i in range(len(bins) - 1):
                mask = index == i
                radius.append(0.5 * (bins[i] + bins[i + 1]))
                if np.sum(mask) > 5:
                    v = vel[mask]
                    mean = np.mean(v, axis=0)
                    sigma.append(np.sqrt(np.mean(np.sum(v * v, axis=1)) - np.sum(mean * mean)))
                else:
                    sigma.append(np.nan)
            return {"r": np.array(radius), "sigma": np.array(sigma)}

        def scale_height(self, pos, mass, centre, rmin, rmax, nbins=None):
            R = self.cylindrical_radius(pos, centre)
            bins = self.radial_bins(rmin, rmax, nbins)
            index = self.bin_indices(R, bins)

            hz, radius = [], []
            z = pos[:, 2] - centre[2]
            for i in range(len(bins) - 1):
                mask = index == i
                radius.append(0.5 * (bins[i] + bins[i + 1]))
                if np.any(mask):
                    hz.append(np.sqrt(np.average(z[mask] ** 2, weights=mass[mask])))
                else:
                    hz.append(np.nan)
            return {"r": np.array(radius), "scale_height": np.array(hz)}

        def vertical_velocity_dispersion(self, pos, vel, centre, rmin, rmax, nbins=None):
            R = self.cylindrical_radius(pos, centre)
            bins = self.radial_bins(rmin, rmax, nbins)
            index = self.bin_indices(R, bins)

            sigma_z, radius = [], []
            vz = vel[:, 2]
            for i in range(len(bins) - 1):
                mask = index == i
                radius.append(0.5 * (bins[i] + bins[i + 1]))
                if np.sum(mask) > 5:
                    v = vz[mask]
                    mean = np.mean(v, axis=0)
                    sigma_z.append(np.sqrt(np.mean(v * v) - np.sum(mean * mean)))
                else:
                    sigma_z.append(np.nan)
            return {"r": np.array(radius), "sigma_z": np.array(sigma_z)}

        def plot(self, profile, ykey, xlabel="Radius", ylabel=None,
                 xlog=True, ylog=True, label=None):
            plt.figure()
            plt.plot(profile["r"], profile[ykey], label=label)
            if xlog:
                plt.xscale("log")
            if ylog:
                plt.yscale("log")
            plt.xlabel(xlabel)
            if ylabel:
                plt.ylabel(ylabel)
            if label:
                plt.legend()
            plt.show()

    # ------------------------------------------------------------------
    # Fallback: GriddingTools (ported from gridding_tools.py)
    # numba-JIT'd ngp_assign/cic_assign replaced with vectorized NumPy.
    # Voronoi/nearest-neighbour slice helpers omitted (unused here).
    # ------------------------------------------------------------------
    class GriddingTools:
        def __init__(self):
            pass

        @staticmethod
        def ngp_assign(grid, coords, values):
            """Nearest-Grid-Point assignment (vectorized). Particles that
            fall outside the grid are skipped rather than corrupting memory."""
            idx = np.floor(coords).astype(np.int64)
            nx, ny, nz = grid.shape
            valid = (
                (idx[:, 0] >= 0) & (idx[:, 0] < nx) &
                (idx[:, 1] >= 0) & (idx[:, 1] < ny) &
                (idx[:, 2] >= 0) & (idx[:, 2] < nz)
            )
            np.add.at(grid, (idx[valid, 0], idx[valid, 1], idx[valid, 2]), values[valid])
            return grid

        @staticmethod
        def cic_assign(grid, coords, values):
            """Cloud-In-Cell assignment (vectorized), periodic in all axes,
            matching the wrap-around behaviour of the original implementation."""
            nx, ny, nz = grid.shape
            i = np.floor(coords[:, 0]).astype(np.int64)
            j = np.floor(coords[:, 1]).astype(np.int64)
            k = np.floor(coords[:, 2]).astype(np.int64)
            fx = coords[:, 0] - i
            fy = coords[:, 1] - j
            fz = coords[:, 2] - k

            corners = [
                (i, j, k, (1 - fx) * (1 - fy) * (1 - fz)),
                (i + 1, j, k, fx * (1 - fy) * (1 - fz)),
                (i, j + 1, k, (1 - fx) * fy * (1 - fz)),
                (i, j, k + 1, (1 - fx) * (1 - fy) * fz),
                (i + 1, j, k + 1, fx * (1 - fy) * fz),
                (i, j + 1, k + 1, (1 - fx) * fy * fz),
                (i + 1, j + 1, k, fx * fy * (1 - fz)),
                (i + 1, j + 1, k + 1, fx * fy * fz),
            ]
            for ii, jj, kk, w in corners:
                np.add.at(grid, (ii % nx, jj % ny, kk % nz), values * w)
            return grid

        def smooth_to_grid(self, positions, values, grid_size, grid_limits,
                            method="NGP", sigma=1.0, filter_sigma=None):
            """Assign particle values to a 3D grid. method: 'NGP', 'CIC', or
            'Gaussian' (CIC + Gaussian smoothing)."""
            dim = len(grid_size)
            grid = np.zeros(grid_size, dtype=float)

            spacing = [
                (grid_limits[2 * i + 1] - grid_limits[2 * i]) / grid_size[i]
                for i in range(dim)
            ]
            coords = [(positions[:, i] - grid_limits[2 * i]) / spacing[i] for i in range(dim)]
            coords = np.stack(coords, axis=1)

            method = method.upper()
            if method == "NGP":
                self.ngp_assign(grid, coords, values)
            elif method == "CIC":
                self.cic_assign(grid, coords, values)
            elif method == "GAUSSIAN":
                self.cic_assign(grid, coords, values)
                from scipy.ndimage import gaussian_filter
                grid = gaussian_filter(grid, sigma=sigma)
            else:
                raise ValueError(f"Unknown method '{method}'")

            if filter_sigma is not None:
                from scipy.ndimage import gaussian_filter
                grid = gaussian_filter(grid, sigma=filter_sigma)

            return grid

        def plot_3d_slice(self, grid_3d, grid_limits,
                           slice_axis='z', slice_index=None,
                           slice_width=None, slice_average=True,
                           mode='slice', projection='mean',
                           title="3D Grid Slice", cmap='viridis', figsize=(12, 4)):
            nx, ny, nz = grid_3d.shape

            if mode == 'slice':
                if slice_axis == 'z':
                    if slice_index is None:
                        slice_index = nz // 2
                    if slice_width is None:
                        slice_width = 1
                    start_idx = max(slice_index - slice_width // 2, 0)
                    end_idx = min(slice_index + slice_width // 2 + 1, nz)
                    data = (grid_3d[:, :, start_idx:end_idx].mean(axis=2) if slice_average
                            else grid_3d[:, :, start_idx:end_idx].sum(axis=2))
                    extent = [grid_limits[0], grid_limits[1], grid_limits[2], grid_limits[3]]
                    xlabel, ylabel, title_str = 'X', 'Y', f'XY slice (Z={slice_index})'

                elif slice_axis == 'y':
                    if slice_index is None:
                        slice_index = ny // 2
                    if slice_width is None:
                        slice_width = 1
                    start_idx = max(slice_index - slice_width // 2, 0)
                    end_idx = min(slice_index + slice_width // 2 + 1, ny)
                    data = (grid_3d[:, start_idx:end_idx, :].mean(axis=1) if slice_average
                            else grid_3d[:, start_idx:end_idx, :].sum(axis=1))
                    extent = [grid_limits[0], grid_limits[1], grid_limits[4], grid_limits[5]]
                    xlabel, ylabel, title_str = 'X', 'Z', f'XZ slice (Y={slice_index})'

                elif slice_axis == 'x':
                    if slice_index is None:
                        slice_index = nx // 2
                    if slice_width is None:
                        slice_width = 1
                    start_idx = max(slice_index - slice_width // 2, 0)
                    end_idx = min(slice_index + slice_width // 2 + 1, nx)
                    data = (grid_3d[start_idx:end_idx, :, :].mean(axis=0) if slice_average
                            else grid_3d[start_idx:end_idx, :, :].sum(axis=0))
                    extent = [grid_limits[2], grid_limits[3], grid_limits[4], grid_limits[5]]
                    xlabel, ylabel, title_str = 'Y', 'Z', f'YZ slice (X={slice_index})'
                else:
                    raise ValueError("slice_axis must be 'x', 'y', or 'z'")

            elif mode == 'projection':
                if slice_axis == 'z':
                    data = {'mean': grid_3d.mean, 'sum': grid_3d.sum, 'max': grid_3d.max}[projection](axis=2)
                    extent = [grid_limits[0], grid_limits[1], grid_limits[2], grid_limits[3]]
                    xlabel, ylabel, title_str = 'X', 'Y', f'XY projection ({projection} along Z)'
                elif slice_axis == 'y':
                    data = {'mean': grid_3d.mean, 'sum': grid_3d.sum, 'max': grid_3d.max}[projection](axis=1)
                    extent = [grid_limits[0], grid_limits[1], grid_limits[4], grid_limits[5]]
                    xlabel, ylabel, title_str = 'X', 'Z', f'XZ projection ({projection} along Y)'
                elif slice_axis == 'x':
                    data = {'mean': grid_3d.mean, 'sum': grid_3d.sum, 'max': grid_3d.max}[projection](axis=0)
                    extent = [grid_limits[2], grid_limits[3], grid_limits[4], grid_limits[5]]
                    xlabel, ylabel, title_str = 'Y', 'Z', f'YZ projection ({projection} along X)'
                else:
                    raise ValueError("slice_axis must be 'x', 'y', or 'z'")
            else:
                raise ValueError("mode must be 'slice' or 'projection'")

            fig, ax = plt.subplots(figsize=figsize)
            with np.errstate(divide="ignore"):
                im = ax.imshow(np.log10(data.T), origin='lower', extent=extent, cmap=cmap, aspect='auto')
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(title_str)
            fig.suptitle(title, fontsize=14)
            plt.colorbar(im, ax=ax)
            plt.tight_layout()
            return fig, ax

        def plot_3d_projections(self, grid_3d, grid_limits,
                                 mode='projection', projection='sum',
                                 cmap='viridis', figsize=(12, 4), title=None,
                                 slice_axis='z', slice_index=None,
                                 slice_width=None, slice_average=True):
            fig, axes = plt.subplots(1, 3, figsize=figsize)

            for ax, axis, lbl in zip(axes, ['z', 'y', 'x'], ['XY', 'XZ', 'YZ']):
                idx = None
                if mode == 'slice' and slice_index is not None:
                    idx = slice_index.get(axis, None) if isinstance(slice_index, dict) else slice_index

                _, single_ax = self.plot_3d_slice(
                    grid_3d, grid_limits,
                    slice_axis=axis,
                    mode=mode,
                    projection=projection,
                    cmap=cmap,
                    figsize=(5, 5),
                )

                im = single_ax.images[0]
                ax.imshow(im.get_array(), origin='lower', extent=im.get_extent(),
                          cmap=im.get_cmap() if hasattr(im, 'get_cmap') else cmap, aspect='auto')
                ax.set_title(lbl)
                ax.set_xlabel(single_ax.get_xlabel())
                ax.set_ylabel(single_ax.get_ylabel())
                plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                plt.close(single_ax.figure)

            if title:
                fig.suptitle(title, fontsize=14)

            plt.tight_layout()
            return fig, axes


def read_snapshot_h5py(path, ptype):
    """Read positions, velocities, mass, and box size for one particle type
    directly from an Arepo HDF5 snapshot, without analysistools."""
    import h5py

    with h5py.File(path, "r") as f:
        header = f["Header"].attrs

        boxsize = np.asarray(header["BoxSize"], dtype=float)
        if boxsize.ndim == 0:
            boxsize = np.full(3, float(boxsize))
        elif boxsize.size == 1:
            boxsize = np.full(3, float(boxsize.reshape(-1)[0]))

        group_key = f"PartType{ptype}"
        if group_key not in f:
            raise KeyError(f"'{group_key}' not present in {path}")
        grp = f[group_key]

        if "Coordinates" not in grp:
            raise KeyError(f"'{group_key}/Coordinates' not present in {path}")
        pos = grp["Coordinates"][()]

        vel = grp["Velocities"][()] if "Velocities" in grp else None

        if "Masses" in grp:
            mass = grp["Masses"][()]
        else:
            mass_table = np.asarray(header.get("MassTable", None), dtype=float) \
                if header.get("MassTable", None) is not None else None
            if mass_table is not None and ptype < len(mass_table) and mass_table[ptype] > 0:
                mass = np.full(pos.shape[0], mass_table[ptype])
            else:
                raise ValueError(
                    f"No per-particle masses and no usable MassTable entry "
                    f"for ptype {ptype} in {path}"
                )

    return pos, vel, mass, boxsize


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute AGORA disc diagnostics from one or more Arepo snapshots"
    )
    parser.add_argument(
        "--snapshot", required=True, nargs="+",
        help="One or more Arepo HDF5 snapshots. Radial profiles are overlaid on "
             "one figure; gridded images are plotted independently per snapshot.",
    )
    parser.add_argument("--output", default="disc_diagnostics.png", help="Output figure for overlaid radial profiles")
    parser.add_argument("--ptype", type=int, default=0, help="Particle type to analyse (default: gas)")
    parser.add_argument("--nbins", type=int, default=100, help="Number of radial bins")
    parser.add_argument("--rmin", type=float, default=0.1, help="Minimum radius [kpc]")
    parser.add_argument("--rmax", type=float, default=50.0, help="Maximum radius [kpc]")
    parser.add_argument(
        "--centre", nargs=3, type=float, default=None, metavar=("X", "Y", "Z"),
        help="Disc centre, applied to every snapshot. Default: each snapshot's own box centre",
    )
    parser.add_argument(
        "--mask-radius", type=float, default=None,
        help="If set, restrict particles to within this cylindrical radius (and "
             "--mask-zheight above/below the midplane) of the centre before "
             "profiling. Speeds things up and excludes off-disc contamination. "
             "Default: off (use the full particle set).",
    )
    parser.add_argument(
        "--mask-zheight", type=float, default=2.5,
        help="Half-height [kpc] of the mask box above/below the midplane "
             "(only used with --mask-radius). Default: 2.5",
    )
    parser.add_argument("--ngrid", type=int, default=512, help="Grid resolution per axis")
    parser.add_argument(
        "--smooth-sigma", type=float, default=2.0,
        help="Gaussian smoothing width [grid cells] applied to the CIC-deposited "
             "grid before plotting. Default: 2.0",
    )
    parser.add_argument("--lgrid", type=float, default=10.0, help="Grid half-extent in x/y [kpc]")
    parser.add_argument(
        "--lzgrid", type=float, default=None,
        help="Grid half-extent in z [kpc]. Default: lgrid/4",
    )
    parser.add_argument(
        "--grid-output", default="disc_grid_{name}.png",
        help="Output figure template for the per-snapshot grid projection. "
             "May contain {name} (snapshot filename stem) and/or {index}. "
             "If neither placeholder is present, the snapshot name is appended "
             "automatically so multiple snapshots don't overwrite each other.",
    )
    parser.add_argument("--save-data", default=None, help="Optional .npz path to save the computed profiles (all snapshots, one file)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    return parser.parse_args()


def snapshot_name(path):
    return os.path.splitext(os.path.basename(path))[0]


def format_output_path(template, name, index):
    if "{name}" in template or "{index}" in template:
        return template.format(name=name, index=index)
    base, ext = os.path.splitext(template)
    return f"{base}_{name}{ext}"


def load_particles(snapshot_path, args):
    if HAVE_ANALYSISTOOLS:
        snap = SnapshotTools(snapfileformat="HDF5", convention="Arepo")
        try:
            data = snap.read_snapshot(snapshot_path)
        except (OSError, IOError) as exc:
            log.error("Could not read snapshot '%s': %s", snapshot_path, exc)
            sys.exit(1)

        boxsize = data.box_size
        num_part = data.num_part_total
        if not (0 <= args.ptype < len(num_part)):
            log.error("ptype %d out of range (snapshot has %d particle types)", args.ptype, len(num_part))
            sys.exit(1)

        ij = np.append(0, np.cumsum(num_part, dtype=np.int32))
        pos = data.pos[ij[args.ptype]:ij[args.ptype + 1]]
        vel = data.vel[ij[args.ptype]:ij[args.ptype + 1]]
        mass = data.mass[ij[args.ptype]:ij[args.ptype + 1]]
    else:
        try:
            pos, vel, mass, boxsize = read_snapshot_h5py(snapshot_path, args.ptype)
        except (OSError, IOError, KeyError, ValueError) as exc:
            log.error("Could not read snapshot '%s': %s", snapshot_path, exc)
            sys.exit(1)

        if vel is None:
            log.error("'PartType%d/Velocities' not present in %s (needed for kinematic profiles)",
                       args.ptype, snapshot_path)
            sys.exit(1)

    centre = np.ones(3) * boxsize / 2 if args.centre is None else np.array(args.centre)

    log.info("Loaded %d particles of type %d from %s", len(pos), args.ptype, snapshot_path)
    log.info("Centre = %s", centre)

    return pos, vel, mass, centre


def apply_mask(pos, vel, mass, centre, args):
    """Restrict to a cylindrical region around the centre. Optional: improves
    speed and avoids off-disc contamination when enabled via --mask-radius."""
    if args.mask_radius is None:
        return pos, vel, mass

    dxy = pos[:, :2] - centre[:2]
    r = np.hypot(dxy[:, 0], dxy[:, 1])
    dz = np.abs(pos[:, 2] - centre[2])

    mask = (r <= args.mask_radius) & (dz <= args.mask_zheight)
    log.info("Mask kept %d / %d particles", mask.sum(), len(pos))

    return pos[mask], vel[mask], mass[mask]


def compute_profiles(pos, vel, mass, centre, args):
    if args.rmin >= args.rmax:
        log.error("--rmin (%g) must be less than --rmax (%g)", args.rmin, args.rmax)
        sys.exit(1)
    if args.nbins <= 0:
        log.error("--nbins must be positive")
        sys.exit(1)

    profiler = ProfileTools(numbins=args.nbins)

    sigma = profiler.surface_density(pos, mass, centre, args.rmin, args.rmax)

    inner_rmax = min(args.rmax, 10)
    inner_nbins = min(args.nbins, 30)

    hz = profiler.scale_height(pos, mass, centre, args.rmin, inner_rmax, nbins=inner_nbins)
    sigmaz = profiler.vertical_velocity_dispersion(pos, vel, centre, args.rmin, inner_rmax, nbins=inner_nbins)

    return sigma, hz, sigmaz


def plot_overlaid_profiles(results, output):
    """results: list of dicts with keys name, sigma, hz, sigmaz."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    for res in results:
        axes[0].plot(res["sigma"]["r"], res["sigma"]["density"], label=res["name"])
        axes[1].plot(res["hz"]["r"], res["hz"]["scale_height"], label=res["name"])
        axes[2].plot(res["sigmaz"]["r"], res["sigmaz"]["sigma_z"], label=res["name"])

    axes[0].set_xlabel("R [kpc]")
    axes[0].set_ylabel(r"$\Sigma(R)$")
    axes[0].set_yscale("log")

    axes[1].set_xlabel("R [kpc]")
    axes[1].set_ylabel("h(z) [kpc]")

    axes[2].set_xlabel("R [kpc]")
    axes[2].set_ylabel(r"$\sigma_z$ [km/s]")

    if len(results) > 1:
        axes[0].legend(fontsize="small")

    fig.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)
    log.info("Saved %s", output)

def plot_grid(pos, mass, centre, args, output_path):
    """Smoothed-grid projection plot for a single snapshot. Always run,
    regardless of --mask-radius (uses the full particle set, like the
    source notebook)."""
    gridding = GriddingTools()

    lzgrid = args.lzgrid if args.lzgrid is not None else args.lgrid / 4
    grid_size = args.ngrid * np.ones(3, dtype=np.int64)
    grid_limits = np.array([
        centre[0] - args.lgrid, centre[0] + args.lgrid,
        centre[1] - args.lgrid, centre[1] + args.lgrid,
        centre[2] - lzgrid, centre[2] + lzgrid,
    ])

    box_mask = (
        (pos[:, 0] >= grid_limits[0]) & (pos[:, 0] < grid_limits[1]) &
        (pos[:, 1] >= grid_limits[2]) & (pos[:, 1] < grid_limits[3]) &
        (pos[:, 2] >= grid_limits[4]) & (pos[:, 2] < grid_limits[5])
    )

    smoothed_grid = gridding.smooth_to_grid(
        positions=pos[box_mask],
        values=mass[box_mask],
        grid_size=grid_size,
        grid_limits=grid_limits,
        method="Gaussian",
        sigma=args.smooth_sigma,
    )

    fig, _ = gridding.plot_3d_projections(smoothed_grid, grid_limits, projection="mean")
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    log.info("Saved %s", output_path)


def save_data(results, path):
    arrays = {}
    for res in results:
        name = res["name"]
        arrays[f"{name}__sigma_r"] = res["sigma"]["r"]
        arrays[f"{name}__sigma_density"] = res["sigma"]["density"]
        arrays[f"{name}__hz_r"] = res["hz"]["r"]
        arrays[f"{name}__hz_scale_height"] = res["hz"]["scale_height"]
        arrays[f"{name}__sigmaz_r"] = res["sigmaz"]["r"]
        arrays[f"{name}__sigmaz_sigma_z"] = res["sigmaz"]["sigma_z"]
    np.savez(path, **arrays)
    log.info("Saved profile data to %s", path)

def main():
    args = parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    log.info("Using %s for snapshot I/O and tools",
              "analysistools" if HAVE_ANALYSISTOOLS else "h5py + local ProfileTools/GriddingTools fallback")

    results = []

    for index, snapshot_path in enumerate(args.snapshot):
        name = snapshot_name(snapshot_path)

        pos, vel, mass, centre = load_particles(snapshot_path, args)

        grid_out = format_output_path(args.grid_output, name, index)
        plot_grid(pos, mass, centre, args, grid_out)

        pos_m, vel_m, mass_m = apply_mask(pos, vel, mass, centre, args)
        sigma, hz, sigmaz = compute_profiles(pos_m, vel_m, mass_m, centre, args)

        results.append({"name": name, "sigma": sigma, "hz": hz, "sigmaz": sigmaz})

    plot_overlaid_profiles(results, args.output)

    if args.save_data:
        save_data(results, args.save_data)


if __name__ == "__main__":
    main()
