#!/usr/bin/env python3
"""
AGORA disc diagnostics command-line script.

Computes radial surface density, scale height, and vertical velocity
dispersion profiles for a chosen particle type, for one or more Arepo HDF5
snapshots. Profiles from multiple snapshots are overlaid on a single set of
axes; the smoothed-grid projection plot is always produced, one image per
snapshot.

Also reads the gas cells' metal mass (`PassiveScalars`, METALS_INDEX==0)
and, if present, dust mass (`DustMass`, only written when the run was
built with `DUST`) directly from PartType0 -- these are gas-only Arepo
scalar fields, independent of whatever `--ptype` is requested for the
main mass/kinematic profiles. For each snapshot this produces the same
smoothed-grid XY/XZ/YZ projection plots as gas mass gets, plus a second
overlaid-profile figure of metal and dust surface density vs. radius,
mirroring the gas surface-density panel. The dust grid/profile is skipped
(with a log message, not an error) if the snapshot has no DustMass field.

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
    # (gas, metals, and dust -- dust grid/profile auto-skipped if the run
    # has no DustMass field)
    python agora_disc_diagnostics.py \
        --snapshot snap_100.hdf5 snap_150.hdf5 snap_200.hdf5 \
        --output disc_profiles.png \
        --grid-output "disc_grid_{name}.png" \
        --metals-grid-output "disc_grid_metals_{name}.png" \
        --dust-grid-output "disc_grid_dust_{name}.png" \
        --metals-dust-output "disc_metals_dust_profiles.png" \
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

        @staticmethod
        def _cubic_spline_kernel(q):
            """Normalized 3D cubic spline shape function (Monaghan &
            Lattanzio 1985), compact support q<2. `q` may be an array."""
            q = np.asarray(q)
            w = np.where(
                q < 1.0, 1.0 - 1.5 * q ** 2 + 0.75 * q ** 3,
                np.where(q < 2.0, 0.25 * (2.0 - q) ** 3, 0.0),
            )
            return w

        def sph_assign(self, grid, positions, values, smoothing_lengths, grid_limits, grid_size):
            """Per-particle SPH-kernel (cubic spline) deposition -- a
            dependency-light (non-numba) port of
            analysistools.gridding_tools.GriddingTools.sph_assign. Uses each
            particle's own smoothing length rather than a single global
            smoothing scale; mass-conserving per particle (weights
            renormalized to sum to `values[p]` over that particle's own
            affected cells); non-periodic (kernel support is clipped to the
            grid, not wrapped); correct for anisotropic grid spacing (works
            in physical coordinates throughout). Each particle's local
            bounding box (~(4h/spacing)**3 cells) is handled with a small
            vectorized numpy computation -- there's no global grid-cells
            array (infeasible at typical --ngrid resolutions), so this is
            still a Python-level loop over particles, just a cheap one per
            particle rather than a full triple-nested loop.
            """
            nx, ny, nz = grid_size
            xmin, xmax, ymin, ymax, zmin, zmax = grid_limits
            dx, dy, dz = (xmax - xmin) / nx, (ymax - ymin) / ny, (zmax - zmin) / nz

            for p in range(len(values)):
                h = smoothing_lengths[p]
                x, y, z = positions[p]

                if h <= 0:
                    i, j, k = (int(np.floor((x - xmin) / dx)),
                               int(np.floor((y - ymin) / dy)),
                               int(np.floor((z - zmin) / dz)))
                    if 0 <= i < nx and 0 <= j < ny and 0 <= k < nz:
                        grid[i, j, k] += values[p]
                    continue

                rmax = 2.0 * h
                imin = max(int(np.floor((x - xmin - rmax) / dx)), 0)
                imax = min(int(np.ceil((x - xmin + rmax) / dx)), nx - 1)
                jmin = max(int(np.floor((y - ymin - rmax) / dy)), 0)
                jmax = min(int(np.ceil((y - ymin + rmax) / dy)), ny - 1)
                kmin = max(int(np.floor((z - zmin - rmax) / dz)), 0)
                kmax = min(int(np.ceil((z - zmin + rmax) / dz)), nz - 1)
                if imin > imax or jmin > jmax or kmin > kmax:
                    continue

                cx = xmin + (np.arange(imin, imax + 1) + 0.5) * dx
                cy = ymin + (np.arange(jmin, jmax + 1) + 0.5) * dy
                cz = zmin + (np.arange(kmin, kmax + 1) + 0.5) * dz
                ddx, ddy, ddz = np.meshgrid(cx - x, cy - y, cz - z, indexing="ij")
                r = np.sqrt(ddx ** 2 + ddy ** 2 + ddz ** 2)
                w = self._cubic_spline_kernel(r / h)
                wsum = w.sum()

                if wsum <= 0:
                    i, j, k = (int(np.floor((x - xmin) / dx)),
                               int(np.floor((y - ymin) / dy)),
                               int(np.floor((z - zmin) / dz)))
                    if 0 <= i < nx and 0 <= j < ny and 0 <= k < nz:
                        grid[i, j, k] += values[p]
                    continue

                grid[imin:imax + 1, jmin:jmax + 1, kmin:kmax + 1] += values[p] * (w / wsum)

            return grid

        def smooth_to_grid(self, positions, values, grid_size, grid_limits,
                            method="NGP", sigma=1.0, filter_sigma=None,
                            smoothing_lengths=None):
            """Assign particle values to a 3D grid. method: 'NGP', 'CIC',
            'Gaussian' (CIC + Gaussian smoothing), or 'SPH' (per-particle
            kernel deposition, requires smoothing_lengths)."""
            dim = len(grid_size)
            grid = np.zeros(grid_size, dtype=float)

            method = method.upper()

            if method == "SPH":
                if smoothing_lengths is None:
                    raise ValueError("method='SPH' requires smoothing_lengths")
                if dim != 3:
                    raise ValueError("method='SPH' only supports 3D grids")
                self.sph_assign(grid, positions, values, smoothing_lengths, grid_limits, grid_size)
                if filter_sigma is not None:
                    from scipy.ndimage import gaussian_filter
                    grid = gaussian_filter(grid, sigma=filter_sigma)
                return grid

            spacing = [
                (grid_limits[2 * i + 1] - grid_limits[2 * i]) / grid_size[i]
                for i in range(dim)
            ]
            coords = [(positions[:, i] - grid_limits[2 * i]) / spacing[i] for i in range(dim)]
            coords = np.stack(coords, axis=1)

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

        @staticmethod
        def _collapse(grid_3d, grid_limits, slice_axis='z', slice_index=None,
                      slice_width=None, slice_average=True, mode='projection',
                      projection='mean', normalize_by_area=False):
            """Local port of analysistools.GriddingTools._collapse -- see
            its docstring for normalize_by_area's units caveats."""
            nx, ny, nz = grid_3d.shape

            if mode == 'slice':
                if slice_axis == 'z':
                    n = nz
                    if slice_index is None:
                        slice_index = n // 2
                    if slice_width is None:
                        slice_width = 1
                    start_idx = max(slice_index - slice_width // 2, 0)
                    end_idx = min(slice_index + slice_width // 2 + 1, n)
                    sub = grid_3d[:, :, start_idx:end_idx]
                    data = sub.mean(axis=2) if slice_average else sub.sum(axis=2)
                    extent = [grid_limits[0], grid_limits[1], grid_limits[2], grid_limits[3]]
                    xlabel, ylabel, title_str = 'X', 'Y', f'XY slice (Z={slice_index})'
                elif slice_axis == 'y':
                    n = ny
                    if slice_index is None:
                        slice_index = n // 2
                    if slice_width is None:
                        slice_width = 1
                    start_idx = max(slice_index - slice_width // 2, 0)
                    end_idx = min(slice_index + slice_width // 2 + 1, n)
                    sub = grid_3d[:, start_idx:end_idx, :]
                    data = sub.mean(axis=1) if slice_average else sub.sum(axis=1)
                    extent = [grid_limits[0], grid_limits[1], grid_limits[4], grid_limits[5]]
                    xlabel, ylabel, title_str = 'X', 'Z', f'XZ slice (Y={slice_index})'
                elif slice_axis == 'x':
                    n = nx
                    if slice_index is None:
                        slice_index = n // 2
                    if slice_width is None:
                        slice_width = 1
                    start_idx = max(slice_index - slice_width // 2, 0)
                    end_idx = min(slice_index + slice_width // 2 + 1, n)
                    sub = grid_3d[start_idx:end_idx, :, :]
                    data = sub.mean(axis=0) if slice_average else sub.sum(axis=0)
                    extent = [grid_limits[2], grid_limits[3], grid_limits[4], grid_limits[5]]
                    xlabel, ylabel, title_str = 'Y', 'Z', f'YZ slice (X={slice_index})'
                else:
                    raise ValueError("slice_axis must be 'x', 'y', or 'z'")

            elif mode == 'projection':
                collapse_axis = {'z': 2, 'y': 1, 'x': 0}.get(slice_axis)
                if collapse_axis is None:
                    raise ValueError("slice_axis must be 'x', 'y', or 'z'")
                data = {'mean': grid_3d.mean, 'sum': grid_3d.sum, 'max': grid_3d.max}[projection](axis=collapse_axis)

                if slice_axis == 'z':
                    extent = [grid_limits[0], grid_limits[1], grid_limits[2], grid_limits[3]]
                    xlabel, ylabel, title_str = 'X', 'Y', f'XY projection ({projection} along Z)'
                elif slice_axis == 'y':
                    extent = [grid_limits[0], grid_limits[1], grid_limits[4], grid_limits[5]]
                    xlabel, ylabel, title_str = 'X', 'Z', f'XZ projection ({projection} along Y)'
                else:
                    extent = [grid_limits[2], grid_limits[3], grid_limits[4], grid_limits[5]]
                    xlabel, ylabel, title_str = 'Y', 'Z', f'YZ projection ({projection} along X)'

                if normalize_by_area:
                    pixel_area = ((extent[1] - extent[0]) / data.shape[0]) * \
                                 ((extent[3] - extent[2]) / data.shape[1])
                    data = data / pixel_area
            else:
                raise ValueError("mode must be 'slice' or 'projection'")

            return data, extent, xlabel, ylabel, title_str

        @staticmethod
        def _plot_panel(ax, data, extent, xlabel, ylabel, title_str, cmap='viridis',
                         vmin=None, vmax=None):
            with np.errstate(divide='ignore', invalid='ignore'):
                plotted = np.log10(data.T)
            im = ax.imshow(plotted, origin='lower', extent=extent, cmap=cmap,
                            aspect='auto', vmin=vmin, vmax=vmax)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(title_str)
            return im

        @staticmethod
        def _shared_clim(data_list):
            with np.errstate(divide='ignore', invalid='ignore'):
                all_vals = np.concatenate([np.log10(d).ravel() for d in data_list])
            finite = all_vals[np.isfinite(all_vals)]
            if finite.size == 0:
                return None, None
            return float(finite.min()), float(finite.max())

        def plot_3d_slice(self, grid_3d, grid_limits,
                           slice_axis='z', slice_index=None,
                           slice_width=None, slice_average=True,
                           mode='slice', projection='mean',
                           title="3D Grid Slice", cmap='viridis', figsize=(12, 4),
                           vmin=None, vmax=None, normalize_by_area=False, cbar_label=None):
            data, extent, xlabel, ylabel, title_str = self._collapse(
                grid_3d, grid_limits, slice_axis=slice_axis, slice_index=slice_index,
                slice_width=slice_width, slice_average=slice_average, mode=mode,
                projection=projection, normalize_by_area=normalize_by_area,
            )
            fig, ax = plt.subplots(figsize=figsize)
            im = self._plot_panel(ax, data, extent, xlabel, ylabel, title_str,
                                   cmap=cmap, vmin=vmin, vmax=vmax)
            fig.suptitle(title, fontsize=14)
            cbar = plt.colorbar(im, ax=ax)
            if cbar_label is not None:
                cbar.set_label(cbar_label)
            plt.tight_layout()
            return fig, ax

        def plot_3d_projections(self, grid_3d, grid_limits,
                                 mode='projection', projection='sum',
                                 cmap='viridis', figsize=(12, 4), title=None,
                                 slice_axis='z', slice_index=None,
                                 slice_width=None, slice_average=True,
                                 vmin=None, vmax=None, normalize_by_area=False,
                                 cbar_label=None):
            fig, axes = plt.subplots(1, 3, figsize=figsize)

            panels = []
            for axis in ['z', 'y', 'x']:
                idx = slice_index
                if mode == 'slice' and isinstance(slice_index, dict):
                    idx = slice_index.get(axis, None)
                panels.append(self._collapse(
                    grid_3d, grid_limits, slice_axis=axis, slice_index=idx,
                    slice_width=slice_width, slice_average=slice_average, mode=mode,
                    projection=projection, normalize_by_area=normalize_by_area,
                ))

            if vmin is None or vmax is None:
                auto_vmin, auto_vmax = self._shared_clim([p[0] for p in panels])
                vmin = auto_vmin if vmin is None else vmin
                vmax = auto_vmax if vmax is None else vmax

            for ax, lbl, (data, extent, xlabel, ylabel, _) in zip(axes, ['XY', 'XZ', 'YZ'], panels):
                im = self._plot_panel(ax, data, extent, xlabel, ylabel, lbl,
                                       cmap=cmap, vmin=vmin, vmax=vmax)
                cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                if cbar_label is not None:
                    cbar.set_label(cbar_label)

            if title:
                fig.suptitle(title, fontsize=14)

            plt.tight_layout()
            return fig, axes

        def plot_ratio_projections(self, numerator_3d, denominator_3d, grid_limits,
                                   projection='sum', cmap='viridis', figsize=(12, 4),
                                   title=None, vmin=None, vmax=None,
                                   min_denominator=0.0, cbar_label=None):
            if numerator_3d.shape != denominator_3d.shape:
                raise ValueError("numerator_3d and denominator_3d must have the same shape")

            fig, axes = plt.subplots(1, 3, figsize=figsize)

            panels = []
            for axis in ['z', 'y', 'x']:
                num, extent, xlabel, ylabel, _ = self._collapse(
                    numerator_3d, grid_limits, slice_axis=axis, mode='projection',
                    projection=projection,
                )
                den, _, _, _, _ = self._collapse(
                    denominator_3d, grid_limits, slice_axis=axis, mode='projection',
                    projection=projection,
                )
                with np.errstate(divide='ignore', invalid='ignore'):
                    ratio = np.where(den > min_denominator, num / den, np.nan)
                panels.append((ratio, extent, xlabel, ylabel))

            if vmin is None or vmax is None:
                auto_vmin, auto_vmax = self._shared_clim([p[0] for p in panels])
                vmin = auto_vmin if vmin is None else vmin
                vmax = auto_vmax if vmax is None else vmax

            for ax, lbl, (data, extent, xlabel, ylabel) in zip(axes, ['XY', 'XZ', 'YZ'], panels):
                im = self._plot_panel(ax, data, extent, xlabel, ylabel, lbl,
                                       cmap=cmap, vmin=vmin, vmax=vmax)
                cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                if cbar_label is not None:
                    cbar.set_label(cbar_label)

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


def read_gas_scalars_h5py(path):
    """Read gas (PartType0) positions, mass, metal mass, dust mass (if
    present), and each cell's smoothing length directly from an Arepo HDF5
    snapshot. Always reads via h5py, even when analysistools is available,
    since PassiveScalars/DustMass are Arepo-specific fields not exposed by
    a generic snapshot reader.

    Smoothing length is each Voronoi cell's effective (equal-volume-sphere)
    radius, (3*Volume/(4*pi))**(1/3), computed from Density (not a native
    SPH kernel length -- Arepo is a moving-mesh code -- but the natural
    per-cell analogue used for SPH-kernel-style grid deposition, see
    GriddingTools.sph_assign in analysistools). None if the snapshot has no
    Density field (e.g. an IC file with ADDBACKGROUNDGRID but no run yet).

    Returns (pos, mass, metal_mass, dust_mass, smoothing_length); dust_mass
    is None if the run wasn't built with DUST.
    """
    import h5py

    with h5py.File(path, "r") as f:
        if "PartType0" not in f:
            raise KeyError(f"'PartType0' not present in {path}")
        grp = f["PartType0"]

        pos = grp["Coordinates"][()]
        mass = grp["Masses"][()]

        metal_mass = None
        if "PassiveScalars" in grp:
            passive_scalars = grp["PassiveScalars"][()]
            metallicity = passive_scalars[:, 0] if passive_scalars.ndim == 2 else passive_scalars
            metal_mass = metallicity * mass

        dust_mass = grp["DustMass"][()] if "DustMass" in grp else None

        smoothing_length = None
        if "Density" in grp:
            volume = mass / grp["Density"][()]
            smoothing_length = (3.0 * volume / (4.0 * np.pi)) ** (1.0 / 3.0)

    return pos, mass, metal_mass, dust_mass, smoothing_length


SOLAR_MASS_G = 1.989e33
PC_PER_KPC = 1.0e3
PC2_PER_KPC2 = PC_PER_KPC ** 2


def surface_density_value_scale(path):
    """Factor to pre-multiply a mass-like `values` array (gas/metal/dust
    mass, in code mass units) by before SPH-depositing it and area-
    normalizing (see plot_grid), so the resulting column density comes out
    directly in Msun/pc^2 -- rather than in [code mass unit] / [grid_limits
    unit]^2, which GriddingTools' normalize_by_area has no way to know is
    physically meaningful. Assumes grid_limits/positions are in kpc (true
    for this example's snapshots, UnitLength_in_cm = 1 kpc).

    Reads UnitMass_in_g from the snapshot's own /Header (written by every
    Arepo run, not just this fork) rather than assuming a fixed value, so
    this keeps working if the example's units ever change.
    """
    import h5py

    with h5py.File(path, "r") as f:
        unit_mass_in_g = f["/Header"].attrs["UnitMass_in_g"]

    return (unit_mass_in_g / SOLAR_MASS_G) / PC2_PER_KPC2


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
        help="Gaussian smoothing width [grid cells] for the CIC+Gaussian "
             "grid-projection fallback -- only used when a snapshot has no "
             "Density field (e.g. an unrun IC) or --ptype isn't gas, so "
             "per-cell smoothing lengths aren't available. Otherwise each "
             "gas cell's own effective radius is used (SPH-kernel "
             "deposition), and this option has no effect. Default: 2.0",
    )
    parser.add_argument("--lgrid", type=float, default=10.0, help="Grid half-extent in x/y [kpc]")
    parser.add_argument(
        "--lzgrid", type=float, default=None,
        help="Grid half-extent in z [kpc]. Default: lgrid/4",
    )
    parser.add_argument(
        "--grid-output", default="disc_grid_{name}.png",
        help="Output figure template for the per-snapshot gas mass grid projection. "
             "May contain {name} (snapshot filename stem) and/or {index}. "
             "If neither placeholder is present, the snapshot name is appended "
             "automatically so multiple snapshots don't overwrite each other.",
    )
    parser.add_argument(
        "--metals-grid-output", default="disc_grid_metals_{name}.png",
        help="Output figure template for the per-snapshot gas metal-mass grid "
             "projection. Same {name}/{index} substitution as --grid-output.",
    )
    parser.add_argument(
        "--dust-grid-output", default="disc_grid_dust_{name}.png",
        help="Output figure template for the per-snapshot gas dust-mass grid "
             "projection. Same {name}/{index} substitution as --grid-output. "
             "Skipped if the snapshot has no DustMass field.",
    )
    parser.add_argument(
        "--dtg-grid-output", default="disc_grid_dust_to_gas_{name}.png",
        help="Output figure template for the per-snapshot dust-to-gas mass "
             "ratio grid projection. Same {name}/{index} substitution as "
             "--grid-output. Skipped if the snapshot has no DustMass field.",
    )
    parser.add_argument(
        "--metals-dust-output", default="disc_metals_dust_profiles.png",
        help="Output figure for overlaid metal and dust surface-density radial "
             "profiles (mirrors --output for gas). Dust panel is omitted if no "
             "snapshot has a DustMass field.",
    )
    parser.add_argument(
        "--grid-vmin", type=float, default=None,
        help="Fixed lower colour limit (log10 scale) shared across all three "
             "panels of every grid projection (gas/metal/dust mass, and "
             "dust-to-gas ratio). Default: computed per-plot from that "
             "plot's own data range.",
    )
    parser.add_argument(
        "--grid-vmax", type=float, default=None,
        help="Fixed upper colour limit (log10 scale), see --grid-vmin.",
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


def compute_mask(pos, centre, args):
    """Boolean mask selecting particles within a cylindrical region around
    the centre (radius --mask-radius, half-height --mask-zheight). Returns
    an all-True mask if --mask-radius is unset."""
    if args.mask_radius is None:
        return np.ones(len(pos), dtype=bool)

    dxy = pos[:, :2] - centre[:2]
    r = np.hypot(dxy[:, 0], dxy[:, 1])
    dz = np.abs(pos[:, 2] - centre[2])

    return (r <= args.mask_radius) & (dz <= args.mask_zheight)


def apply_mask(pos, vel, mass, centre, args):
    """Restrict to a cylindrical region around the centre. Optional: improves
    speed and avoids off-disc contamination when enabled via --mask-radius."""
    mask = compute_mask(pos, centre, args)
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


def compute_metals_dust_profiles(pos, metal_mass, dust_mass, centre, args):
    """Radial surface-density profiles for gas metal mass and (if present)
    dust mass, mirroring compute_profiles' gas surface-density panel.
    pos/metal_mass/dust_mass are expected already masked (see apply_mask)."""
    profiler = ProfileTools(numbins=args.nbins)

    sigma_metals = profiler.surface_density(pos, metal_mass, centre, args.rmin, args.rmax)

    sigma_dust = None
    if dust_mass is not None:
        sigma_dust = profiler.surface_density(pos, dust_mass, centre, args.rmin, args.rmax)

    return sigma_metals, sigma_dust


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


def plot_metals_dust_profiles(results, output):
    """results: list of dicts with keys name, sigma_metals, sigma_dust
    (sigma_dust may be None for a given snapshot). Mirrors
    plot_overlaid_profiles' gas surface-density panel."""
    have_dust = any(res["sigma_dust"] is not None for res in results)

    ncols = 2 if have_dust else 1
    fig, axes = plt.subplots(1, ncols, figsize=(5 * ncols + 2, 4))
    axes = np.atleast_1d(axes)

    for res in results:
        axes[0].plot(res["sigma_metals"]["r"], res["sigma_metals"]["density"], label=res["name"])
        if have_dust and res["sigma_dust"] is not None:
            axes[1].plot(res["sigma_dust"]["r"], res["sigma_dust"]["density"], label=res["name"])

    axes[0].set_xlabel("R [kpc]")
    axes[0].set_ylabel(r"$\Sigma_{\rm metals}(R)$")
    axes[0].set_yscale("log")
    if len(results) > 1:
        axes[0].legend(fontsize="small")

    if have_dust:
        axes[1].set_xlabel("R [kpc]")
        axes[1].set_ylabel(r"$\Sigma_{\rm dust}(R)$")
        axes[1].set_yscale("log")

    fig.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)
    log.info("Saved %s", output)


def plot_grid(pos, values, centre, args, output_path, title=None, smoothing_lengths=None,
              value_scale=1.0, cbar_label=None, vmin=None, vmax=None):
    """Smoothed-grid projection plot for a single quantity/snapshot. Always
    run on the full particle set (regardless of --mask-radius), like the
    source notebook. `values` is whatever per-particle quantity is being
    projected (mass, metal mass, or dust mass), in code mass units.

    With `smoothing_lengths` (each cell's own effective radius, see
    read_gas_scalars_h5py), uses GriddingTools' per-particle SPH-kernel
    deposition instead of a single global Gaussian sigma applied to a
    CIC-binned grid -- a fixed sigma either over-smooths the dense inner
    disc or leaves particle-scale gaps in the sparser outer disc, since
    the AGORA IC's cell size varies by orders of magnitude across the box.
    Falls back to the old CIC+Gaussian path if smoothing_lengths is None
    (e.g. no Density field, or a particle type other than gas).

    Plots a proper column density (sum along the projection axis, divided
    by pixel area -- not a bare mean of raw code-unit values, which isn't
    a meaningful physical quantity). `value_scale` (see
    surface_density_value_scale) converts `values` from code mass units to
    Msun/pc^2 once area-normalized; pass 1.0 (default) to leave the
    colourbar in code units instead.
    """
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
    scaled_values = values[box_mask] * value_scale

    if smoothing_lengths is not None:
        smoothed_grid = gridding.smooth_to_grid(
            positions=pos[box_mask],
            values=scaled_values,
            grid_size=grid_size,
            grid_limits=grid_limits,
            method="SPH",
            smoothing_lengths=smoothing_lengths[box_mask],
        )
    else:
        smoothed_grid = gridding.smooth_to_grid(
            positions=pos[box_mask],
            values=scaled_values,
            grid_size=grid_size,
            grid_limits=grid_limits,
            method="Gaussian",
            sigma=args.smooth_sigma,
        )

    try:
        fig, _ = gridding.plot_3d_projections(
            smoothed_grid, grid_limits, projection="sum", normalize_by_area=True,
            title=title, vmin=vmin, vmax=vmax, cbar_label=cbar_label,
        )
    except TypeError:
        # the installed analysistools.GriddingTools.plot_3d_projections may
        # be an older version without vmin/vmax/normalize_by_area/cbar_label
        # -- degrade gracefully rather than failing the run
        fig, _ = gridding.plot_3d_projections(smoothed_grid, grid_limits, projection="sum", title=title)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    log.info("Saved %s", output_path)


def plot_dust_to_gas_ratio_grid(pos, gas_mass, dust_mass, smoothing_lengths, centre, args,
                                output_path, title=None, vmin=None, vmax=None):
    """Dust-to-gas mass ratio grid projection: gas mass and dust mass are
    each SPH-deposited and projected (summed) independently, then divided
    -- not divided cell-by-cell before deposition/projection, which would
    incorrectly weight the ratio by wherever gas happens to be sparse. No
    value_scale/area normalization needed: it cancels in the ratio (see
    GriddingTools.plot_ratio_projections), so this is dimensionless
    regardless of the mass unit gas_mass/dust_mass are given in.
    """
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

    gas_grid = gridding.smooth_to_grid(
        positions=pos[box_mask], values=gas_mass[box_mask], grid_size=grid_size,
        grid_limits=grid_limits, method="SPH", smoothing_lengths=smoothing_lengths[box_mask],
    )
    dust_grid = gridding.smooth_to_grid(
        positions=pos[box_mask], values=dust_mass[box_mask], grid_size=grid_size,
        grid_limits=grid_limits, method="SPH", smoothing_lengths=smoothing_lengths[box_mask],
    )

    try:
        fig, _ = gridding.plot_ratio_projections(
            dust_grid, gas_grid, grid_limits, projection="sum", title=title,
            vmin=vmin, vmax=vmax, cbar_label="log$_{10}$ Dust-to-gas ratio",
        )
    except TypeError:
        fig, _ = gridding.plot_ratio_projections(dust_grid, gas_grid, grid_limits, projection="sum", title=title)
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
        if res.get("sigma_metals") is not None:
            arrays[f"{name}__sigma_metals_r"] = res["sigma_metals"]["r"]
            arrays[f"{name}__sigma_metals_density"] = res["sigma_metals"]["density"]
        if res.get("sigma_dust") is not None:
            arrays[f"{name}__sigma_dust_r"] = res["sigma_dust"]["r"]
            arrays[f"{name}__sigma_dust_density"] = res["sigma_dust"]["density"]
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

        # Gas-only metal/dust scalars and per-cell smoothing length,
        # independent of --ptype (metals/dust/Density only exist on
        # PartType0). Read before the gas-mass grid below so its own
        # smoothing_lengths are available too. Reused for the grid
        # projections (unmasked, mirroring the gas mass grid) and the
        # radial surface density profiles (masked, mirroring the gas sigma
        # profile).
        try:
            gas_pos, gas_mass, metal_mass, dust_mass, gas_h = read_gas_scalars_h5py(snapshot_path)
        except (OSError, IOError, KeyError) as exc:
            log.error("Could not read gas scalars from '%s': %s", snapshot_path, exc)
            sys.exit(1)

        # gas_h aligns 1:1 with `pos` only when --ptype is gas (both read
        # PartType0/Coordinates in the same order); smoothing length has no
        # meaning for other particle types.
        mass_smoothing_lengths = gas_h if args.ptype == 0 else None

        value_scale = surface_density_value_scale(snapshot_path)
        sigma_label = r"$\Sigma$ [M$_\odot$ pc$^{-2}$]"

        grid_out = format_output_path(args.grid_output, name, index)
        plot_grid(pos, mass, centre, args, grid_out, title=f"Gas mass ({name})",
                  smoothing_lengths=mass_smoothing_lengths, value_scale=value_scale,
                  cbar_label=sigma_label, vmin=args.grid_vmin, vmax=args.grid_vmax)

        sigma_metals = sigma_dust = None

        if metal_mass is not None:
            metals_grid_out = format_output_path(args.metals_grid_output, name, index)
            plot_grid(gas_pos, metal_mass, centre, args, metals_grid_out, title=f"Metal mass ({name})",
                      smoothing_lengths=gas_h, value_scale=value_scale,
                      cbar_label=sigma_label, vmin=args.grid_vmin, vmax=args.grid_vmax)
        else:
            log.warning("No PassiveScalars field in %s -- skipping metals grid/profile", snapshot_path)

        if dust_mass is not None:
            dust_grid_out = format_output_path(args.dust_grid_output, name, index)
            plot_grid(gas_pos, dust_mass, centre, args, dust_grid_out, title=f"Dust mass ({name})",
                      smoothing_lengths=gas_h, value_scale=value_scale,
                      cbar_label=sigma_label, vmin=args.grid_vmin, vmax=args.grid_vmax)

            dtg_grid_out = format_output_path(args.dtg_grid_output, name, index)
            plot_dust_to_gas_ratio_grid(
                gas_pos, gas_mass, dust_mass, gas_h, centre, args, dtg_grid_out,
                title=f"Dust-to-gas ratio ({name})", vmin=args.grid_vmin, vmax=args.grid_vmax,
            )
        else:
            log.info("No DustMass field in %s -- skipping dust/dust-to-gas grid/profile (run not built with DUST?)", snapshot_path)

        if metal_mass is not None:
            gas_mask = compute_mask(gas_pos, centre, args)
            gas_pos_m = gas_pos[gas_mask]
            metal_mass_m = metal_mass[gas_mask]
            dust_mass_m = dust_mass[gas_mask] if dust_mass is not None else None
            sigma_metals, sigma_dust = compute_metals_dust_profiles(
                gas_pos_m, metal_mass_m, dust_mass_m, centre, args
            )

        pos_m, vel_m, mass_m = apply_mask(pos, vel, mass, centre, args)
        sigma, hz, sigmaz = compute_profiles(pos_m, vel_m, mass_m, centre, args)

        results.append({
            "name": name, "sigma": sigma, "hz": hz, "sigmaz": sigmaz,
            "sigma_metals": sigma_metals, "sigma_dust": sigma_dust,
        })

    plot_overlaid_profiles(results, args.output)

    if any(res["sigma_metals"] is not None for res in results):
        plot_metals_dust_profiles(results, args.metals_dust_output)
    else:
        log.warning("No metals data in any snapshot -- skipping %s", args.metals_dust_output)

    if args.save_data:
        save_data(results, args.save_data)


if __name__ == "__main__":
    main()
