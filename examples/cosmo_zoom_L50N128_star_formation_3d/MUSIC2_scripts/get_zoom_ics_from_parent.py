#!/usr/bin/env python3

"""
Identify a halo, trace its particles back to the initial conditions,
and determine the corresponding Lagrangian patch.

Can be used from the command line:

    python get_zoom_ics_from_parent.py \
        --snapshot snap_020.hdf5 \
        --halocat fof_subhalo_tab_020.hdf5 \
        --ics L50_N128_z99.hdf5 \
        --mass-min 95 \
        --mass-max 105

or imported and called from a notebook:

    import get_zoom_ics_from_parent as halo_zoom

    result = halo_zoom.main([
        "--snapshot", "snap_020.hdf5",
        "--halocat", "fof_subhalo_tab_020.hdf5",
        "--ics", "L50_N128_z99.hdf5",
        "--mass-min", "95",
        "--mass-max", "105",
    ])

The returned dictionary contains the target halo properties, selected
particle indices/IDs, and Lagrangian patch centre and extent.

CBP, 10/09/2026
"""

import argparse
from pathlib import Path

import analysistools as at
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

def get_halo_properties(halos, mass_range):
    """
    Extract halo properties for halos in the requested mass range.

    Parameters
    ----------
    halos
        analysistools halo catalogue.
    mass_range : tuple
        (minimum, maximum) halo mass. Either limit may be None.

    Returns
    -------
    tuple
        masses, radii, positions, halo IDs
    """
    selected = halos.select(mass=mass_range)

    return (
        selected["mass"],
        selected["radius"],
        selected["pos"],
        selected["halo_id"],
    )

def select_target_halo(halos, mass_range, index=0):
    """
    Select one halo from the requested mass range.

    Parameters
    ----------
    halos
        analysistools halo catalogue.
    mass_range : tuple
        Halo mass range.
    index : int
        Index of the selected halo within the mass-selected catalogue.

    Returns
    -------
    dict
        Halo properties.
    """
    mass, radius, pos, halo_id = get_halo_properties(halos, mass_range)

    if len(mass) == 0:
        raise ValueError(
            f"No halos found in mass range {mass_range}."
        )

    if index >= len(mass):
        raise IndexError(
            f"Halo index {index} requested, but only {len(mass)} "
            f"halos are available in the selected mass range."
        )

    return {
        "mass": mass[index],
        "radius": radius[index],
        "pos": pos[index],
        "id": halo_id[index],
    }

def minimum_image_displacement(pos, centre, boxsize):
    """
    Return displacement from centre using the minimum-image convention.
    """
    dpos = pos - centre
    dpos -= boxsize * np.round(dpos / boxsize)
    return dpos


def unwrap_positions(pos, centre, boxsize):
    """
    Put positions into the periodic image nearest to centre.
    """
    return centre + minimum_image_displacement(pos, centre, boxsize)

def get_halo_particle_mask(snap, halo_props, fcut=5.0, spherical=True):
    """
    Identify particles within fcut times the halo radius.

    Periodic boundaries are handled using the minimum-image convention.

    Parameters
    ----------
    snap
        analysistools snapshot.
    halo_props : dict
        Dictionary containing ``pos`` and ``radius``.
    fcut : float
        Radius of the region in units of the halo radius.
    spherical : bool
        If True, select a sphere. If False, select a cube.

    Returns
    -------
    ndarray
        Boolean mask selecting particles in the region.
    """
    dpos = minimum_image_displacement(
        snap.pos,
        halo_props["pos"],
        snap.meta["boxsize"],
    )

    radius = fcut * halo_props["radius"]

    if spherical:
        mask = np.sum(dpos**2, axis=1) < radius**2
    else:
        mask = np.all(np.abs(dpos) < radius, axis=1)

    return mask

def get_lagrangian_patch(ics_snap, particle_ids, f_padding=1.25):
    """
    Determine the bounding box of the selected particles in the ICs.

    The IC coordinates are unwrapped before calculating the bounding box,
    so that patches crossing a periodic boundary are treated correctly.

    Parameters
    ----------
    ics_snap
        analysistools IC snapshot.
    particle_ids
        Particle IDs corresponding to the selected halo particles.
    f_padding : float
        Multiplicative padding applied to the patch dimensions.

    Returns
    -------
    centre : ndarray
        Centre of the Lagrangian patch in box units.
    extent : ndarray
        Patch dimensions in box units.
    positions : ndarray
        Unwrapped IC particle positions in simulation units.
    """
    positions = ics_snap.pos[particle_ids]
    boxsize = ics_snap.meta["boxsize"]

    # Obtain a reference point and unwrap all particles around it.
    reference = positions[0]

    positions = unwrap_positions(
        positions,
        reference,
        boxsize,
    )

    minimum = positions.min(axis=0)
    maximum = positions.max(axis=0)

    centre = 0.5 * (minimum + maximum)
    extent = f_padding * (maximum - minimum)

    centre_box = centre / boxsize
    extent_box = extent / boxsize

    return centre_box, extent_box, positions

def overplot_radii(
    ax,
    halos,
    mass_range,
    projection="xy",
    color="white",
):
    """
    Overplot virial-radius circles for halos in a mass range.
    """
    _, radii, positions, _ = get_halo_properties(
        halos,
        mass_range,
    )

    projection_mapping = {
        "xy": (0, 1),
        "yx": (1, 0),
        "yz": (1, 2),
        "zy": (2, 1),
        "xz": (0, 2),
        "zx": (2, 0),
    }

    projection = projection.lower()

    if projection not in projection_mapping:
        raise ValueError(
            "Projection must be one of: "
            "'xy', 'yx', 'yz', 'zy', 'xz', 'zx'."
        )

    idx_x, idx_y = projection_mapping[projection]

    for radius, pos in zip(radii, positions):
        ax.add_patch(
            Circle(
                (pos[idx_x], pos[idx_y]),
                radius,
                fill=False,
                color=color,
            )
        )

    return ax

def plot_box_projection(
    snap,
    halos,
    mass_range,
    output=None,
    gridsize=512,
):
    """
    Plot full simulation-box projections with halo radii overplotted.
    """
    fig, axes = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=(16, 8),
    )

    # x-y
    axes[0].hexbin(
        snap.pos[:, 0],
        snap.pos[:, 1],
        gridsize=gridsize,
        bins="log",
        cmap="plasma",
    )
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    axes[0].set_aspect("equal")

    overplot_radii(
        axes[0],
        halos,
        mass_range,
        projection="xy",
    )

    # x-z
    axes[1].hexbin(
        snap.pos[:, 0],
        snap.pos[:, 2],
        gridsize=gridsize,
        bins="log",
        cmap="plasma",
    )
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("z")
    axes[1].set_aspect("equal")

    overplot_radii(
        axes[1],
        halos,
        mass_range,
        projection="xz",
    )

    fig.tight_layout()

    if output is not None:
        fig.savefig(output, dpi=300)

    return fig, axes

def plot_halo_projection(
    snap,
    ics_snap,
    particle_mask,
    halo_props,
    output=None,
):
    """
    Plot the selected halo particles at z=0 and in the ICs.
    """
    particle_indices = np.flatnonzero(particle_mask)
    particle_ids = snap.pids[particle_indices]

    # Evolved positions: unwrap relative to halo centre.
    halo_positions = unwrap_positions(
        snap.pos[particle_indices],
        halo_props["pos"],
        snap.meta["boxsize"],
    )

    # IC positions: unwrap relative to first selected particle.
    ic_positions = ics_snap.pos[particle_ids]
    ic_positions = unwrap_positions(
        ic_positions,
        ic_positions[0],
        ics_snap.meta["boxsize"],
    )

    fig, axes = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(16, 16),
    )

    halo_id = halo_props["id"]

    # Evolved x-y
    axes[0, 0].scatter(
        halo_positions[:, 0],
        halo_positions[:, 1],
        s=1,
    )
    axes[0, 0].set_title(
        f"Halo {halo_id} - z={snap.meta['redshift']}"
    )
    axes[0, 0].set_xlabel("x")
    axes[0, 0].set_ylabel("y")
#    axes[0, 0].set_aspect("equal")

    # Evolved x-z
    axes[0, 1].scatter(
        halo_positions[:, 0],
        halo_positions[:, 2],
        s=1,
    )
    axes[0, 1].set_title(
        f"Halo {halo_id} - z={snap.meta['redshift']}"
    )
    axes[0, 1].set_xlabel("x")
    axes[0, 1].set_ylabel("z")
#    axes[0, 1].set_aspect("equal")

    # IC x-y
    axes[1, 0].scatter(
        ic_positions[:, 0],
        ic_positions[:, 1],
        s=1,
    )
    axes[1, 0].set_title(
        f"Halo {halo_id} - z={ics_snap.meta['redshift']}"
    )
    axes[1, 0].set_xlabel("x")
    axes[1, 0].set_ylabel("y")
#    axes[1, 0].set_aspect("equal")

    # IC x-z
    axes[1, 1].scatter(
        ic_positions[:, 0],
        ic_positions[:, 2],
        s=1,
    )
    axes[1, 1].set_title(
        f"Halo {halo_id} - z={ics_snap.meta['redshift']}"
    )
    axes[1, 1].set_xlabel("x")
    axes[1, 1].set_ylabel("z")
#    axes[1, 1].set_aspect("equal")

    fig.tight_layout()

    if output is not None:
        fig.savefig(output, dpi=300)

    return fig, axes

def analyse(
    snapshot_filename,
    halocat_filename,
    ics_snapshot_filename,
    mass_range=(95, 105),
    halo_index=0,
    fcut=5.0,
    f_padding=1.25,
    spherical=True,
    output_dir=".",
    make_plots=True,
):
    """
    Run the complete halo-to-Lagrangian-patch analysis.

    This is the main programmatic interface for notebooks and scripts.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading particle data from {snapshot_filename}")
    snap = at.load(snapshot_filename).preload()

    print(f"Redshift: {snap.meta['redshift']}")
    print(f"BoxSize: {snap.meta['boxsize']}")

    print(f"Reading catalogue data from {halocat_filename}")
    halos = at.load(halocat_filename)

    print(f"Reading ICs particle data from {ics_snapshot_filename}")
    ics_snap = at.load(ics_snapshot_filename).preload()

    halo_props = select_target_halo(
        halos,
        mass_range,
        index=halo_index,
    )

    print()
    print("Target halo")
    print("-----------")
    print(f"ID:     {halo_props['id']}")
    print(f"Mass:   {halo_props['mass']}")
    print(f"Radius: {halo_props['radius']}")
    print(f"Pos:    {halo_props['pos']}")

    particle_mask = get_halo_particle_mask(
        snap,
        halo_props,
        fcut=fcut,
        spherical=spherical,
    )

    particle_indices = np.flatnonzero(particle_mask)
    particle_ids = snap.pids[particle_indices]

    print()
    print("Particle selection")
    print("------------------")
    print(f"Selection radius: {fcut} Rvir")
    print(f"Geometry:         {'sphere' if spherical else 'cube'}")
    print(f"Particles:        {len(particle_ids)}")

    patch_centre, patch_extent, ic_positions = get_lagrangian_patch(
        ics_snap,
        particle_ids,
        f_padding=f_padding,
    )

    print()
    print("Lagrangian patch")
    print("----------------")
    print(f"Centre (box units): {patch_centre}")
    print(f"Extent (box units): {patch_extent}")

    box_output = (
        output_dir
        / "Box_projection.png"
    )

    halo_output = (
        output_dir
        / f"Halo_projection_{halo_props['id']}.png"
    )

    if make_plots:
        plot_box_projection(
            snap,
            halos,
            mass_range,
            output=box_output,
        )

        plot_halo_projection(
            snap,
            ics_snap,
            particle_mask,
            halo_props,
            output=halo_output,
        )

    return {
        "halo": halo_props,
        "particle_mask": particle_mask,
        "particle_indices": particle_indices,
        "particle_ids": particle_ids,
        "patch_centre": patch_centre,
        "patch_extent": patch_extent,
        "ic_positions": ic_positions,
        "box_output": box_output if make_plots else None,
        "halo_output": halo_output if make_plots else None,
    }

def build_parser():
    """
    Construct the command-line argument parser.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Identify a halo and determine its Lagrangian patch "
            "in the initial conditions."
        )
    )

    parser.add_argument(
        "--snapshot",
        required=True,
        help="Evolved snapshot filename.",
    )

    parser.add_argument(
        "--halocat",
        required=True,
        help="Halo catalogue filename.",
    )

    parser.add_argument(
        "--ics",
        required=True,
        help="Initial conditions snapshot filename.",
    )

    parser.add_argument(
        "--mass-min",
        type=float,
        required=True,
        help="Minimum halo mass.",
    )

    parser.add_argument(
        "--mass-max",
        type=float,
        default=None,
        help="Maximum halo mass.",
    )

    parser.add_argument(
        "--halo-index",
        type=int,
        default=0,
        help="Index of halo within selected mass range.",
    )

    parser.add_argument(
        "--fcut",
        type=float,
        default=5.0,
        help="Particle selection radius in units of Rvir.",
    )

    parser.add_argument(
        "--padding",
        type=float,
        default=1.25,
        help="Padding factor for the Lagrangian patch.",
    )

    parser.add_argument(
        "--cube",
        action="store_true",
        help="Use a cubic rather than spherical particle selection.",
    )

    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory for output plots.",
    )

    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Do not generate plots.",
    )

    return parser


def main(argv=None):
    """
    Command-line entry point.

    ``argv=None`` means parse arguments from sys.argv.

    Supplying a list makes the function convenient to call from a notebook.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    return analyse(
        snapshot_filename=args.snapshot,
        halocat_filename=args.halocat,
        ics_snapshot_filename=args.ics,
        mass_range=(args.mass_min, args.mass_max),
        halo_index=args.halo_index,
        fcut=args.fcut,
        f_padding=args.padding,
        spherical=not args.cube,
        output_dir=args.output_dir,
        make_plots=not args.no_plots,
    )

if __name__ == "__main__":
    main()
