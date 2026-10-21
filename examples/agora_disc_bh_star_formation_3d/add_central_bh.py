#!/usr/bin/env python3
"""Add one central black hole (PartType5) to an Arepo HDF5 initial-conditions file.

The AGORA disc ICs contain gas, halo, disc and bulge only. This copies the input
file and appends a single PartType5 particle at the box centre (the AGORA discs are
constructed about BoxSize/2), at rest in the construction frame, with the next unused
ParticleID. NumPart_ThisFile / NumPart_Total in /Header are updated to match.

    python add_central_bh.py IN.hdf5 OUT.hdf5 [--mass-msun 1e6] [--unit-mass-g G] [--set-gas-u U]

--unit-mass-g   overrides UnitMass_in_g, both when converting the BH mass and in the written header.
                (The AGORA medium-res IC declares 1.989e43 g although its masses are in 1e9 Msun.)
--set-gas-u     sets InternalEnergy of every gas cell with non-zero mass to U (code units, (km/s)^2),
                used to give the medium-res IC the same thermal state as the low-res one.
"""
import argparse
import shutil

import h5py
import numpy as np

MSUN_G = 1.989e33


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("infile")
    ap.add_argument("outfile")
    ap.add_argument("--mass-msun", type=float, default=1.0e6, help="BH mass in solar masses (default 1e6)")
    ap.add_argument("--unit-mass-g", type=float, default=None)
    ap.add_argument("--set-gas-u", type=float, default=None)
    a = ap.parse_args()

    shutil.copyfile(a.infile, a.outfile)
    with h5py.File(a.outfile, "r+") as f:
        hd = f["Header"].attrs
        if "PartType5" in f:
            raise SystemExit("PartType5 already present in %s; refusing to add a second black hole" % a.infile)

        if a.unit_mass_g is not None:
            print("UnitMass_in_g: header %g -> %g" % (hd["UnitMass_in_g"], a.unit_mass_g))
            hd["UnitMass_in_g"] = a.unit_mass_g
        mass = a.mass_msun * MSUN_G / hd["UnitMass_in_g"]  # code mass units
        centre = 0.5 * hd["BoxSize"]
        types = [t for t in range(5) if "PartType%d" % t in f]
        next_id = max(int(f["PartType%d/ParticleIDs" % t][:].max()) for t in types) + 1

        pos_dtype = f["PartType1/Coordinates"].dtype
        id_dtype = f["PartType1/ParticleIDs"].dtype
        if next_id > np.iinfo(id_dtype).max:
            raise SystemExit("next ParticleID %d overflows %s" % (next_id, id_dtype))

        g = f.create_group("PartType5")
        g.create_dataset("Coordinates", data=np.full((1, 3), centre, dtype=pos_dtype))
        g.create_dataset("Velocities", data=np.zeros((1, 3), dtype=pos_dtype))
        g.create_dataset("Masses", data=np.array([mass], dtype=f["PartType1/Masses"].dtype))
        g.create_dataset("ParticleIDs", data=np.array([next_id], dtype=id_dtype))

        for key in ("NumPart_ThisFile", "NumPart_Total"):
            n = hd[key].copy()
            n[5] += 1
            hd[key] = n

        if a.set_gas_u is not None:
            m0 = f["PartType0/Masses"][:]
            u = f["PartType0/InternalEnergy"][:]
            u[m0 > 0] = a.set_gas_u
            f["PartType0/InternalEnergy"][...] = u
            print("gas InternalEnergy set to %g in %d cells" % (a.set_gas_u, (m0 > 0).sum()))

        # diagnostics: where is the disc centre of mass relative to the BH?
        # stellar disc: PartType2 in the low-res IC, PartType1 in the medium-res IC (the halo is the other)
        mt = [f["PartType%d/Masses" % t][0] for t in (1, 2)]
        disc = "PartType1" if mt[0] < mt[1] else "PartType2"
        m = f[disc + "/Masses"][:].astype(np.float64)
        x = f[disc + "/Coordinates"][:].astype(np.float64)
        com = (x * m[:, None]).sum(0) / m.sum()
        print("BH: ID=%d mass=%.6g code units (%.3g Msun) at (%g, %g, %g), v=0" % (next_id, mass, a.mass_msun, centre, centre, centre))
        print("stellar-disc (%s) COM offset from BH [code length]: %s" % (disc, np.array2string(com - centre, precision=3)))


if __name__ == "__main__":
    main()
