""" @package ./examples/bh_jet_box_3d/add_bh.py
Turns a snapshot of the relaxed, black-hole-free atmosphere into the initial conditions of the jet run by adding the black hole at the
centre of the box (the relaxation step of Maragkakis et al. 2026, section 4.5).

usage: python3 add_bh.py <relaxed_snapshot.hdf5> <new_IC.hdf5> [bh_mass_in_code_units]

The gas is copied as it is (positions, velocities, masses, internal energies, IDs), including any motions left over from the relaxation.
The header time is reset to 0, so the jet run starts at TimeBegin = 0; the BH gets the next free ID and is at rest at the box centre
(BoxSize/2, BoxSize/2, LONG_Z*BoxSize/2). LONG_Z must match Config.sh (2 in the example). The BH mass defaults to 1 (1e6 Msun).
"""

import sys
import numpy as np
import h5py

LONG_Z = 2

src, dst = sys.argv[1], sys.argv[2]
bh_mass = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0

FloatType = np.float64
IntType = np.int32

snap = h5py.File(src, "r")
head = snap["Header"].attrs
box = float(head["BoxSize"])
gas = snap["PartType0"]

ids = gas["ParticleIDs"][:]
if "PartType5" in snap and len(snap["PartType5"]["ParticleIDs"]) > 0:
    raise SystemExit("the snapshot already contains a black hole")
n = len(ids)

print("examples/bh_jet_box_3d/add_bh.py: %s (%d gas cells, t = %g) -> %s" % (src, n, head["Time"], dst))

IC = h5py.File(dst, "w")
header = IC.create_group("Header")
header.attrs.create("NumPart_ThisFile", np.array([n, 0, 0, 0, 0, 1], dtype=IntType))
header.attrs.create("NumPart_Total", np.array([n, 0, 0, 0, 0, 1], dtype=IntType))
header.attrs.create("NumPart_Total_HighWord", np.zeros(6, dtype=IntType))
header.attrs.create("MassTable", np.zeros(6, dtype=FloatType))
header.attrs.create("Time", 0.0)
header.attrs.create("Redshift", 0.0)
header.attrs.create("BoxSize", box)
header.attrs.create("NumFilesPerSnapshot", 1)
for k, v in (("Omega0", 0.0), ("OmegaB", 0.0), ("OmegaLambda", 0.0), ("HubbleParam", 1.0)):
    header.attrs.create(k, float(head[k]) if k in head else v)
for k in ("Flag_Sfr", "Flag_Cooling", "Flag_StellarAge", "Flag_Metals", "Flag_Feedback"):
    header.attrs.create(k, 0)
header.attrs.create("Flag_DoublePrecision", 1)

p0 = IC.create_group("PartType0")
p0.create_dataset("ParticleIDs", data=ids.astype(np.int64))
for name in ("Coordinates", "Masses", "Velocities", "InternalEnergy"):
    p0.create_dataset(name, data=gas[name][:].astype(FloatType))

p5 = IC.create_group("PartType5")
p5.create_dataset("ParticleIDs", data=np.array([int(ids.max()) + 1], dtype=np.int64))
p5.create_dataset("Coordinates", data=np.array([[0.5 * box, 0.5 * box, 0.5 * box * LONG_Z]], dtype=FloatType))
p5.create_dataset("Masses", data=np.array([bh_mass], dtype=FloatType))
p5.create_dataset("Velocities", data=np.zeros([1, 3], dtype=FloatType))

IC.close()
snap.close()
print("  added a BH of mass %g at (%g, %g, %g)" % (bh_mass, 0.5 * box, 0.5 * box, 0.5 * box * LONG_Z))
