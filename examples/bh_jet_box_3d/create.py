""" @package ./examples/bh_jet_box_3d/create.py
Creates the initial conditions for the BH jet box: a uniform hot-gas medium in a periodic box that is
BoxSize x BoxSize x 2*BoxSize, with one black hole at its centre.

usage: python3 create.py <simulation_directory> [cells_per_side_in_x_and_y]

Code units (match param.txt): 1 kpc, 1 km/s, 1e6 Msun.
"""

import sys
import numpy as np
import h5py

simulation_directory = str(sys.argv[1])
CellsPerSide = int(sys.argv[2]) if len(sys.argv) > 2 else 32

print("examples/bh_jet_box_3d/create.py: creating ICs in directory " + simulation_directory)

FloatType = np.float64
IntType = np.int32

""" physical constants and units; must agree with UnitLength/Mass/Velocity in param.txt """
UnitLength_in_cm = 3.086e21
UnitMass_in_g = 1.989e39
UnitVelocity_in_cm_per_s = 1.0e5
PROTONMASS = 1.6726e-24
BOLTZMANN = 1.3806e-16
SOLAR_MASS = 1.989e33
GAMMA = 5.0 / 3.0

""" box: BoxSize (x, y) and 2*BoxSize (z, from LONG_Z=2 in Config.sh) """
BoxSize = FloatType(64.0)
LONG_Z = 2
Lx, Ly, Lz = BoxSize, BoxSize, BoxSize * LONG_Z

""" gas: uniform, at rest, hot (a stand-in for a group/cluster atmosphere; there is no gravity to hold a profile) """
n_H = 1.0e-2             # hydrogen number density [cm^-3]
X_H = 0.76               # hydrogen mass fraction
T = 1.0e7                # temperature [K]
mu = 0.6                 # mean molecular weight

rho_cgs = n_H * PROTONMASS / X_H
rho = rho_cgs * UnitLength_in_cm**3 / UnitMass_in_g                                       # code density
utherm = BOLTZMANN * T / (mu * PROTONMASS) / (GAMMA - 1.0) / UnitVelocity_in_cm_per_s**2  # code specific energy

""" grid: CellsPerSide in x and y, LONG_Z times as many in z, so the cells are cubic """
Nx, Ny, Nz = CellsPerSide, CellsPerSide, CellsPerSide * LONG_Z
NumberOfCells = Nx * Ny * Nz
dx = Lx / Nx

x1 = (np.arange(Nx) + 0.5) * dx
y1 = (np.arange(Ny) + 0.5) * dx
z1 = (np.arange(Nz) + 0.5) * dx
xx, yy, zz = np.meshgrid(x1, y1, z1, indexing="ij")
Pos = np.zeros([NumberOfCells, 3], dtype=FloatType)
Pos[:, 0] = xx.reshape(NumberOfCells)
Pos[:, 1] = yy.reshape(NumberOfCells)
Pos[:, 2] = zz.reshape(NumberOfCells)

Mass = np.full(NumberOfCells, rho * dx**3, dtype=FloatType)
Velocity = np.zeros([NumberOfCells, 3], dtype=FloatType)
Uthermal = np.full(NumberOfCells, utherm, dtype=FloatType)

""" black hole at the centre of the box (the jet code launches along +-z from here).
    Its mass only matters for Bondi accretion, which BH_JET_FIXED_POWER ignores. """
BhMass = FloatType(1.0)            # 1e6 Msun
BhPos = np.array([[0.5 * Lx, 0.5 * Ly, 0.5 * Lz]], dtype=FloatType)
BhVel = np.zeros([1, 3], dtype=FloatType)

""" write *.hdf5 file """
IC = h5py.File(simulation_directory + "/IC.hdf5", "w")

header = IC.create_group("Header")
part0 = IC.create_group("PartType0")
part5 = IC.create_group("PartType5")

NumPart = np.array([NumberOfCells, 0, 0, 0, 0, 1], dtype=IntType)
header.attrs.create("NumPart_ThisFile", NumPart)
header.attrs.create("NumPart_Total", NumPart)
header.attrs.create("NumPart_Total_HighWord", np.zeros(6, dtype=IntType))
header.attrs.create("MassTable", np.zeros(6, dtype=FloatType))
header.attrs.create("Time", 0.0)
header.attrs.create("Redshift", 0.0)
header.attrs.create("BoxSize", Lx)   # Arepo takes the base box size; the z extent comes from LONG_Z
header.attrs.create("NumFilesPerSnapshot", 1)
header.attrs.create("Omega0", 0.0)
header.attrs.create("OmegaB", 0.0)
header.attrs.create("OmegaLambda", 0.0)
header.attrs.create("HubbleParam", 1.0)
header.attrs.create("Flag_Sfr", 0)
header.attrs.create("Flag_Cooling", 0)
header.attrs.create("Flag_StellarAge", 0)
header.attrs.create("Flag_Metals", 0)
header.attrs.create("Flag_Feedback", 0)
header.attrs.create("Flag_DoublePrecision", 1)

part0.create_dataset("ParticleIDs", data=np.arange(1, NumberOfCells + 1, dtype=np.int64))
part0.create_dataset("Coordinates", data=Pos)
part0.create_dataset("Masses", data=Mass)
part0.create_dataset("Velocities", data=Velocity)
part0.create_dataset("InternalEnergy", data=Uthermal)

part5.create_dataset("ParticleIDs", data=np.array([NumberOfCells + 1], dtype=np.int64))
part5.create_dataset("Coordinates", data=BhPos)
part5.create_dataset("Masses", data=np.array([BhMass], dtype=FloatType))
part5.create_dataset("Velocities", data=BhVel)

IC.close()

""" summary, and the values in param.txt that depend on the grid """
cellmass_msun = rho * dx**3 * UnitMass_in_g / SOLAR_MASS
print("  box                : %g x %g x %g kpc, %d x %d x %d = %d gas cells (dx = %g kpc) + 1 BH" % (Lx, Ly, Lz, Nx, Ny, Nz, NumberOfCells, dx))
print("  gas                : n_H = %g cm^-3, T = %g K, rho = %g, u = %g (code units)" % (n_H, T, rho, utherm))
print("  cell mass          : %g code = %g Msun  -> ReferenceGasPartMass in param.txt should be %.3g" % (rho * dx**3, cellmass_msun, rho * dx**3))
print("  total gas energy   : %g code" % (NumberOfCells * rho * dx**3 * utherm))
