""" @package ./examples/bh_jet_stratified_3d/create.py
Creates the initial conditions for the BH jet in a stratified atmosphere: an isothermal beta-model gas distribution in a periodic box that is
BoxSize x BoxSize x 2*BoxSize, with one black hole at its centre.

usage: python3 create.py <simulation_directory> [cells_per_side_in_x_and_y] [nobh]

  The density is eq. 13 of Maragkakis et al. (2026): rho = rho_0 [1 + (r/r_c)^2]^(-3 beta / 2) with rho_0 = 1e-25 g/cm^3, beta = 0.38, r_c = 10 kpc, at
  T = 1e7 K. The cells have equal mass and are drawn from the profile by inverse-transform sampling inside the sphere that encloses the box, keeping the
  cells that fall in the box; their number is about cells_per_side x cells_per_side x 2*cells_per_side (default 32). The gas is in hydrostatic equilibrium
  only with the acceleration of eq. 14: build with STATICBETA (see Config.sh).
  nobh               leave out the black hole (for the relaxation run; add_bh.py puts it back)

Code units (match param.txt): 1 kpc, 1 km/s, 1e6 Msun.
"""

import sys
import numpy as np
import h5py

simulation_directory = str(sys.argv[1])
CellsPerSide = int(sys.argv[2]) if len(sys.argv) > 2 else 32
with_bh = not (len(sys.argv) > 3 and sys.argv[3] == "nobh")

print("examples/bh_jet_stratified_3d/create.py: creating ICs in directory " + simulation_directory)

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

""" gas, at rest, T = 1e7 K, mu = 0.6 """
T = 1.0e7                # temperature [K]
mu = 0.6                 # mean molecular weight
utherm = BOLTZMANN * T / (mu * PROTONMASS) / (GAMMA - 1.0) / UnitVelocity_in_cm_per_s**2  # code specific energy

Nx, Ny, Nz = CellsPerSide, CellsPerSide, CellsPerSide * LONG_Z   # sets the number of cells: about Nx * Ny * Nz
centre = np.array([0.5 * Lx, 0.5 * Ly, 0.5 * Lz])

""" isothermal beta model, eq. 13: rho = rho_0 [1 + (r/r_c)^2]^(-3 beta / 2), rho_0 = 1e-25 g/cm^3, beta = 0.38, r_c = 10 kpc.
    Equal-mass cells: draw radii from the enclosed-mass profile of the sphere that encloses the box, with isotropic directions,
    and keep the cells that lie in the box. The number kept is about Nx*Ny*Nz (the cell mass is the mass in the box over that). """
rho0_cgs = 1.0e-25
beta, rc = 0.38, 10.0                                                       # rc in kpc (= code length)
rho0 = rho0_cgs * UnitLength_in_cm**3 / UnitMass_in_g
rng = np.random.default_rng(12345)
Rs = 0.5 * np.sqrt(Lx**2 + Ly**2 + Lz**2)                                   # sphere enclosing the box

r_tab = np.linspace(0.0, Rs, 40001)
shell = 4.0 * np.pi * r_tab**2 * rho0 * (1.0 + (r_tab / rc) ** 2) ** (-1.5 * beta)
M_tab = np.concatenate(([0.0], np.cumsum(0.5 * (shell[1:] + shell[:-1]) * np.diff(r_tab))))   # M(<r)

nf = 4 * CellsPerSide                                                       # fine grid for the mass inside the box
fx = (np.arange(nf) + 0.5) * Lx / nf
fz = (np.arange(nf * LONG_Z) + 0.5) * Lx / nf
gx, gy, gz = np.meshgrid(fx - centre[0], fx - centre[1], fz - centre[2], indexing="ij")
M_box = np.sum(rho0 * (1.0 + (gx**2 + gy**2 + gz**2) / rc**2) ** (-1.5 * beta)) * (Lx / nf) ** 3
cellmass = M_box / (Nx * Ny * Nz)

n_sphere = int(round(M_tab[-1] / cellmass))
r = np.interp(rng.random(n_sphere) * M_tab[-1], M_tab, r_tab)
d = rng.normal(size=(n_sphere, 3))
d /= np.linalg.norm(d, axis=1)[:, None]
P_all = centre + r[:, None] * d
inbox = (P_all[:, 0] >= 0) & (P_all[:, 0] < Lx) & (P_all[:, 1] >= 0) & (P_all[:, 1] < Ly) & (P_all[:, 2] >= 0) & (P_all[:, 2] < Lz)
Pos = P_all[inbox].astype(FloatType)
NumberOfCells = len(Pos)
Mass = np.full(NumberOfCells, cellmass, dtype=FloatType)
rho = rho0
summary = "beta model, rho_0 = %g g/cm^3, beta = %g, r_c = %g kpc" % (rho0_cgs, beta, rc)

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
if with_bh:
    part5 = IC.create_group("PartType5")

NumPart = np.array([NumberOfCells, 0, 0, 0, 0, 1 if with_bh else 0], dtype=IntType)
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

if with_bh:
    part5.create_dataset("ParticleIDs", data=np.array([NumberOfCells + 1], dtype=np.int64))
    part5.create_dataset("Coordinates", data=BhPos)
    part5.create_dataset("Masses", data=np.array([BhMass], dtype=FloatType))
    part5.create_dataset("Velocities", data=BhVel)

IC.close()

""" summary, and the values in param.txt that depend on the grid """
print("  box                : %g x %g x %g kpc, %d gas cells%s" % (Lx, Ly, Lz, NumberOfCells, " + 1 BH" if with_bh else ""))
print("  gas                : %s, T = %g K, u = %g (code units)" % (summary, T, utherm))
print("  cell mass          : %g code = %g Msun  -> ReferenceGasPartMass in param.txt should be %.3g" % (cellmass, cellmass * UnitMass_in_g / SOLAR_MASS, cellmass))
print("  total gas energy   : %g code" % (NumberOfCells * cellmass * utherm))
