# Modules for OzSTAR/NT (from the AGORA_LOWRES_TEST job script). Sourced by the .sbatch files.
module load gcc/11.2.0 openmpi/4.1.1 hdf5/1.12.1 fftw/3.3.10 gsl/2.7

# Grackle must be findable by the Makefile: set GRACKLE_DIR in Makefile.systype (repository root), or point to an install here.
# export GRACKLE_DIR=/fred/oz009/cpower/software/grackle_solas

# compare_snapshots.py needs python with numpy + h5py. EDIT: a module providing both, or a venv to activate.
# module load python/<version> scipy-bundle/<version> h5py/<version>
# source /fred/oz009/cpower/venv/bin/activate
