#include <stdio.h>
#include "arepoconfig.h"
#ifdef HAVE_HDF5
#include <hdf5.h>

hid_t my_H5Acreate(hid_t loc_id, const char *attr_name, hid_t type_id, hid_t space_id, hid_t acpl_id);
hid_t my_H5Screate(H5S_class_t type);
herr_t my_H5Aclose(hid_t attr_id, const char *attr_name);
herr_t my_H5Awrite(hid_t attr_id, hid_t mem_type_id, const void *buf, const char *attr_name);
herr_t my_H5Sclose(hid_t dataspace_id, H5S_class_t type);

herr_t my_H5Tclose(hid_t type_id);

void write_compile_time_options_in_hdf5(hid_t handle)
{
hid_t hdf5_dataspace, hdf5_attribute;
double val;
hid_t atype = H5Tcopy(H5T_C_S1);
H5Tset_size(atype, 1);
hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "METALS" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "METALS");
my_H5Aclose(hdf5_attribute, "METALS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "USE_GRACKLE" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "USE_GRACKLE");
my_H5Aclose(hdf5_attribute, "USE_GRACKLE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GRACKLE_CHEMISTRY" , H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 3;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "GRACKLE_CHEMISTRY");
my_H5Aclose(hdf5_attribute, "GRACKLE_CHEMISTRY");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "AGORA_SF" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "AGORA_SF");
my_H5Aclose(hdf5_attribute, "AGORA_SF");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "STARS" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "STARS");
my_H5Aclose(hdf5_attribute, "STARS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "STAR_PARTICLES" , H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 1;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "STAR_PARTICLES");
my_H5Aclose(hdf5_attribute, "STAR_PARTICLES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "WINDS" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "WINDS");
my_H5Aclose(hdf5_attribute, "WINDS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SUPERNOVAE" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SUPERNOVAE");
my_H5Aclose(hdf5_attribute, "SUPERNOVAE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "REGULARIZE_MESH_CM_DRIFT" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REGULARIZE_MESH_CM_DRIFT");
my_H5Aclose(hdf5_attribute, "REGULARIZE_MESH_CM_DRIFT");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED");
my_H5Aclose(hdf5_attribute, "REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "REGULARIZE_MESH_FACE_ANGLE" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REGULARIZE_MESH_FACE_ANGLE");
my_H5Aclose(hdf5_attribute, "REGULARIZE_MESH_FACE_ANGLE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "REFINEMENT_SPLIT_CELLS" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REFINEMENT_SPLIT_CELLS");
my_H5Aclose(hdf5_attribute, "REFINEMENT_SPLIT_CELLS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "REFINEMENT_MERGE_CELLS" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REFINEMENT_MERGE_CELLS");
my_H5Aclose(hdf5_attribute, "REFINEMENT_MERGE_CELLS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "REFINEMENT_VOLUME_LIMIT" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REFINEMENT_VOLUME_LIMIT");
my_H5Aclose(hdf5_attribute, "REFINEMENT_VOLUME_LIMIT");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "NODEREFINE_BACKGROUND_GRID" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "NODEREFINE_BACKGROUND_GRID");
my_H5Aclose(hdf5_attribute, "NODEREFINE_BACKGROUND_GRID");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "COOLING" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "COOLING");
my_H5Aclose(hdf5_attribute, "COOLING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "USE_SFR" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "USE_SFR");
my_H5Aclose(hdf5_attribute, "USE_SFR");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SELFGRAVITY" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SELFGRAVITY");
my_H5Aclose(hdf5_attribute, "SELFGRAVITY");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "HIERARCHICAL_GRAVITY" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "HIERARCHICAL_GRAVITY");
my_H5Aclose(hdf5_attribute, "HIERARCHICAL_GRAVITY");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "CELL_CENTER_GRAVITY" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "CELL_CENTER_GRAVITY");
my_H5Aclose(hdf5_attribute, "CELL_CENTER_GRAVITY");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GRAVITY_NOT_PERIODIC" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GRAVITY_NOT_PERIODIC");
my_H5Aclose(hdf5_attribute, "GRAVITY_NOT_PERIODIC");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "ALLOW_DIRECT_SUMMATION" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "ALLOW_DIRECT_SUMMATION");
my_H5Aclose(hdf5_attribute, "ALLOW_DIRECT_SUMMATION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "DIRECT_SUMMATION_THRESHOLD" , H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 500;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "DIRECT_SUMMATION_THRESHOLD");
my_H5Aclose(hdf5_attribute, "DIRECT_SUMMATION_THRESHOLD");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "NSOFTTYPES" , H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 4;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "NSOFTTYPES");
my_H5Aclose(hdf5_attribute, "NSOFTTYPES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "MULTIPLE_NODE_SOFTENING" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "MULTIPLE_NODE_SOFTENING");
my_H5Aclose(hdf5_attribute, "MULTIPLE_NODE_SOFTENING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "ADAPTIVE_HYDRO_SOFTENING" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "ADAPTIVE_HYDRO_SOFTENING");
my_H5Aclose(hdf5_attribute, "ADAPTIVE_HYDRO_SOFTENING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "TREE_BASED_TIMESTEPS" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "TREE_BASED_TIMESTEPS");
my_H5Aclose(hdf5_attribute, "TREE_BASED_TIMESTEPS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "DOUBLEPRECISION" , H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 1;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "DOUBLEPRECISION");
my_H5Aclose(hdf5_attribute, "DOUBLEPRECISION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "NGB_TREE_DOUBLEPRECISION" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "NGB_TREE_DOUBLEPRECISION");
my_H5Aclose(hdf5_attribute, "NGB_TREE_DOUBLEPRECISION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "PROCESS_TIMES_OF_OUTPUTLIST" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "PROCESS_TIMES_OF_OUTPUTLIST");
my_H5Aclose(hdf5_attribute, "PROCESS_TIMES_OF_OUTPUTLIST");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "HAVE_HDF5" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "HAVE_HDF5");
my_H5Aclose(hdf5_attribute, "HAVE_HDF5");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "DEBUG" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "DEBUG");
my_H5Aclose(hdf5_attribute, "DEBUG");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "OVERRIDE_PEANOGRID_WARNING" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OVERRIDE_PEANOGRID_WARNING");
my_H5Aclose(hdf5_attribute, "OVERRIDE_PEANOGRID_WARNING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

my_H5Tclose(atype);
}
#endif
