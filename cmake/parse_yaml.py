#!/usr/bin/env python3
import sys, yaml

if len(sys.argv) != 3:
    print("Usage: python parse_yaml.py <input_yaml> <build_dir>")
    sys.exit(1)

def yaml_to_cmake(config_file: str, build_dir: str)-> None:
    output = f"{build_dir}/generated_options.cmake"
    with open(config_file) as stream:
        result = yaml.safe_load(stream)
        with open(output, "w") as f:
            for k,v in result.items():
                if v is True or v is False:
                    if v is True: 
                        v = "ON"
                        add_definitions = f"add_definitions(-D{k.upper()})"
                    else: 
                        v = "OFF"
                        add_definitions = ""
                    cmake_option = f"option(AREPO_ENABLE_{k.upper()} {v})"
                else:
                    cmake_option = f"set(AREPO_{k.upper()} {v})"
                    add_definitions = f"add_definitions(-D{k.upper()}={v})"

                f.write(cmake_option + "\n")
                f.write(add_definitions + "\n")

def print_hdf5_attribute(k, v) -> str:
    outstring : str = ""
    outstring += "hdf5_dataspace = my_H5Screate(H5S_SCALAR);\n";
    if v is True or v is False: 
        if v is True:
            outstring += f"hdf5_attribute = my_H5Acreate(handle, \"{k.upper()}\" , atype, hdf5_dataspace, H5P_DEFAULT);\n"
            outstring += f"my_H5Awrite(hdf5_attribute, atype, \"\", \"$subfields[0]\");\n"
    else:
        outstring += f"hdf5_attribute = my_H5Acreate(handle, \"{k.upper()}\" , H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);\n"
        outstring += f"val = {v};\n"
        outstring += f"my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, \"$subfields[0]\");\n"
    outstring += f"my_H5Aclose(hdf5_attribute, \"{k.upper()}\");\n"
    outstring += f"my_H5Sclose(hdf5_dataspace, H5S_SCALAR);\n\n"
    return outstring


def yaml_to_config_src_files(config_file : str, build_dir: str) -> None:
    arepoconfig_content = "#ifndef AREPO_CONFIG_H\n#define AREPO_CONFIG_H\n"
    compile_time_info_content = """#include <stdio.h>
void output_compile_time_options(void)
{
    printf(
"""
    compile_time_info_hdf5_content = """#include <stdio.h>
#include "arepoconfig.h"

#ifdef HAVE_HDF5
#include <hdf5.h>

hid_t my_H5Acreate(hid_t loc_id, const char *attr_name, hid_t type_id, hid_t space_id, hid_t acpl_id);
hid_t my_H5Screate(H5S_class_t type);
herr_t my_H5Aclose(hid_t attr_id, const char *attr_name);
herr_t my_H5Awrite(hid_t attr_id, hid_t mem_type_id, const void *buf, const char *attr_name);
herr_t my_H5Sclose(hid_t dataspace_id, H5S_class_t type);
herr_t my_H5Tclose(hid_t type_id);


void write_compile_time_options_in_hdf5(hid_t handle) {
hid_t hdf5_dataspace, hdf5_attribute;
double val;
hid_t atype = H5Tcopy(H5T_C_S1);
H5Tset_size(atype, 1);
"""

    with open(config_file) as stream:
        result = yaml.safe_load(stream)
    
    with open(f"{build_dir}/arepoconfig.h", 'w') as arepoconfig_file, \
         open(f"{build_dir}/compile_time_info.c", 'w') as compile_time_info_file, \
         open(f"{build_dir}/compile_time_info_hdf5.c", 'w') as compile_time_info_hdf5_file:

        arepoconfig_file.write(arepoconfig_content)
        compile_time_info_file.write(compile_time_info_content)
        compile_time_info_hdf5_file.write(compile_time_info_hdf5_content)

        print("Processing configuration options...")
        for k,v in result.items():
            compile_time_info_hdf5_file.write(print_hdf5_attribute(k, v))
            if v is True or v is False:
                if v is True:
                    arepoconfig_file.write(f"#define {k.upper()}\n")
                    compile_time_info_file.write(f"\"   {k.upper()} \\n\"\n")
            else:
                arepoconfig_file.write(f"#define {k.upper()} {v}\n")
                compile_time_info_file.write(f"\"   {k.upper()}={v} \\n\"\n")
        arepoconfig_file.write("#endif\n")
        compile_time_info_file.write(");\n")
        compile_time_info_file.write("}\n")

        compile_time_info_hdf5_file.write("my_H5Tclose(atype);\n")
        compile_time_info_hdf5_file.write("}\n")
        compile_time_info_hdf5_file.write("#endif\n")

if __name__ == "__main__":
    yaml_to_cmake(sys.argv[1], sys.argv[2])
    yaml_to_config_src_files(sys.argv[1], sys.argv[2])
