#!/usr/bin/env python3
import sys, yaml

if len(sys.argv) != 3:
    print("Usage: python parse_yaml.py <input_yaml> <output_cmake>")
    sys.exit(1)

def yaml_to_cmake(input: str, output: str)-> None:
    with open(input) as stream:
        result = yaml.safe_load(stream)
        with open(output, "w") as f:
            for k,v in result.items():
                if v is True:
                    v = "ON"
                    cmake_option = f"option(AREPO_ENABLE_{k.upper()} {v})"
                elif v is False:
                    v = "OFF"
                    cmake_option = f"option(AREPO_ENABLE_{k.upper()} {v})"
                else:
                    cmake_option = f"set(AREPO_{k.upper()} {v})"

                f.write(cmake_option + "\n")

if __name__ == "__main__":
    yaml_to_cmake(sys.argv[1], sys.argv[2])
    # try:
    #     result = yaml.safe_load(stream)
    #     for k,v in result.items():
    #         print(f"AREPO_ENABLE_({k.upper()}=\"{k} option\" {v})")
    # except yaml.YAMLError as exc:
    #     result = exc
    #     print(result)
