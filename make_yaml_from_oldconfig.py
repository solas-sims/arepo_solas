#!/usr/bin/env python3
import re
import os
import sys
import argparse

def make_yaml(
        input_file: str,
        output_file: str=""
        ) -> None:
    if not os.path.isfile(input_file) and os.access(input_file, os.R_OK):
        raise ValueError(f"File {input_file} is not a readable file.")        
    print(f"Processing file {input_file}")
    if output_file == "":
        output_file=input_file+".yaml"

    # all config patterns either are an on switch or =some value 
    pattern1 = r"^[A-Z][A-Z0-9_]*"
    pattern2 = r"^[A-Z][A-Z0-9_]*+="
    pattern3 = r"#[A-Z][A-Z_]*"

    # Compile the regex for efficiency
    regex1 = re.compile(pattern1)
    regex2 = re.compile(pattern2)
    regex3 = re.compile(pattern3)

    # --- READ AND PROCESS THE FILE ---
    ofile = open(output_file, 'w')
    with open(input_file, 'r') as ifile:
        for line_number, line in enumerate(ifile, start=1):
            # Strip any leading/trailing whitespace from the line
            stripped_line = line.strip()

            # Check if the stripped line matches our pattern
            if  "#!/bin/bash" in stripped_line:
                continue 
            if regex2.search(stripped_line):
                newline = stripped_line.lower().replace("=", ": " , 1)
                ofile.write(newline+'\n')
            elif regex1.search(stripped_line):
                size = len(regex1.search(stripped_line).group())
                newline = stripped_line[:size].lower()
                newline += ": ON " + stripped_line[size:]
                ofile.write(newline+'\n')
            elif regex3.search(stripped_line):
                size = len(regex3.search(stripped_line).group())
                newline = stripped_line[1:size].lower()
                if len(stripped_line) > size:
                    if stripped_line[size] == "=":
                        newline = "#" + newline + stripped_line[size:]
                    else:
                        newline += ": OFF " + stripped_line[size:]
                else:
                    newline += ": OFF "
                ofile.write(newline+'\n')
            else:
                ofile.write(stripped_line+'\n')

# --- RUN THE SCRIPT ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",
                        "--input_file",
                        type=str,
                        help="Input file",
                        required=True,
                        )
    parser.add_argument("-o",
                        "--output_file",
                        type=str,
                        help="Output file [default will be <input_file>.yaml]",
                        default="",
                        )
    args = parser.parse_args()
    make_yaml(args.input_file, args.output_file)

