#!/usr/bin/env python

import os
from pathlib import Path
import sys

dir_path = Path(r"{}".format(os.path.dirname(os.path.realpath(__file__)))).parent

try:
    data_dir = sys.argv[1]
except IndexError:
    print("[DataExplorer] Error! Please specify the organism name as a command line argument as `data/organism_name`")

def get_data_files(data_dir):
    return ["{}/{}".format(data_dir, name) for name in os.listdir(data_dir)]

def replace_whitespace(data_dir):
    files = get_data_files("{}/{}".format(dir_path, data_dir))
    delim = " "

    for fname in files:
        if "protein.sequences" in fname:
            with open(fname, "r+") as file:
                if line.strip() == "": continue

                lines = {}
                protein_id = ""
                
                for line in file:
                    if line.startswith(">"):
                        protein_id = line[line.index(">") + 1:].strip()
                        lines[protein_id] = ""
                    else:
                        lines[protein_id] += line.strip()
            with open(fname, "w+") as file:
                for key in lines:
                    if line.strip() == "": continue

                    try:
                        file.write("{} {}\n".format(key, lines[key]))
                    except e:
                        print("[DataExplorer] Error! {}. Continuing...".format(e))

            continue

        lines = []
        with open(fname, "r+") as file:
            for line in file:
                if line.strip() == "": continue

                try:
                    lines.append(delim.join(line.split()))
                except e:
                    print("[DataExplorer] Error! {}. Continuing...".format(e))

        with open(fname, "w+") as file:
            for line in lines:
                if line.strip() == "": continue

                try:
                    file.write("{}\n".format(line))
                except Exception as e:
                    print("[DataExplorer] Error! {}. Continuing...".format(e))

replace_whitespace(data_dir)