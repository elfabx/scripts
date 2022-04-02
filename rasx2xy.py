#!/usr/bin/env python3

import sys
import zipfile
import argparse

parser = argparse.ArgumentParser(
        description="Extract diffraction profile from .rasx file as .xy")
parser.add_argument("-p", action="store_true", help="print result to stdout instead of saving")
parser.add_argument("filename",help="name of rasx file")

args = parser.parse_args()

if not zipfile.is_zipfile(args.filename):
    print("Error:",args.filename,"is not recognised as a rasx file.")
    sys.exit(1)

with zipfile.ZipFile(args.filename) as rasx:
    try:
        with rasx.open('Data0/Profile0.txt') as xy:
            lines = xy.readlines()
    except:
        print("Error: could not find profile data in ", args.filename,".", sep="")
        sys.exit(1)

if (args.p):
    for line in lines:
        f = line.rstrip().decode().split()
        print(f[0],f[1],sep=" ")
else:
    fno = args.filename
    if fno.endswith((".rasx", ".RASX", ".Rasx")):
        fno = fno[:-5]
    fno = fno + ".xy"
    try:
        with open(fno,"w") as of:
            for line in lines:
                f = line.rstrip().decode().split()
                print(f[0],f[1],sep=" ", file=of)
    except:
        print("Error: could not write ", fno, ".", sep="")
        sys.exit(1)

