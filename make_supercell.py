#!/usr/bin/env python3

import os
import sys
import argparse
import logging
import json
from pymatgen import Structure


def make_supercell(argv):
    parser = argparse.ArgumentParser(description='A tool to create supercells of poscars')
    parser.add_argument("--scd", help='supercell dimension as 3 space-seperated numbers ie. "2 2 2".',
                        nargs="+", type=int, default=[2, 2, 2])
    parser.add_argument("-i", "--input", help='Set input POSCAR file, Default POSCAR [in current dir.]',
                        default="POSCAR", type=str)
    args = parser.parse_args()

    os.makedirs(os.curdir + '/sup' + str(args.scd[0]) + str(args.scd[1]) + str(args.scd[2]), exist_ok=True)
    if len(args.scd) != 3:
        print("Supercell dimensions =! 3. Has form [A, B, C]")
        sys.exit()

    print(args.scd)

    struc = Structure.from_file(args.input)
    struc.make_supercell(args.scd)
    struc.to(filename=(os.curdir + '/sup' + str(args.scd[0]) + str(args.scd[1]) + str(args.scd[2]) + '/POSCAR'))


if __name__ == "__main__":
    make_supercell(sys.argv[1:])
