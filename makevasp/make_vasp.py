#!/usr/bin/env python3

import os
import argparse
import sys
from vaspmake_utils import get_potcar, get_simple_incar, incar_add_perfopt, gen_gamma_kpoint

__author__ = 'Bud Macaulay, Andrea Silva (LamaKing) for the potcar util'
__version__ = '0.2'


""" A wrapper for the vaspmake_utils method of generating input files. Care the INCAR is very sparce,
    it's quite heavy but makes an incredibly simple INCAR, KPOINTS and POTCAR
"""


def make_vasp(argv):
    # optional arguments is all cause using argparse for anything else may cause issues

    parser = argparse.ArgumentParser(description="A wrapper for the vaspmake_utils, uses some quite simple"
                                                 "the file flags disable making of "
                                                 "said file")
    # Positional arguments
    parser.add_argument("poscar", default="POSCAR", type=str, nargs='?',
                        help='Set the input file default = POSCAR')
    # Optional arguments
    parser.add_argument('-p', '--ppdir', dest="ppdir",
                        help='set custom ppdirectory else it will use a random one')
    parser.add_argument('--potcar', action='store_false',
                        help='Whether or not to make a very simple POTCAR file')
    parser.add_argument('--incar',
                        action='store_false', dest='incar',
                        help='Whether or not to make a very simple INCAR file')
    parser.add_argument('--kpoints',
                        action='store_false', dest='kpoints',
                        help='Whether or not to make a very simple KPOINT file')

    parser.add_argument('--cpn',
                        default=40, dest='cpn',
                        help='The cores per node on the desired submit cluster (default=40 for IRIDIS5)')

    args = parser.parse_args(argv)

    if args.potcar:
        if not args.ppdir:
            print("ppdirectory not supplied, using default - this may crash if you're not Bud")
            with open(os.curdir + "/POTCAR", "w+") as outfile:
                outfile.write(get_potcar(input_poscar=args.poscar))
        else:
            with open(os.curdir + "/POTCAR", "w+") as outfile:
                outfile.write(get_potcar(input_poscar=args.poscar))
    else:
        pass

    if args.kpoints:
        with open(os.curdir + "/KPOINTS", "w+") as outfile:
            outfile.write(gen_gamma_kpoint())

    if args.incar:
        s = get_simple_incar(args.poscar, magmom=True, u=True)
        s = s + incar_add_perfopt(cpn=40)
        with open(os.curdir + "/INCAR", "w+") as outfile:
            outfile.write(s)

if __name__ == "__main__":
    make_vasp(sys.argv[1:])