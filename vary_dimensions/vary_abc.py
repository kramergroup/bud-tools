#!/usr/bin/env python3
import vary_dimensions
import os, sys
import numpy as np
from pymatgen.core import Lattice, Structure

"""Wrapper to allow a user to make numerous slightly varied unit cells of an input structure both for use both in onetep
 and vasp should allow users to be less worried about whether the ISIF=3 is honest between the two codes"""

__version__ = '0.1'
__author__ = 'BudMacaulay'


def vary_abc(argv):

    #  -------------------------------------------------------------------------------
    # Argument parser
    #  -------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description=vary_abc.__doc__,
                                     epilog="19/03/2021")
    # Positional arguments
    parser.add_argument('filename',
                        default="POSCAR",
                        type=str, nargs='?',
                        help='set input xml file. Default POSCAR;')
    # Optional args
    parser.add_argument('--num',
                        type=int, dest="num", default=5,
                        help='total number of seperate inputs to make; defaults to 5')

    parser.add_argument('--range',
                        type=float, dest='range', default=0.10,
                        help='The fractional range in which the dimension will be varied by; +/-')

    parser.add_argument('-a', '--a',
                        action='store_true',
                        help='Vary a?; defaults to false')
    parser.add_argument('-b', '--b',
                        action='store_true',
                        help='Vary b?; defaults to false')
    parser.add_argument('-c', '--c',
                        action='store_true',
                        help='Vary c?; defaults to false')
    args = parser.parse_args()

    #  Small bit of math here just to ensure that the function is called the right number of times.
    l = np.linspace(-args.range,args.range,args.num)

    count = 0
    for i in l:
        new_struc = vary_lattice_dimensions(args.filename, varyamount=i, a=args.a, b=args.b, c=args.c)
        os.makedirs("/".join(args.filename.split()[:-1]) + str(count))
        new_struc.to("/".join(args.filename.split()[:-1]) + str(count) + "/POSCAR")
        count += 1

if __name__ == "__main__":
    check_conv(sys.argv[1:])
