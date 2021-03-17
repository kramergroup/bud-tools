#!/usr/bin/env python3
import argparse
import os
import sys

def c_value(e):  # Sort via C distance from 0
    return e.split()[2]


def b_value(e):  # Sort via B distance from 0
    return e.split()[1]


def a_value(e):  # Sort via A distance from 0
    return e.split()[0]


def bulk_substitution(argv):
    parser = argparse.ArgumentParser(description='A tool to Replace a central bulk atom'
                                                 'Should only be used with c-vaccuum surfaces as may cause some issues')
    parser.add_argument("-i", "--input", help='Set input POSCAR file, Default POSCAR [in current dir.]',
                        default="POSCAR", type=str)
    parser.add_argument("-c", "--current", help='The current species in the POSCAR that is wanting to be replaced',
                        type=str)
    parser.add_argument("-r", "--replace", help='The species that is desired to be changed', type=str)

    parser.add_argument('-d', dest='direction', type=str, default='c', help="direction of vacuum")

    args = parser.parse_args()

    # We do our own file parsing because i believe that pymatgens Structure method is a bit dodge. Maybe not idk.
    # If you read this and can make this use pymatgen in a simplier way do it.

    l = []
    with open(args.input) as f:
        for line in f:
            l.append(line)

    speccount = list(zip(l[5].split(), l[6].split()))
    if args.current not in dict(speccount):  # Check current is present
        print("Current element {} not found in POSCAR, have you made a typo".format(args.current))
        sys.exit()
    else:
        print("Current element {0} is found in POSCAR - {1} atoms".format(args.current, dict(speccount)['Co']))

    # Check if Selective dynamics line is present if so, shift all species list up 1
    if "S" in l[7]:
        print("Selective dynamics line found, maintaining properties")
        l_co = l[9:]
    else:
        l_co = l[8:]

    l__ = []
    count = 0  # l of l split - this is primitive but w/e
    for i in range(0, len(speccount)):
        print(i)
        l__.append(l_co[count:count + int(speccount[i][1])])
        count += int(speccount[i][1])

    l_chn = l__[l[5].split().index(args.current)]

    if args.direction == 'c':
        l_chn.sort(key=c_value)
    elif args.direction == 'b':
        l_chn.sort(key=b_value)
    elif args.direction == 'a':
        l_chn.sort(key=a_value)
    else:
        print('UNRECOGNISED Vaccuum direction "{}" from -d flag defaulting to c'.format(args.direction))
        l_chn.sort(key=c_value)

    # Note the middle most atom and change it

    s1 = l_chn[len(l_chn )//2]
    s1_a = s1.split("\n")[0] + '#Changed to {} \n'.format(args.replace)
    l.append(s1_a)
    # Strip the original line from the list, take away one from the species count
    vals = [int(t[-1]) for t in speccount]
    l[5] = l[5].split("\n")[0] + " {} \n".format(args.replace)
    l.remove(s1)
    vals[l[5].split().index(args.current)] = vals[l[5].split().index(args.current)] - 1
    vals.append(1)
    l[6] = ' '.join([str(b) for b in vals]) + '\n'
    l[0] = l[0].split("\n")[0] + " # edited with surface_substitution\n"

    os.makedirs(os.curdir + "/bulk_{}_to_{}/".format(args.current, args.replace), exist_ok=True)
    with open(os.curdir + "/bulk_{}_to_{}/POSCAR".format(args.current, args.replace), "w+") as outfile:
        for line in l:
            outfile.write(line)


if __name__ == "__main__":
    bulk_substitution(sys.argv[1:])
