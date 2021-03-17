#!/usr/bin/env python3
import argparse
import sys
import os


def c_value(e):  # Sort via C distance from 0
    return e.split()[2]


def b_value(e):  # Sort via B distance from 0
    return e.split()[1]


def a_value(e):  # Sort via A distance from 0
    return e.split()[0]


def the_key(tol, val):
    return val // tol


def chunks(l, n):
    """Yield n number of sequential chunks from l."""
    d, r = divmod(len(l), n)
    for i in range(n):
        si = (d + 1) * (i if i < r else r) + d * (0 if i < r else i - r)
        yield l[si:si + (d + 1 if i < r else d)]


def dynamic_slab(argv):
    parser = argparse.ArgumentParser(description='A tool to implement partial selective dynamics or to convert a full '
                                                 'poscar into selective dynamics')
    parser.add_argument("-i", "--input", help='Set input POSCAR file, Default POSCAR [in current dir.]',
                        default="POSCAR", type=str)
    parser.add_argument("-l", "--layers", help='The initial number of layers within the slab',
                        type=int)
    parser.add_argument("-w", "--wanted", help='The desired number of layers to be dyanmic on each side', type=int)

    parser.add_argument('-d', dest='direction', type=str, default='c', help="direction of vacuum")

    args = parser.parse_args()

    # Some initial warnings:
    if args.layers % 2 == 0:
        print('total number of layers {} specified by user is non-even, resulting slab will be non-sym'.format(
            args.layers))

    # Once again doing our own file parsing because we dont trust pymatgen and it seems to get slower every day
    l = []
    with open(args.input) as f:
        for line in f:
            l.append(line)

    print("Making a poscar with central {} layers fixed and edge {} relaxed on each side"
          .format(args.layers - (2 * args.wanted), args.wanted))

    if "S" in l[7]:
        print("Selective dynamics line found, the current Ts and Fs will be overwritten - soz")
        l_chn = l[9:]
        l_top = l[1:9]
    else:
        print("Selective dynamics line not found - proceeding")
        l_chn = l[8:]
        l_top = l[1:8]
        l_top.insert(-1, "Selective dynamics \n")

    if args.direction == 'c':
        l_chn.sort(key=c_value)
    elif args.direction == 'b':
        l_chn.sort(key=b_value)
    elif args.direction == 'a':
        l_chn.sort(key=a_value)
    else:
        print('UNRECOGNISED Vaccuum direction "{}" from -d flag defaulting to c'.format(args.direction))

    # ll = [[]] # Since the list is sorted we can do some clever stuff here (not to check all elements together.)
    # v = 0
    # k = 0
    # while v < len(l_chn) -1:
    #    float(l_chn[0].split()[2])
    #    if math.isclose(float(l_chn[v].split()[2]), float(l_chn[v+1].split()[2])):
    #        ll[k].append(float(l_chn[v].split()[2]))
    #    else:
    #        # print(str(eucdis[v][0]) + ' is not near ' + str(eucdis[v + 1][0]))
    #        ll.append([])
    #        k += 1
    #    v += 1
    # if len(ll) > args.layers:
    #    print("Initial dyna-test gave 15 layers (unique {} values, proceeding to do clever stuff)".format())

    # First attempt at doing smart stuff:
    # if len(ll) % args.layers == 0:
    #    print("Since {0} is a multiple of {1} doing cool method 1".format(len(ll), args.layers))

    ll = []
    for i in l_chn:
        ll.append("{}".format(" ".join(i.split()[0:3])))
    ll = list(chunks(ll, args.wanted))

    lfin = []
    for i in ll[0:args.wanted]:
        for k in i:
            lfin.append("{} T T T\n".format(k))
    for i in ll[args.layers:len(ll) - args.wanted + 1]:
        for k in i:
            lfin.append("{} F F F\n".format(k))
    for i in ll[-args.wanted:]:
        for k in i:
            lfin.append("{} T T T\n".format(k))

    lfin = l_top + lfin
    lfin[0] = lfin[0].split("\n")[0] + " # Made dynamic via dyanmic_slab\n"

    os.makedirs(os.curdir + "/dynamic_slab/", exist_ok=True)
    with open(os.curdir + "/dynamic_slab/POSCAR", "w+") as outfile:
        for line in l:
            outfile.write(line)


if __name__ == "__main__":
    dynamic_slab(sys.argv[1:])
