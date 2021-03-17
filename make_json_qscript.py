#!/usr/bin/env python3
import argparse
import sys
import os


def make_json_qscript(argv):
    parser = argparse.ArgumentParser(description='A tool to generate a relevant qscript file and a json file in a file '
                                                 'containing a POSCAR, it will try to add a node per atom up '
                                                 'to a max of 120')
    # Positional Arguments
    parser.add_argument("qscript", type=int, default=2, nargs='?',
                        help='The qscript to use [of which there are currently three options '
                             '0: qscript off, 1: MICHAEL, 2: IRIDIS5')

    parser.add_argument('-j', '--schedulertype', type=str, default='sbatch',
                        help='The scheduler type on the HPC in which there are two options "sbatch" and "qsub"')

    parser.add_argument('-t', '--time', type=str, default=20,
                        help='A flag to set the amount of time desired [defaults to 20hrs.]')

    num= 40
    args = parser.parse_args(argv)
    try:
        line = open("{}/POSCAR".format(os.curdir), "r").readlines()[6]
        num = sum([int(s) for s in line.split() if s.isdigit()])
        if num > 120:
            num = 120
    except IOError:
        print("No POSCAR found in directory defaulting cores")


    if args.qscript == 1:
        num = round(num/24) * 24
        with open(os.curdir + "/qscript", "w+") as outfile:
            outfile.write("#!/bin/bash -l\n"
                          "#$ -S /bin/bash\n"
                          "\n"
                          "#$ -l h_rt={0}:00:00 #Made with make_json_qscript \n"
                          "\n"
                          "#Request 1 gigabyte of RAM per process.\n"
                          "\n"
                          "#Request 10G gigabyte of TMPDIR space per node\n"
                          "#$ -l tmpfs=10G"
                          "\n"
                          "#$ -N {1}"
                          "\n"
                          "#Select the MPI parallel environment and 24 processes.\n"
                          "#$ -pe mpi {2}\n"
                          "\n"
                          "#$ -cwd\n"
                          "\n"
                          "#add to Faraday_MSM\n"
                          "#$ -P Gold\n"
                          "#$ -A Faraday_MSM_kra\n"
                          "\n"
                          "gerun $HOME/vasp/vasp_std 2>&1 > print-out\n"
                          "\n"
                          "rm WAVECAR CHGCAR CHG\n".format(args.time, os.getcwd().split('/')[-1], num))

            with open(os.curdir + "/local.json", "w+") as outfile:
                outfile.write('{\n'
                              '\t"sub_cmd": "' + args.schedulertype + '",\n' \
                                                                      '\t"hpc_fld":" ~/' + os.getcwd().split('/')[
                                  -2] + ',\n' \
                                        '\t"hostname": "michael",\n' \
                                        '\t"env_setup": "~/env.sh"\n' \
                                        '}')

    elif args.qscript == 2:
        num = round(num/40) * 40
        with open(os.curdir + "/qscript", "w+") as outfile:
            outfile.write("#!/bin/bash\n"
                          "# Made with make_json_qscript\n"
                          "\n"
                          "#  -- SLURM INFO --\n"
                          "\n"
                          "#SBATCH --ntasks={0}\n"
                          "#SBATCH --job-name={1}\n"
                          "\n"
                          "#SBATCH --time={2}:00:00\n"
                          "\n"
                          "# -- DIREC + MODULES --\n"
                          "\n"
                          "cd $SLURM_SUBMIT_DIR\n"
                          "module purge\n"
                          "module load intel-mpi/2017.4.239 intel-mkl/2017.4.239\n"
                          "\n"
                          "# -- RUN VASP --\n"
                          "\n"
                          "mpirun $HOME/vasp/vasp_std 2>&1 > print-out\n"
                          "\n"
                          "# -- CLEANUP --\n"
                          "rm WAVECAR CHG CHGCAR\n".format(num, os.getcwd().split('/')[-1], args.time))

        with open(os.curdir + "/local.json", "w+") as outfile:
            outfile.write('{\n' \
                          '\t"sub_cmd": "' + args.schedulertype + '",\n' \
                                                                  '\t"hpc_fld": "~/' + os.getcwd().split('/')[
                              -2] + '",\n' \
                                    '\t"hostname": "iridis5",\n' \
                                    '\t"env_setup": "~/env.sh"\n' \
                                    '}')
    else:
        with open(os.curdir + "/locar.json", "w+") as outfile:
            outfile.write('{\n' \
                          '\t"sub_cmd": "' + args.schedulertype + '",\n' \
                                                                  '\t"hpc_fld": "~/' + os.getcwd().split('/')[
                              -2] + '",\n' \
                                    '\t"hostname": "iridis5",\n' \
                                    '\t"env_setup": "~/env.sh"\n' \
                                    '}')


if __name__ == "__main__":
    make_json_qscript(sys.argv[1:])
