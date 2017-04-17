#!/usr/bin/env python

"""Script with main function for automatically running Camargue in Batches.

This script can be used to invoke Camargue, writing results to file, depending
on success or failure.

The -p, --prob option is for a straightforward, single invocation. For example,
    ./cmr_boo_yay.py --prob d493
is like running
    ./camargue problems/d493.tsp,
but with output logged to file, and suffix changed upon completion.
    d493.yay is successful solution,
    d493.boo is interruption (by timeout) or failure.

The -f, --file option is used for batch execution from a list of problems. Say
probs.txt is a file with instance names ("lin318", "d493", "pr1002", etc.) on
each line. So, since "d493" is the 2nd line of probs.txt,
    ./cmr_boo_yay.py -f 2 probs.txt
would be equivalent to running with "--prob d493" as above.
"""

import argparse
import sys
import os
import subprocess
import time

def tsp_path(inst_name):
    """Returns the string problems/inst_name.tsp, throwing if non-existant"""
    result = "problems/" + inst_name + ".tsp"
    if not os.path.isfile(result):
        raise IOError("Can't find problem at " + result)
    return result

def log_tuple(inst_name):
    """Returns a tuple of files for logging a run on inst_name"""
    return (inst_name + ".run", inst_name + ".yay", inst_name + ".boo")

def run_camargue(inst_name, seed=None):
    """Runs Camargue on inst_name with a time limit if specified"""

    fpath = tsp_path(inst_name)
    (run_f, yay_f, boo_f) = log_tuple(inst_name)

    cmr_cmd = ["./camargue", "-T", fpath]
    if seed is not None:
        cmr_cmd.append("-s" + str(seed))

    print "At time %s, running: " % time.ctime()
    print cmr_cmd

    retcode = 1

    with open(run_f, 'a') as write_run:
        retcode = subprocess.call(cmr_cmd, stdout=write_run, stderr=write_run)

    if retcode == 0:
        print "Completed at %s, written to %s" % (time.ctime(), yay_f)
        os.rename(run_f, yay_f)
    else:
        print "Completed at %s, written to %s" % (time.ctime(), boo_f)
        os.rename(run_f, boo_f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument("-p", "--prob", nargs=1,
                       help="Specify problem name with no tsp suffix")

    group.add_argument("-f", "--file", nargs=2,
                       metavar=('LINENUM', 'FILENAME'),
                       help="Get problem name from line LINENUM of FILENAME")

    parser.add_argument("-s", "--seed", nargs=1,
                        help="Fix the random seed for execution")

    if len(sys.argv) == 1:
        parser.print_help()
        print
        print __doc__
        exit(0)

    args = parser.parse_args()

    pname = None

    if args.prob:
        pname = args.prob[0]
    elif args.file:
        LINE_NUM = int(args.file[0])
        FNAME = args.file[1]
        with open(FNAME) as read_probs:
            lines = read_probs.readlines()
            pname = lines[LINE_NUM - 1].strip()

    if args.seed:
        run_camargue(pname, int(args.seed[0]))
    else:
        run_camargue(pname)
