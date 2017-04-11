#!/usr/bin/env python

"""Usage: ./tour_to_sol.py filename.tour
Converts input to filename.sol in same directory.
"""

import sys
import os.path
import re

arg_count = len(sys.argv)

if arg_count != 2:
    print "Wrong number of arguments specified"
    print __doc__
    exit(1)

(filepath, tour_fname) = os.path.split(sys.argv[1])

split_tname = tour_fname.rsplit(".", 1)

if len(split_tname) != 2 or split_tname[-1] != "tour":
    print "Tour file %s doesn't appear to have right format/suffix" % tour_fname
    print __doc__
    exit(1)

sol_fname = tour_fname.rsplit(".", 1)[0] + ".sol"
write_fname = None

if len(filepath) == 0:
    write_fname = sol_fname
else:
    write_fname = "/".join(["%s" % elem for elem in [filepath, sol_fname]])

if os.path.isfile(write_fname):
    print "Output %s already exists, please move to prevent overwriting" % \
        write_fname
    exit(1)

dimension = None

with open(sys.argv[1]) as read_tour:
    for line in read_tour:
        if 'DIMENSION' in line:
            dimension = int(re.sub('[^0-9]', '', line))
            break

    if dimension is None:
        print "Couldn't find DIMENSION line in tour file"
        exit(1)

    with open(write_fname, "w") as write_sol:
        write_sol.write("%d\n" % dimension)
        num_scanned = 0
        for line in read_tour:
            if re.search('^[0-9]', line) is not None:
                num_scanned += 1
                node = int(re.sub('[^0-9]', '', line)) - 1
                write_sol.write("%d " % node)
                if (num_scanned % 10) == 0:
                    write_sol.write("\n")

        if num_scanned != dimension:
            print "Warning: dimension mismatch in output file!\n" \
                "Tour file had dimension %d, but only " \
                "%d nodes written to output" % (dimension, num_scanned)

print "Converted %s to %s" % (sys.argv[1], write_fname)
