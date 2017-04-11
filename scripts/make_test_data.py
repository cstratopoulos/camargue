#!/usr/bin/env python

"""Test generation script for Camargue

This script will use Concorde to generate tours/LP solutions that are used in
running Catch test cases.

For each TSP instance in test_data/main_examples.txt, Concorde will be called
to solve over the subtour polytope, saving the LP solution in
test_data/subtour_lp and the tour in test_data/tours.

For the few examples in blossom_examples.txt, Concorde will be called to solve
over the blossom polytope, saving the result in test_data/blossom_lp.

To run this script, make sure:
    Concorde is installed/built and pointed to in camargue/externals/
    A folder containing .tsp files of TSPLIB instances present (or linked to)
    in camargue/problems.
"""

import argparse
import sys
import re
import os
import subprocess

parser = argparse.ArgumentParser()

parser.add_argument("-L", "--large",
                    help="Include large (> 15,000 node) instances",
                    action="store_true")

args = parser.parse_args()

if not os.path.isdir("test_data"):
    print "Can't find test_data dir in %s" % os.getcwd()
    print "Please call this script from the camargue directory"
    exit(1)

if not os.path.isdir("problems"):
    print "Can't find problems folder or symlink in %s" % os.getcwd()
    print "Please link \"problems\" to a folder of TSPLIB examples"
    exit(1)

cc_exec = "externals/concorde/TSP/concorde"

if not os.path.isfile(cc_exec):
    print "Can't find concorde executable %s" % cc_exec
    print "Please build/link concorde into the externals folder"
    exit(1)

def tsp_file(prob):
    """Maps instance names to paths to TSPLIB files in a problems/ folder
    Returns a string path, or None if not found.
    For example, tsp_file("dantizg42") returns "problems/dantzig42.tsp"
    Throws an exception if the problems/ folder doesn't exist
    Returns None if the instance name is not downloaded"""
    if not os.path.isdir("problems"):
        print "Can't find problems folder or symlink in %s" % os.getcwd()
        print "Please link \"problems\" to a folder of TSPLIB examples"
        exit(1)

    result = "problems/" + prob + ".tsp"
    if not(os.path.isfile(result)):
        print "%s doesn't appear to be in your TSPLIB problems folder" % prob
        return None
    else:
        return result

ex_txt = "test_data/main_examples.txt"
tour_dir = "test_data/tours/"
sublp_dir = "test_data/subtour_lp/"

print "Getting subtour LP examples/tours...."
with open(ex_txt) as read_examples:
    seed = "-s99"
    flags = "-Ix"
    for prob in read_examples:
        size = int(re.sub('[^0-9]', '', prob))
        if size > 15000 and not args.large:
            print "Skipping large instance %s" % prob
            continue
        prob = re.sub('[^0-9a-zA-Z]', '', prob)
        lp_dest = sublp_dir + prob + ".sub.x"
        tour_raw = prob + ".sol"
        tour_dest = tour_dir + tour_raw
        tspfile = tsp_file(prob)
        if tspfile is None:
            continue

        if os.path.isfile(lp_dest) and os.path.isfile(tour_dest):
            print "Already have data for %s, skipping" % prob
            continue

        retcode = 1

        with open(os.devnull, 'w') as write_null:
            subprocess.call([cc_exec, seed, flags, "-X", lp_dest, tspfile],
                            stdout = write_null, stderr = write_null)

            if os.path.isfile(tour_raw) and os.path.isfile(lp_dest):
                try:
                    os.rename(tour_raw, tour_dest)
                except OSError:
                    print "Couldn't move tour file %s to test_data/tours" % \
                        tour_raw
                    pass
                else:
                    print "Created subtour lp data for %s" % prob
            else:
                print "Subtour data creation for %s failed" % prob

bloss_txt = "test_data/blossom_examples.txt"

if not os.path.isfile(bloss_txt):
    print "Can't find %s, exiting without generating blossom lp examples" % \
        bloss_txt
    exit(1)

bloss_dir = "test_data/blossom_lp/"

if not os.path.isdir(bloss_dir):
    try:
        os.makedirs(bloss_dir)
    except OSError:
        print "Couldn't make directory for blossom lp solutions"
        exit(1)

print "Getting blossom LP examples/tours...."
with open(bloss_txt) as read_examples:
    seed = "-s99"
    flags = "-ix"
    for prob in read_examples:
        prob = re.sub('[^0-9a-zA-Z]', '', prob)
        lp_dest = bloss_dir + prob + ".2m.x"
        if os.path.isfile(lp_dest):
            print "Already have data for %s, skipping" % prob
            continue
        tspfile = tsp_file(prob)
        if tspfile is None:
            continue
        with open(os.devnull, 'w') as write_null:
            subprocess.call([cc_exec, seed, flags, "-X", lp_dest, tspfile],
                            stdout = write_null, stderr = write_null)

        if os.path.isfile(lp_dest):
            print "Created blossom LP data for %s" % prob
        else:
            print "Error creating blossom LP data for %s" % prob

        try:
            os.remove(prob + ".sol")
        except OSError:
            print "Error removing %s tour cluttering pwd, oops!" % prob
            pass

exit(0)
