#!/usr/bin/env python

""" The main Camargue install script.
Accepts command line options to configure the desired Camargue installation"""

import argparse
import sys
import re
import os

sys.path.append('scripts/')
import manage_externals
import check_cc
import tsp_header

from check_cc import ConfigExcept

parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)

group.add_argument("-B", "--bare",
                   help="Bare install: no externals besides Concorde/CPLEX",
                   action="store_true")

group.add_argument("-F", "--full",
                   help="Full install with all external dependencies",
                   action="store_true")

group.add_argument("-W", "--which", nargs='+',
                   help="Specify precisely Which externals to use")

def print_Which():
    print
    print "The option -W, --which takes 1-3 arguments:"
    print "\tcatch\t\tCatch unit testing"
    print "\tsafemir\t\tSafe GMI implementation"
    print "\tomp\t\tOpenMP parallelism for some separation routines"

if len(sys.argv) == 1:
    parser.print_help()
    print_Which()
    exit(1)

args = parser.parse_args()

want_catch = False
want_gmi = False
want_omp = False

if args.full:
    print "Full install selected"
    want_catch = True
    want_gmi = True
    want_omp = True
elif args.which:
    for s in args.which:
        if s == "catch":
            print "Catch unit tests selected"
            want_catch = True
        elif s == "safemir":
            print "SafeGMI implementation selected"
            want_gmi = True
        elif s == "omp":
            print "OpenMP configuration selected"
            want_omp = True;
        else:
            parser.print_help()
            print_Which()
            print "%s does not correspond to an external dependency" % s
            exit(1)
else:
    print "Bare install selected."

print "Performing basic checks/edits...."

try:
    check_cc.smoke_test()
    check_cc.CPX_defs()
    tsp_header.do_edits()
except ConfigExcept as ce:
    print "Error %s in configuration checks" % str(ce)
    print "Please verify Makefile definitions/symlinks to Concorde"
    exit(1)
except IOError as ie:
    print "IO error %s in configuration checks" % str(ie)
    print "There was a problem reading/writing some files"
    exit(1)
except OSError as oe:
    print "OS error %s in configuration checks" % str(oe)
    print "There was a problem moving/renaming some files"
    exit(1)

got_catch = False
got_gmi = False
got_omp = False

if want_catch:
    try:
        manage_externals.get_Catch()
    except Exception as e:
        print "%s trying to download Catch, building without it" % str(e)
        print "Please try again later or maybe do a manual download"
        pass
    else:
        got_catch = True

if want_gmi:
    try:
        manage_externals.get_GMI()
    except Exception as e:
        print "%s downloading/extracting safe Gomory code, building \
without it" % str(e)
        print "Please try again later or maybe do a manual download"
        pass
    else:
        got_gmi = True

    if got_gmi:
        try:
            manage_externals.clean_GMI()
        except Exception as e:
            print "%s sanitizing GMI code, fatal error" % str(e)
            print "These errors can severely break project compilation, so \
disabling GMI in configuration"
            print "If you still want an install with safeGMI, consider \
editing the source manually based on scripts/edit_safemir.py"
            got_gmi = False
            pass

if want_omp:
    try:
        got_omp = check_cc.omp_test()
    except Exception as e:
        print "%s checking for OpenMP support, building without it" % str(e)
        pass
    else:
        print "Checking Makefile for -fopenmp flag"
        if got_omp:
            needs_change = False
            try:
                with open("Makefile") as read_make:
                    for line in read_make:
                        if re.match("FOMP +:=", line):
                            spline = line.split()
                            if len(spline) != 3 or spline[-1] != "-fopenmp":
                                needs_change = True
                                break
            except IOError:
                print "Couldn't check Makefile, building without it"
                print "You can manually add FOMP := -fopenmp if you like"
                pass
            else:
                if not needs_change:
                    print "Makefile already specifies flag"
                else:
                    print "Editing Makefile to add flag"
                    try:
                        with open("Makefile") as read_make, \
                             open("Makefile.new", "w") as new_make:
                            for line in read_make:
                                if re.match("FOMP +:=", line):
                                    new_make.write("FOMP        := -fopenmp\n")
                                else:
                                    new_make.write(line)
                    except IOError:
                        print "Error adding OMP flag to Makefile, \
building without it"
                        got_omp = False
                        pass
                    else:
                        try:
                            os.rename("Makefile.new", "Makefile")
                        except OSError:
                            print "Couldn't write -fopenmp to Makefile, \
        building without it"
                            print "To build with OpenMP, just edit the FOMP line in \
        Makefile to say \"FOMP := -fopenmp\""
                            got_omp = False
                            pass
                        else:
                            print "Makefile edited."


print
print "...Install script completed with the following choices:"
macs = {"CMR_HAVE_CATCH" : int(got_catch),
        "CMR_HAVE_OMP" : int(got_omp),
        "CMR_HAVE_SAFEGMI" : int(got_gmi),
        "CMR_DO_TESTS": 0,
        "CMR_USE_OMP" : int(got_omp and want_omp)}
for k, v in macs.items():
    print "%s %d" % (k, v)

print "Writing preferences to camargue/includes/_cfg_prefs.txt"
try:
    with open("includes/_cfg_prefs.txt", "w") as write_prefs:
        for k, v in macs.items():
            write_prefs.write("%s %d\n" % (k, v))
except IOError:
    print "Error writing preferences to file."
    print "For future use, you can create _cfg_prefs.txt with the data printed \
above"
    print "Run scripts/gen_config.py with no arguments for info"
    pass

import gen_config

try:
    gen_config.make_hpp(macs, "includes/")
except Exception as e:
    print "%s generating config.hpp" % str(e)
    exit(1)
else:
    print "Looks like everything went ok, compile with make"
