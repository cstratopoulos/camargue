#!/usr/bin/env python

"""Functions for generating and working with config.hpp files.

The functions in this script are used by the Camargue install script, and its
Makefile, to generate and modify config.hpp files."""

import os
import re

from check_cc import ConfigExcept

# valid keys for a configuration dictionary
valid_keys = {"CMR_USE_OMP",
              "CMR_HAVE_CATCH",
              "CMR_HAVE_OMP",
              "CMR_DO_TESTS",
              "CMR_HAVE_SAFEGMI"}

# preferences for a bare install with no unit tests or externals
bare_prefs = {"CMR_USE_OMP" : 0,
              "CMR_HAVE_CATCH" : 0,
              "CMR_HAVE_OMP" : 0,
              "CMR_DO_TESTS" : 0,
              "CMR_HAVE_SAFEGMI" : 0}

def valid_prefs(pref_dict):
    """Returns true if pref_dict is a dictionary of valid preferences.

    A dictionary of valid preferences consists of all the keys in valid_keys,
    with all values zero or one.

    The dictionary must also not DO/USE any external it does not HAVE"""

    if set(pref_dict.keys()) != valid_keys:
        print "Preferences dictionary must have precisely these valid keys:"
        print valid_keys
        print "Actual keys:"
        print pref_dict.keys()
        return False

    valuset = set(pref_dict.values())
    if not valuset <= {0, 1}:
        print "Dictionary contains nonbinary values: "
        print valuset
        return False

    if pref_dict["CMR_USE_OMP"] and not pref_dict["CMR_HAVE_OMP"]:
        print "Preference dictionary uses OpenMP but doesn't have it"
        return False

    if pref_dict["CMR_DO_TESTS"] and not pref_dict["CMR_HAVE_CATCH"]:
        print "Preference dictionary wants unit tests but doesn't have Catch"
        return False

    return True

def test_dict(pref_dict):
    """Returns a dictionary with the same preferences as pref_dict, but built
    for running unit tests

    Throws an exception if pref_dict is not valid, or if
    pref_dict["CMR_HAVE_OMP"] is false
"""

    if not(valid_prefs(pref_dict)):
        raise \
            ConfigExcept("Tried to generate unit test prefs from invalid prefs")

    if not pref_dict["CMR_HAVE_CATCH"]:
        raise ConfigExcept("Tried to generate unit test dict without Catch")

    result = pref_dict
    result["CMR_DO_TESTS"] = 1
    return result

def txt_grab_prefs(cfg_txt):
    """Returns a preference dictionary from a plain text file with keys in
    valid_keys and binary values. Throws an exception if not valid."""

    if not os.path.isfile(cfg_txt):
        raise ConfigExcept(cfg_txt + " is not a file in " + os.getcwd())
    else:
        print "Grabbing keys from %s" % cfg_txt

    ret_dict = dict()
    with open(cfg_txt) as read_cfg:
        for line in read_cfg:
            split_line = line.split()
            param = str(split_line[0])
            val = int(split_line[-1])
            ret_dict[param] = val

    if not valid_prefs(ret_dict):
        raise ConfigExcept("Invalid preferences read from " + cfg_txt)
    else:
        return ret_dict

def header_grab_prefs(cfg_hpp):
    """Returns a preference dictionary from a C++ header file consisting of
    #include/#define preprocessor directives

    This function is less clear-cut than txt_grab_prefs, since users may disable
    macros by comments, undefs, a mix of both, or, god forbid, they might
    delete a relevant line entirely. Moreover, there may be
    additional macros present which are not governed by this file. Thus, the
    approach is only to set lines which are explicitly treated by a #define
    macro, either just defining something, or defining it to one or zero.

    This function throws if the dictionary so defined is not valid.
    """

    if not os.path.isfile(cfg_hpp):
        raise ConfigExcept(cfg_hpp + " is not a file in " + os.getcwd())
    else:
        print "Grabbing keys from %s" % cfg_hpp

    ret_dict = bare_prefs
    with open(cfg_hpp) as read_cfg:
        for line in read_cfg:
            if re.match("#define", line) or re.match("#undef", line):
                param = None
                for par in valid_keys:
                    if par in line:
                        param = par
                        break

                if param is not None:
                    sp_line = line.split()
                    macro = sp_line[0]
                    if macro == "#undef":
                        ret_dict[param] = 0
                    else:
                        if len(sp_line) == 3:
                            ret_dict[param] = int(sp_line[2])
                        else:
                            ret_dict[param] = 1

    if not valid_prefs(ret_dict):
        raise ConfigExcept("Couldn't get valid prefs from " + cfg_hpp)
    else:
        return ret_dict


def make_hpp(prefs_dict=bare_prefs, incdir="includes/"):
    """Generates a config.hpp file from prefs_dict and config.blank.

    Throws if prefs_dict is invalid, or if config.blank can't be found, or if
    the new file can't be written."""

    cfg_blank = incdir + "config.blank"
    cfg_hpp = incdir + "config.hpp"
    if not os.path.isfile(cfg_blank):
        raise ConfigExcept("Can't find " + cfg_blank)

    if not valid_prefs(prefs_dict):
        raise ConfigExcept("Tried to build config.hpp from invalid prefs")

    if not os.path.isfile(cfg_hpp):
        print "No config.hpp found, creating one from blank"
    else:
        print "%s already exists" % cfg_hpp
        already_same = False
        already_prefs = None
        try:
            already_prefs = header_grab_prefs(cfg_hpp)
        except ConfigExcept:
            pass
        else:
            already_same = already_prefs == prefs_dict

        if already_same:
            print "Existing file already matches desired prefs, \
no need for new one."
            return None
        else:
            cfg_bak = incdir + "config.bak"
            print "Moving to %s to prevent overwrite" % cfg_bak
            try:
                os.rename(cfg_hpp, cfg_bak)
            except OSError:
                print "Couldn't create backup config.hpp"
                print "Please delete or rename manually and try again"
                raise

    polite_warn = ["//>>> config.hpp generated by scripts/gen_config.py\n",
                   "//>>> For best results, please use that script \
instead of manually editing\n"]

    print "Reading from blank and prefs to generate config.hpp...."
    try:
        with open(cfg_blank) as read_blank, open(cfg_hpp, "w") as write_hpp:
            for msg in polite_warn:
                write_hpp.write(msg)

            for line in read_blank:
                if ">>>" in line:
                    continue

                if re.match("#define", line) or re.match("#undef", line):
                    param = None
                    for par in valid_keys:
                        if par in line:
                            param = par
                            break
                    if param is not None:
                        parval = prefs_dict[param]
                        if parval == 0:
                            write_hpp.write("#undef %s\n" % param)
                        else:
                            write_hpp.write("#define %s 1\n" % param)
                    else:
                        write_hpp.write(line)

                else:
                    write_hpp.write(line)
    except IOError as ie:
        raise ConfigExcept(str(ie) + " writing to new config.hpp")
    else:
        print "Generated %s" % cfg_hpp

main_help = \
["Running gen_config as main will automatically generate a config.hpp",
 "The default directory is:",
 os.getcwd() + "/includes/",
 "We look there for prefs in _cfg_prefs.txt, with barebones prefs if not found",
 "A _cfg_prefs.txt should be generated automatically by cmr_install.py",
 "See there for the format if you want to write one manually",
 "--test requires Catch downloaded and a _cfg_prefs.txt with CMR_HAVE_CATCH 1"
]

def extra_parse_help(parser):
    """Prints the help for parser, plus the main_help strings above, and exits
    """

    parser.print_help()
    for msg in main_help:
        print "%s" % msg

    exit(0)

if __name__ == "__main__":
    import argparse
    import sys


    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument("-T", "--test",
                        help="config.hpp for compiling Catch unit tests",
                        action="store_true")
    group.add_argument("-S", "--solve",
                       help="config.hpp for compiling the Camargue solver",
                       action="store_true")
    parser.add_argument("-I", "--incdir", nargs=1,
                        help=\
                        "Relative or absolute path to directory with \
config.blank and _cfg_prefs.txt")

    if len(sys.argv) == 1:
        extra_parse_help(parser)

    args = parser.parse_args()

    incdir = "includes/"

    if args.incdir is not None:
        incdir = args.incdir[0]

    print "Looking for config.blank and _cfg_prefs.txt in %s" % incdir

    cfg_txt = incdir + "_cfg_prefs.txt"
    prefs = bare_prefs

    try:
        prefs = txt_grab_prefs(incdir + "_cfg_prefs.txt")
    except ConfigExcept as ce:
        print "%s trying to read _cfg_prefs.txt" % str(ce)
        print "Just using bare preferences"
        pass


    if args.solve:
        print "Generating config.hpp for Camargue solver..."
    else:
        print "Generating config.hpp for Catch unit tests....."
        try:
            prefs = test_dict(prefs)
        except ConfigExcept as ce:
            print "%s trying to generate config.hpp for unit tests" % str(ce)
            print "catch.hpp needs to be downloaded and _cfg_prefs.txt needs \
to specify its usage"
            exit(1)

    try:
        make_hpp(prefs, incdir)
    except OSError as oe:
        print str(e)
        exit(1)
    except IOError as ie:
        print str(ie)
        exit(1)
    else:
        exit(0)
