#!/usr/bin/env python

"""Modifies some function/macro definitions in the Cook et al safe MIR code.
This script edits certain source/header files which are included in Camargue,
removing code that is not ISO C++ compliant and/or code that would break the
Camargue compilation.
"""

import re
import os

if not os.path.isdir('externals'):
    print "Error: script called in %s" % os.getcwd()
    print "Please invoke from camargue main directory containing externals/"
    raise OSError("Wrong calling directory")

if not os.path.isdir('externals/safemir/src'):
    print "Can't find safemir/src, please make sure safemir.tar.gz is extracted"
    raise OSError("Missing source directory")

def MIR_full_path(srcfile):
    """Get a path to a file in the safemir/src directory.
    Returns a string, raising an exception if the file can't be found.
    """
    result = "externals/safemir/src/" + srcfile
    if not os.path.isfile(result):
        raise IOError("File to be edited doesn't exist")
    return result

def temp_file_tuple(fname):
    """Get a tuple of filenames for creating temp files in a search/replace
    If called on fname, will return the tuple (fname.new, fname.bak)
    Changes should be written to fname.new, with the old file moved to
    fname.bak, and then fname.new renamed to fname.
    """
    return (fname + ".new", fname + ".bak")

def fix_minmax():
    """Check for definitions of min or max as a macro in global namespace
    Some of the safemir files define 'min' and 'max' as macros, which wreaks
    havoc on anything included afterwards. If these are found, we append an
    undef to the end of the file"""
    print "Checking for issues with min/max macros...."
    for minmax_file in ('safe_mir_dbl.cpp', 'sys_cuts.hpp'):
        def_min = False
        undef_min = False
        def_max = False
        undef_max = False

        open_target = MIR_full_path(minmax_file)

        with open(open_target, "r") as read_minmax:
            for line in read_minmax:
                if '#define min' in line:
                    def_min = True
                elif '#define max' in line:
                    def_max = True
                elif '#undef min' in line:
                    undef_min = True
                elif '#undef max' in line:
                    undef_max = True

        if def_min == undef_min and def_max == undef_max:
            print "Macros in %s look fine" % minmax_file
        else:
            try:
                with open(open_target, "a") as app_minmax:
                    if def_min and not undef_min:
                        app_minmax.write("#undef min\n")
                        print "Undef'd min in %s" % minmax_file
                    if def_max and not undef_max:
                        app_minmax.write("#undef max\n")
                        print "Undef'd max in %s" % minmax_file
            except IOError:
                print "Error editing min/max macros in %s, " \
                    "expect compilation failure" % minmax_file
                raise

def fix_const():
    """Editing const argument qualifiers in util_cuts.hpp
    The header util_cuts.hpp includes several functions which check the
    activity or violation of a cut at some x vector. This x-vector is not
    modified by the routine, so we change the function prototype to declare it
    as const, which agrees with how the lp vector is passed in the calling
    routine in Camargue
    """
    print "Checking for issues with const qualifiers...."
    open_target = MIR_full_path("util_cuts.hpp")
    with open(open_target) as read_util:
        found_bad = False
        for line in read_util:
            if 'number_t* x,' in line and 'const' not in line:
                found_bad = True
                read_util.seek(0)
                break

        if not found_bad:
            print "Const qualifiers already look fine"
        else:
            (new_target, bak_target) = temp_file_tuple(open_target)
            try:
                with open(new_target, "w") as new_util:
                    for line in read_util:
                        if 'number_t* x,' in line and 'const' not in line:
                            newline = line.replace('number_t* x',
                                                   'const number_t* x')
                            new_util.write(newline)
                        else:
                            new_util.write(line)
            except IOError:
                print "Error changing const qualifiers, "\
                    "expect compilation failure"
                raise
            else:
                os.rename(open_target, bak_target)
                os.rename(new_target, open_target)
                print "Changes made, old file at %s" % bak_target

def fix_braced():
    print "Checking for braced groups in malloc macros...."
    open_target = MIR_full_path("slmem.h")
    with open(open_target) as read_mem:
        found_bad = False
        for line in read_mem:
            if '({' in line or '})' in line:
                found_bad = True
                read_mem.seek(0)
                break

        if found_bad is False:
            print "slmem.h macros already look fine"
        else:
            (new_target, bak_target) = temp_file_tuple(open_target)
            try:
                with open(new_target, "w") as new_mem:
                    for line in read_mem:
                        if '({' in line or '})' in line:
                            leftline = line.replace('({', '{')
                            rightline = leftline.replace('})', '}')
                            new_mem.write(rightline)
                        else:
                            new_mem.write(line)
            except IOError:
                print "Error editing macros, project may still compile " \
                    "but with warnings"
                pass
            else:
                os.rename(open_target, bak_target)
                os.rename(new_target, open_target)
                print "Changes made, old file at %s" % bak_target
