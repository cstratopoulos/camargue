#!/usr/bin/env python

"""Modifying Concorde INCLUDE/tsp.h header file.

There are a few places in concorde/INCLUDE/tsp.h where '*new' appears in
function prototypes. In C++, 'new' is a special keyword for memory allocation,
so we edit tsp.h to prevent compilation failures.

Declarations such as

CCtsp_lpcut_in *new

will be replaced with

CCtsp_lpcut_in *new_lpcut_in.

Thus, the changes will keep the header file in a sensible state, and a backup
of the old version will be kept. Moreover, no changes will be made if no such
declarations are found."""

import re
import os

tsp_h = 'externals/concorde/INCLUDE/tsp.h'

def check_extracted():
    """Basic checks to make sure Concorde is extracted and the header exists"""
    if not os.path.isdir('externals'):
        print "Error: no externals/ folder in path, script called from:"
        print os.getcwd()
        print "Please call script from camargue/ folder"
        raise OSError("Called from wrong directory")

    if not os.path.isfile(tsp_h):
        print "Error: can't find Concorde tsp.h"
        print "Create a symlink to the Concorde folder before running"
        raise OSError("Concorde not extracted/symlinked")


def do_edits():
    """Perform the actual edits to the header file"""
    check_extracted()

    print "Checking for 'new' in concorde's tsp.h"
    try:
        with open(tsp_h) as read_tsp:
            found_bad = False
            for line in read_tsp:
                if ' *new,' in line or ' *new)' in line:
                    found_bad = True
                    read_tsp.seek(0)
                    break

            if not found_bad:
                print "tsp.h already looks fine"
                return None

            tsp_h_new = tsp_h + '.new'
            tsp_h_bak = tsp_h + '.bak'
            with open(tsp_h_new, 'w') as write_tsp:
                for line in read_tsp:
                    if ' *new,' in line or ' *new)' in line:
                        newline = re.sub(r"CCtsp([_a-z]+) \*new([),])",
                                             r"CCtsp\1 *new\1\2", line)
                        write_tsp.write(newline)
                    else:
                        write_tsp.write(line)
    except IOError as e:
        print "%s editing tsp.h" % str(e)
        raise IOError("Error editing tsp.h, project may not compile")
    else:
        print "Performed edits to tsp.h"
        try:
            os.rename(tsp_h, tsp_h_bak)
            os.rename(tsp_h_new, tsp_h)
        except OSError as oe:
            print "Edited successfully, but %s renaming edited file" % str(oe)
            print "You can try to manually rename tsp.h.new to tsp.h"
            raise OSError("Error renaming files, project may not compile")
        else:
            print "Edited and renamed, old header is at %s" % tsp_h_bak
