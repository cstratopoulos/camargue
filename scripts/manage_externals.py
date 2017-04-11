#!/usr/bin/env python

import os.path
import urllib
import tarfile

catch_fname = "externals/catch.hpp"
gmi_archive = "externals/safemir.tar.gz"
gmi_dir = "externals/safemir"

catch_url = "https://raw.githubusercontent.com/philsquared/Catch/master/single_include/catch.hpp"
gmi_url = "https://www.informs.org/content/download/257165/2428239/file/safemir090309.tar.gz"

def get_Catch():
    """Downloading the Catch header from GitHub"""
    print "Installing with Catch......."

    if os.path.isfile(catch_fname):
        print "...looks like Catch is already downloaded"
    else:
        try:
            urllib.urlretrieve(catch_url, catch_fname)
        except Exception:
            print "Problem downloading Catch header"
        else:
            print "...downloaded Catch header"

def get_GMI():
    """ The actual downloading/extraction of the safe GMI code"""
    print "Installing with Safe GMI code......"
    if os.path.isdir(gmi_dir):
        print "...looks like Safe GMI implementation is already installed."
    else:
        try:
            urllib.urlretrieve(gmi_url, gmi_archive)
        except Exception:
            print "Problem downloading safemir archive"
            dl_probcount +=1
        else:
            try:
                gmitar = tarfile.open(gmi_archive, "r:gz")
                gmitar.extractall("externals")
                gmitar.close()
                os.remove(gmi_archive)
            except tarfile.TarError:
                print "Error extracting safemir tgz"
            except OSError:
                print "Error removing compressed archive"
                pass
            else:
                print "...downloaded and extracted safemir tgz"

def clean_GMI():
    """Invoking the edit functions to sanitize the GMI code"""
    print "Checking for issues in safe GMI code...."
    import edit_safemir
    try:
        edit_safemir.fix_minmax()
        edit_safemir.fix_const()
        edit_safemir.fix_braced()
        edit_safemir.fix_undeclared()
    except OSError:
        print "Error with some of the directory creation/extraction"
        raise
    except IOError:
        print "Fatal error editing some of the safe GMI code"
        raise
