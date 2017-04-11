#!/usr/bin/env python

"""Some simple checks for compiler specification from a Makefile
The functions in this file can be used to check if a C++ compiler has been
validly specified from a Makefile in the current directory. There are also
functions to test the compiler, to see whether it is a valid C++ compiler at
the most basic level (with a hello world program), and moreover whether it
compiles with some basics of the OpenMP standard.
"""

import re
import os
import subprocess

class ConfigExcept(Exception):
    """Custom exception type for indicating config/setup error"""
    pass


def check_def_line(def_name):
    """Checks for a definition in a Makefile in this directory

    Throws an exception if cwd contains no Makefile. If it does, searches for
    a line of the form
    def_name := [definition]
    defining def_name. Returns [definition] if found, throwing ConfigExcept if
    not found, or if the definition looks clearly invalid."""

    if not os.path.isfile('Makefile'):
        print "Error: no Makefile in path, script called from:"
        print os.getcwd()
        raise ConfigExcept("Makefile doesn't exist in curent directory")

    re_string = def_name + " +:="
    with open("Makefile") as read_make:
        for line in read_make:
            if re.match(re_string, line):
                decomp_line = line.split()
                if len(decomp_line) != 3:
                    print "%s should be defined in the following line:" \
                        % def_name
                    print line
                    raise ConfigExcept("Makefile doesn't define " + def_name)
                else:
                    return decomp_line[-1]

    raise ConfigExcept("Reached EOF searching Makefile for " + def_name)

def get_CC_name():
    """Grabs compiler name from the Makefile"""
    return check_def_line("CC")

def CPX_defs():
    """Checks if CPLEX definitions are specified in Makefile

    Camargue requires setting the path to cplex.h, and the file libcplex.a
    itself. This function will check if both those are specified.
    """
    print "Checking for CPLEX definitions in Makefile...."
    try:
        dir_line = check_def_line("CPXDIR")
    except ConfigExcept as e:
        print "%s looking for cplex.h directory" % str(e)
        raise ConfigExcept("Path to cplex.h dir not specified in Makefile")
    else:
        if not os.path.isfile(dir_line + "cplex.h"):
            raise ConfigExcept("Couldn't find cplex.h in " + dir_line)
        else:
            print "Found directory with cplex.h"

    try:
        lib_line = check_def_line("CPX_LIB")
    except ConfigExcept as e:
        print "%s looking for path to libcplex.a" % str(e)
        raise ConfigExcept("Path to libcplex.a not specified in Makefile")
    else:
        if os.path.split(lib_line)[-1] != "libcplex.a":
            raise ConfigExcept(lib_line + " is not a path to libcplex.a")
        else:
            print "Found libcplex.a"


def hello_world_string():
    """Generates a hello world program"""
    codelines = ["#include <iostream>",
                 "int main()", "{",
                 "std::cout << \"henlo stinky world, \"",
                 "<< \"C++ compiler is functional\" << std::endl;",
                 "return 0;",
                 "}"]
    return "\n".join(codelines)

def test_omp_string():
    """Generates a trivial OpenMP program"""
    codelines = ["#include <iostream>",
                 "#include <omp.h>",
                 "int main()", "{",
                 "std::cout << \"Compiler supports OpenMP, \"",
                 "<< \"max thread count \" << omp_get_max_threads()",
                 "<< std::endl;",
                 "return 0;",
                 "}"]
    return "\n".join(codelines)


def compiler_check(cc_name, code_string, link_option=None):
    """Invoke a chosen compiler on a string corresponding to a program
    Returns true if compilation succeeds, and if program is executed with
    return value of zero. A linker option may optionally be specified as
    an additional argument to the compiler, to be used for example in
    testing OpenMP support.
    Returns false if an error occurred during compilation, or if program
    executes with nonzero retcode.
    Implementation adapted from:
    http://stackoverflow.com/q/4293297/6516346
    http://stackoverflow.com/a/11269627/6516346"""

    cc_check_src = "check_test_prog.cpp"
    cc_check_exec = "check_test.out"

    with open(cc_check_src, "w") as write_src:
        try:
            write_src.write(code_string)
        except IOError:
            print "Couldn't create test program source file"
            raise

    cc_cmd = [cc_name, cc_check_src, "-o", cc_check_exec]
    if link_option is not None:
        cc_cmd.append(link_option)

    result = False
    retcode = 1

    with open(os.devnull, 'w') as write_null:
        retcode = subprocess.call(cc_cmd,
                                  stdout = write_null, stderr = write_null)

    if retcode == 0:
        print "Compilation successful, executing"
        retcode = subprocess.call(["./" + cc_check_exec])
        if retcode == 0:
            result = True
        else:
            print "Program compiled but terminated with an error"
    else:
        print "Compilation check failed"

    try:
        if os.path.isfile(cc_check_src):
            os.remove(cc_check_src)
        if os.path.isfile(cc_check_exec):
            os.remove(cc_check_exec)
    except OSError:
        print "Error deleting temporary program files"
        pass

    return result

def smoke_test(compiler="c++"):
    """A smoke test for the C++ compiler: compile/run a Hello World program"""
    print "Checking if %s is a valid C++ compiler...." % compiler
    result = compiler_check(compiler, hello_world_string())
    if not result:
        print "%s does not seem to be a valid C++ compiler" % compiler
    return result

def omp_test(compiler="c++"):
    """OpenMP test for the C++ compiler: compile/run an OpenMP program"""
    print "Checking if compiler %s supports OpenMP..." % compiler
    result = compiler_check(compiler, test_omp_string(), "-fopenmp")
    if not result:
        print "%s does not support OpenMP" % compiler
    return result

if __name__ == "__main__":
    print "Running check_cc as main, testing Makefile or default c++ compiler"
    compiler = None
    if os.path.isfile("Makefile"):
        try:
            compiler = get_CC_name()
        except ConfigExcept:
            pass
        else:
            print "Grabbed compiler %s from Makefile in pwd" % compiler

    if compiler is None:
        print "No Makefile in pwd, testing default compiler \"c++\""
        compiler = "c++"

    if smoke_test(compiler) is True:
        omp_test(compiler)
