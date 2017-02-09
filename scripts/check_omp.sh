# OpenMP Configuration routine for Camargue
# The script in this file will attempt to check if the user-specified compiler
# supports OpenMP (OMP) parallelism. It makes no changes to the config file
# itself; that is done by config_externals.sh which invokes this as a
# subroutine.
#
# The approach is to check the Makefile for the C++ compiler, and attempt to
# use it to compile a dummy temporary program with the openmp flag added. An
# exit 0 compilation will be used to indicate that the compiler supports OMP.
# exit 1 means that the compiler appears nonfunctional and may have been
# grabbed incorrectly
# exit 2 means the test script ran but there appears to be no OMP support

if [ -f Makefile ]
then
    echo "Found Makefile, attempting to configure for OMP...."
else
    (>&2 echo "Error: no Makefile in main directory")
    exit 1
fi

CC=$(grep '^CC *:=' Makefile | sed 's/.*:= //')

cwords=$(echo "$CC" | wc -w)

if [ "$cwords" -eq 0 ]; then
    (>&2 echo "Script grabbed no characters for compiler variable")
    exit 1
fi

echo "Compiler is $CC"

"$CC" --version 2>&1 > /dev/null

if [ "$?" -ne 0 ]
then
    echo Basic compiler test failed
    exit 1
fi

echo '#include <iostream>
      #include <omp.h>
 int main(void) {
     std::cout << "Compiler supports OMP, max thread count "
     	       << omp_get_max_threads();
}' > tmp_prog.cpp

"$CC" -fopenmp tmp_prog.cpp -o tmp_prog.o 2>/dev/null && ./tmp_prog.o && echo

worked="$?"

rm tmp_prog*

if [ "$worked" -eq 0 ]
then
    exit 0
else
    exit 2
fi
