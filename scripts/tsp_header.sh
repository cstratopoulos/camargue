# Modifying Concorde INCLUDE/tsp.h header file
#
# The keyword 'new' is reserved in C++ as something of a C++ version of malloc.
# There are a few places in the Concorde tsp.h header where new appears as
# the name of a variable in a function prototype. This script will attempt to
# change those in a sensible way, and in such a way that accidentally applying
# the script twice won't start to mangle the header file.
#
# The approach is to replace declarations such as
#
# CCtsp_lpcut_in *new with CCtsp_lpcut_in *new_lpcut_in.
#
# And similarly with lpclique, lpdomino, and skeleton.
# This script should only be run after a softlink to the concorde directory
# is present in externals/

header=externals/concorde/INCLUDE/tsp.h

sed 's/CCtsp\([_a-z]\) \*new)/CCtsp\1 *new\1)/g' "$header" > tmp_header \
    && mv tmp_header "$header"

if [ "$?" -ne 0 ]
then
    echo "Attempt to remove new from concorde/INCLUDES/tsp.h failed"
    exit 1
else
    exit 0
fi


