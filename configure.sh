# Configuration script for Camargue
# The script in this file runs a couple subroutine scripts to configure a fresh
# install of Camargue. It should only have to be run once on a fresh download,
# or it can be re-run if you downloaded a barebones install and later added
# some of the external dependencies. There are some basic safety checks, but
# please make sure you've looked at the README first!

# First check to make sure at Concorde has been placed/softlinked in externals
have_cc=$(ls externals | grep 'concorde' | wc -l)

if [ "$have_cc" -eq 0 ]
then
    echo "Error! Please put a softlink to Concorde in externals."
    exit 1
fi

scripts/tsp_header.sh

if [ "$?" -eq 1 ]
then
    echo "Fatal error in configure.sh"
    exit 1
fi

scripts/config_externals.sh

if [ "$?" -eq 0 ]
then
    echo "Configuration appears successful, run make or make test to compile."
    exit 0
else
    echo "Configure may not compile, or will compile without desired choices."
    exit 1
fi
