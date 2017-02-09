# Script for modifying the safemir code
# This script is meant to make some edits to the safe GMI code of Cook et al.
# Without these, Camargue will almost certainly not compile.
# The modifications are described inline below.

if [ ! -e externals/safemir ]; then
    (>&2 echo "Called edit_safemir but safemir does not appear to exist")
    exit 1
fi

# safe_mir_dbl.cpp, sys_cuts.cpp
#    replacing macro definitions of max and min with safemir_max and
#    safemir_min. Otherwise, the definitions interfere with system headers.
#

SRCDIR=externals/safemir/src

printf 'Changing max/min macros....'
for f in safe_mir_dbl.cpp sys_cuts.hpp; do
    sed 's/ifndef \(m[ai][xn]\)/ifndef safemir_\1/g' "$SRCDIR"/"$f" > tmp_"$f" \
	&& mv tmp_"$f" "$SRCDIR"/"$f"
    if [ "$?" -ne 0 ]; then
	(>&2 echo "Couldn't change macro ifndef in $f")
	exit 1
    fi
    sed 's/ \(m[ai][xn](\)/ safemir_\1/g' "$SRCDIR"/"$f" > tmp_"$f" \
	&& mv tmp_"$f" "$SRCDIR"/"$f"
    if [ "$?" -ne 0 ]; then
	(>&2 echo "Couldn't change definition or usage of macro in $f")
	exit 1
    fi
done

printf 'done\n'

printf 'Changing const qualifiers in util_cuts.hpp....'

sed 's/  \(number_t\* x\)/const \1/g' "$SRCDIR"/util_cuts.hpp > tmp_util_cuts \
    && mv tmp_util_cuts "$SRCDIR"/util_cuts.hpp

if [ "$?" -ne 0 ]; then
    (>&2 echo "Editing failed")
    exit 1
fi

printf 'done\n'

printf 'Removing braced groups from macro in slmem.h....'

sed 's/({/{/g' "$SRCDIR"/slmem.h | sed 's/})/}/g' > tmp_slmem && \
    mv tmp_slmem "$SRCDIR"/slmem.h

if [ "$?" -ne 0 ]; then
    (>&2 echo "Editing failed")
    exit 1
fi

printf 'done\n'

exit 0

