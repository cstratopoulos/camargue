# a script for converting TSPLIB 1-indexed .tour files to the 0-indexed .sol
# file format used by Camargue and Concorde.
# Input: instance_name.tour
#     A TSPLIB tour file.
# Output: instance_name.sol
#     A 0-indexed .sol file whose first line is the node count, with remaining
#     lines as the sequence of nodes to visit.
#

tfile="$1"
infix="$(echo "$1" | sed 's@\.tour@@')"
solfile="$(echo "$infix".sol)"

if [ -f "$solfile" ]
then
    echo "$solfile exists in working dir, please erase or rename to avoid overwriting"
    exit 1
fi

echo "Converting $tfile to $solfile"

dim="$(grep -m1 'DIMENSION' "$tfile" | tr -d -c [0-9])"

echo "$dim" > "$solfile"

grep '^[1-9]' "$tfile" | awk '{print ($1 -1)}' >> "$solfile"
