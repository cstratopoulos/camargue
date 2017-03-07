# Test data generation script for Camargue
# This script will use Concorde to generate examples that are used in running
# Catch test cases. For each TSP instance in test_data/main_examples.txt, it
# will call Concorde to solve over the subtour polytope, saving the LP solution
# in test_data/subtour_lp and saving the tour in test_data/tours. For the few
# files in blossom_examples.txt, it will call Concorde to solve over the
# blossom polytope, saving the result in test_data/bossom_lp.
#
# To run this script, just make sure you have Concorde installed and built and
# pointed to in externals/ as described in the README.

CC_TSP=externals/concorde/TSP

for f in $(cat test_data/main_examples.txt); do
    have_f=$(ls test_data/tours | grep "$f" | wc -l)
    if [ "$have_f" -eq 0 ]; then
	problem=problems/"$f".tsp
	sub_fname="$f".sub.x
	sol_fname="$f".sol
	"$CC_TSP"/concorde -s99 -IxX test_data/subtour_lp/"$sub_fname" "$problem"
	mv "$sol_fname" test_data/tours/"$sol_fname"
	if [ "$?" -ne 0 ]; then
	    echo "Couldn't get subtour_lp/tour for " "$f"
	    exit 1
	fi
    fi
done

if [ ! -d test_data/blossom_lp ]; then
    mkdir test_data/blossom_lp
fi

for f in $(cat test_data/blossom_examples.txt); do
    have_f=$(ls test_data/blossom_lp | grep "$f" | wc -l)
    if [ "$have_f" -eq 0 ]; then
	problem=problems/"$f".tsp
	blos_fname="$f".sub.x
	"$CC_TSP"/concorde -s99 -ixX test_data/blossom_lp/"$blos_fname" "$problem"
	rm "$f".sol
	if [ "$?" -ne 0 ]; then
	    echo "Couldn't get blossom_lp for " "$f"
	    exit 1
	fi
    fi
done

exit 0
