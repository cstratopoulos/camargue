if [ "$#" -ne 1 ]; then
    echo "Bad argument count: $#"
    exit 1
fi

tspname="$1"
opt="$(grep "$tspname" test_data/opts.txt | sed 's@.*:@@g' | tr -dc '[0-9]')"

foundopt="$(echo "$opt" | wc -c)"
if [ "$foundopt" -le 1 ]; then
    echo "Opt length $opt looks invalid"
    exit 1
fi

init_pre="$tspname"-pre
full_pre="$tspname"-full
opt_name="$tspname"."$opt".sol

dir_name=chained-CMR-runs
if [ ! -d "$dir_name" ]; then
    mkdir "$dir_name"
fi

main_pre="$tspname"-main
main_run="$main_pre".run

for trial in $(seq 1 3); do
    if [ -e "$opt_name" ]; then
	echo "Opt tourfile $opt_name present before trial $trial, exiting" >> \
	     "$dir_name"/"$main_run"
	exit 0
    fi
    this_pre="$init_pre$trial"
    this_run="$this_pre".run

    timeout 5m ./camargue -b1 -l"$opt" -GPS problems/"$tspname".tsp \
	    >> "$this_run"
    mv "$this_run" "$dir_name"/"$this_run"
done

if [ -e "$opt_name" ]; then
    echo "Opt tourfile exists after pre-trials, exiting" >> \
	 "$dir_name"/"$main_run"
    exit 0
fi

best_sol="$(ls | grep "$tspname.*sol" | head -n1)"

timeout 1h ./camargue -e1 -b1 -l"$opt" -t"$best_sol" \
	-GS problems/"$tspname".tsp >> "$main_run"
mv "$main_run" "$dir_name"/"$main_run"
