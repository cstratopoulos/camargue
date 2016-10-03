for prob; do
    for qcap in 15 50 100; do
	./cutq.sh "$qcap" "$prob"
    done

    sort -nk3 "$prob"_ratios.txt -o "$prob"_ratios.txt
    sort -nk3 "$prob"_times.txt -o "$prob"_times.txt
    rm "$prob"_q*.txt
done
