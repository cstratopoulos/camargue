for prob; do
    for maxround in 1 2 4 8; do
	./cutsper.sh "$maxround" "$prob"
    done

    for suff in _pivratios.txt _pivtimes.txt _soltimes.txt; do
	sort -nk3 "$prob""$suff" -o "$prob""$suff"
    done
done
