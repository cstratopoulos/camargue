for prob; do
    for maxround in 1 2 4 8 12 16 25 50; do
	./cutsper.sh "$maxround" "$prob"
    done

    for suff in _pivratios.txt _soltimes.txt _failratios.txt; do
	sort -nk3 "$prob""$suff" -o "$prob""$suff"
    done
done
