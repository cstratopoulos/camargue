for prob; do
    for maxround in 15 25 50 100; do
	./cutsper.sh "$maxround" "$prob"
    done

    for suff in _pivratios.txt _soltimes.txt _failratios.txt; do
	sort -nk3 "$prob""$suff" -o "$prob""$suff"
    done
done
