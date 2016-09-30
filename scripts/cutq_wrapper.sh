for args; do
    for qcap in 1 5 25 50 100 150 250; do
	./cutq.sh "$qcap" "$args"
    done
done

for f in *.txt; do
    sort -k5n "$f" -o "$f"
done
