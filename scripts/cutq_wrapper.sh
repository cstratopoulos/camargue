for args; do
    for qcap in 1 15 50 100; do
	./cutq.sh "$qcap" "$args"
    done
done