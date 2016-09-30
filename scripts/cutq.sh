qsize=$1 #first argument is queue size to test
probname=$2 #second argument is problem name with no .tsp suffix
rawtimes="$probname"_q"$qsize".txt
>"$rawtimes"

for i in 1 2 3 4 5; do
    ../PSEP -S -s 1000 -c 1 -q "$qsize" ../problems/"$probname".tsp |
	grep 'blossom sep' >> "$rawtimes"
done

cat "$rawtimes" | tr -d -c '[.0-9\n]' > "$rawtimes"tmp.txt
mv "$rawtimes"tmp.txt "$rawtimes"

sumvar=$(paste -s -d + "$rawtimes")
avg=$(echo "($sumvar)/5" | bc -l)
printf "%s %.3d avg %.6f\n" "$probname" "$qsize" "$avg" \
       >> "$probname"_q_reports.txt

rm "$rawtimes"
