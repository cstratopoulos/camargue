qsize=$1 #first argument is queue size to test
probname=$2 #second argument is problem name with no .tsp suffix
rawdata="$probname"_q"$qsize".txt
rawratios="$probname"_q"$qsize"_ratios.txt
rawtimes="$probname"_q"$qsize"_times.txt
>"$rawdata"
>"$rawratios"
>"$rawtimes"

for i in 1 2 3 4; do
    ../PSEP -S -c 4 -q "$qsize" ../problems/"$probname".tsp |
	grep 'blossom sep' >> "$rawdata"
done

sed 's/.*ratio. //g' "$rawdata" >> "$rawratios"
sed 's/.*sep. //g' "$rawdata" | sed 's/,.*//g' >> "$rawtimes"

sumtimes=$(paste -s -d + "$rawtimes")
sumratios=$(paste -s -d + "$rawratios")

avgtime=$(echo "($sumtimes)/4" | bc -l)
avgratio=$(echo "($sumratios)/4" | bc -l)

printf "%s %.3d %.6f\n" "$probname" "$qsize" "$avgtime" \
       >> "$probname"_times.txt

printf "%s %.3d %.6f\n" "$probname" "$qsize" "$avgratio" \
       >> "$probname"_ratios.txt

rm "$rawratios" "$rawtimes"
