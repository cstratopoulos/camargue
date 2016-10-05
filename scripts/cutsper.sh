# This script is for experimenting with the number of cuts added per round
# That number can be interpreted as cuts of each type per round, i.e.,
# a value of 2 means we would add at most 2 segments and 2 blossoms (if we
# can find that many), and if we can't find any we will try to add at most 2
# safe GMI cuts
#
# Output data of interest:
# (Time related)
#    -Total solution time
#    -Total time pivoting
#    -Ratio of pivot time to solution time
# (Cut related) (later?)
#    -Total number of cuts added
#    -Total number of rounds of separation
#    -Ratio of first to second
#

rsize=$1 #first argument is number of cuts of each type to be added per round
probname=$2 #problem name with no tsp suffix
prefix="$probname"_c"$rsize"
rawdata="$prefix".txt
rawpiv="$prefix"_pivot.txt
rawsol="$prefix"_sol.txt

>"$rawdata"
>"$rawpiv"
>"$rawsol"

for i in 1 2 3 4 5; do
    ../PSEP -S -c "$rsize" ../problems/"$probname".tsp 2>&1 |
	grep 'pivoting\|solve\|failed' >> "$rawdata"
done

cat "$rawdata" | grep 'Total.*solve' >> "$rawsol"
cat "$rawdata" | grep 'pivot' >> "$rawpiv"

pivtimes="$prefix"_piv_times.txt
pivratios="$prefix"_piv_ratios.txt
soltimes="$prefix"_sol_times.txt

>"$pivtimes"
>"$pivratios"
>"$soltimes"

cat "$rawsol" |
    sed 's/.*solve. //g' >> "$soltimes"
cat "$rawpiv" |
    sed 's/.*pivoting. //g' |
    sed 's/,.*//g' >> "$pivtimes"
cat "$rawpiv" |
    sed 's/.*ratio. //g' >> "$pivratios"

rm "$rawsol" "$rawpiv"

sumpivtimes=$(paste -s -d + "$pivtimes")
sumpivratios=$(paste -s -d + "$pivratios")
sumsoltimes=$(paste -s -d + "$soltimes")
failcount=$(cat "$rawdata" | grep 'failed' | wc -l)

rm "$pivtimes" "$pivratios" "$soltimes" "$rawdata"

avgpivtime=$(echo "($sumpivtimes)/5" | bc -l)
avgpivratio=$(echo "($sumpivratios)/5" | bc -l)
avgsoltime=$(echo "($sumsoltimes)/5" | bc -l)
failratio=$(echo "($failcount)/5" | bc -l)

printf "%s %.2d %.6f\n" "$probname" "$rsize" "$avgpivtime" \
       >> "$probname"_pivtimes.txt

printf "%s %.2d %.6f\n" "$probname" "$rsize" "$avgpivratio" \
       >> "$probname"_pivratios.txt

printf "%s %.2d %.6f\n" "$probname" "$rsize" "$avgsoltime" \
       >> "$probname"_soltimes.txt

printf "%s %.2d %.1f\n" "$probname" "$rsize" "$failratio" \
       >> "$probname"_failratios.txt
