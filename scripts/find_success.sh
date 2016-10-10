PID=$$
quad=$1
probname=$2
prefix="$probname"_quad"$quad"_"$PID"
rawdata="$prefix"_raw.txt
>"$rawdata"


../PSEP -Sc2 -u"$quad"  ../problems/"$probname".tsp >> "$rawdata" 2>&1

opt=$(grep -i 'fathomed optimal' "$rawdata" | wc -l)

if [ "$opt" -ge 1 ]
then
    mv "$rawdata" "$prefix"_success.txt
    echo "Successful run!"
    echo -n $'\a'
    echo -n $'\a'
    echo -n $'\a'
else
    rm "$rawdata"
fi
    
