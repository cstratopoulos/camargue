PID=$$
probname=$1
rawdata="$probname"_"$PID"
>"$rawdata".txt

../PSEP -Sc2 -u10  ../problems/"$1".tsp >> "$rawdata".txt 2>&1

opt=$(grep -i 'optimal' "$rawdata".txt | wc -l)

if [ "$opt" -ge 1 ]
then
    mv "$rawdata".txt "$rawdata"_success.txt
    echo "Successful run!"
    echo -n $'\a'
    echo -n $'\a'
    echo -n $'\a'
else
    rm "$rawdata".txt
fi
    
