#this script will run ./PSEP on its first argument, saving the result
#to file if it contains an augmented tour even and/or an inseparable
#integral subtour event
#
#should be safe for running in batches on the same problem by use of process
#id to protect filename

PID=$$
probname=$1
rawdata="$probname"_"$PID"
>"$rawdata".txt

sleep 2s
../PSEP -Sc2 -u10 ../problems/"$probname".tsp >> "$rawdata".txt 2>&1

insep=$(grep -il 'inseparable' "$rawdata".txt | wc -l)
aug=$(grep -i 'augmented' "$rawdata".txt | wc -l)

if [ \( "$insep" -eq 0 \) -a \( "$aug" -eq 0 \) ]
then
    rm "$rawdata".txt
else
    mv "$rawdata".txt "$rawdata"_i"$insep"_a"$aug".txt
fi

