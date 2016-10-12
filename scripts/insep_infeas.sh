#this script will run ./PSEP on its first argument, saving the result
#to file if it contains an inseparable integral subtour or a best tour
#made infeasible upon pivot back
#
#should be safe for running in batches on the same problem by use of process
#id to protect filename

PID=$$
probname=$1
rawdata="$probname"_"$PID"
>"$rawdata".txt

sleep 2s
../PSEP -Sc2 -u5 ../problems/"$probname".tsp >> "$rawdata".txt 2>&1

insep=$(grep -il 'inseparable' "$rawdata".txt | wc -l)
infeas=$(grep -il 'infeasible' "$rawdata".txt | wc -l)

if [ \( "$insep" -eq 0 \) -a \( "$infeas" -eq 0 \) ]
then
    rm "$rawdata".txt
else
    mv "$rawdata".txt "$rawdata"_insep"$insep"_infeas"$infeas".txt
fi

