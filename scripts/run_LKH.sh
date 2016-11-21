prob=$1
opt=$(grep "$prob" opts.txt | cut -f2 -d: | tr -d -C '[0-9\n]')
fname="$prob".par

echo PROBLEM_FILE = problems/"$prob".tsp > "$fname"
echo OPTIMUM = "$opt" >> "$fname"
echo MAX_TRIALS = 10 >> "$fname"
echo ASCENT_CANDIDATES = 25 >> "$fname"
echo TOUR_FILE = "$prob".sol >> "$fname"

./LKH "$fname"
rm "$fname"
