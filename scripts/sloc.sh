#!/bin/sh
cat $* | # cat input files into pipe
sed '/^\ *\/\// d' | #delete // style comments
sed '/^\ *[ \/\\]\*/ d' | #delete /* style comments and \* and /* ' *'
#sed '/^#include/ d' | #delete include statements
wc -l
