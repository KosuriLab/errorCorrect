#!/bin/bash
# parDiff - create diff strings in parallel with bash trickery!
# Nathan Lubock
#
# Outline:
#   create diff string and delete all duplicate reads cols
#   spawn a child proc for each cpu
#   wait for all to finish before submitting next batch
#   && tests success of previous command
# Input:
#   $1 - list of sam files to run
#   $2 - /path/to/reference.fasta
# Output:
#   Diff string with *.diff.clean.sam according to file names in current dir
echo "Generating diff strings and cleaning files"
count=0
numproc=`nproc`
while read f; do
    samtools calmd -e $f $2 2> `basename ${f%.sam}.diff.err` |
    sed 's/\tXS\:i\:-\?[0-9]*//' > `basename ${f%.sam}.diff.clean.sam` &
    let count+=1
    [[ $((count%numproc)) -eq 0 ]] && wait
done < $1
wait # make sure that you wait until this whole loop is finished
