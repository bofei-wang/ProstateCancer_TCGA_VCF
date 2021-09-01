#!/bin/bash
#set -x

file=$1
echo "file name $1"
path1=$2
path2=$3
path3=$4
path4=$5
cat $file | while read line
#for line in $(cat $file|sed -n '50,100p')
do 
    echo "instance:$line"
    seq=$(echo "$path1$line" | tr -d '\r')
    var=$(echo "$path2$line" | tr -d '\r')
    out=$(echo "$path3$line" | tr -d '\r')
    err=$(echo "$path4$line" | tr -d '\r')
    #stdout=$(echo "$path5$line" | tr -d '\r')
    seqf=$seq".fasta"
    varf=$var".var"
    outf=$out".txt"
    errf=$err".err"
    #stdoutf=$stdout".out"
    bsub  -q bio_priority -W 02:00 -F 3000000 -n 5 -o $outf -e $errf /home/oncominer/bbrcexome/bin/provean.sh -q $seqf -v $varf --num_threads 5
    #bsub -q bio_priority -n 5 -o $outf -e $errf /home/oncominer/bbrcexome/bin/provean.sh -q $seqf -v $varf --num_threads 5 
done
 
