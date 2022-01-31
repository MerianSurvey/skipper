#!/bin/bash

day=$1
month=$2
year=$3
fmt_day=$(printf "%02d" $day)
fmt_month=$(printf "%02d" $month)

file=$4
logfilename="../json/"$year$fmt_month$fmt_day"/output.log"

if [ $# -eq 0 ]; then
    python sigjson_s2022a.py
    exit 0
fi

python sigjson_s2022a.py $@ > ../json/output.log && mv ../json/output.log $logfilename && cat $logfilename