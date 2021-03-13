#!/bin/bash

day=$1
fmt_day=$(printf "%02d" $day)

file=$2

if [ $# -eq 0 ]; then
    python sigjson_s2021a.py
    exit 0
fi

python sigjson_s2021a.py $day $2 > ../json/output.log && mv ../json/output.log ../json/202103"$fmt_day"/output.log

cat ../json/202103"$fmt_day"/output.log