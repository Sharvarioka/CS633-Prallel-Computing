#!/bin/bash

rm hostfile.txt > /dev/null 2>&1
touch hostfile.txt
#hosts=()
count=0
for i in $(seq 1 1 39)
do
    ping -c2 csews$i.cse.iitk.ac.in > /dev/null
    if [ $? -eq 0 ]
    then
       # hosts+=("csews$i")
        echo "csews$i" >> hostfile.txt
        count=$((count+1))
        if [ $count == $1 ]
        then
            break
        fi
    fi
done
