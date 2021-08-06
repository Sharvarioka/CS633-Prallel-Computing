#!/bin/bash
make all

for P in 16 36 49 64

do
    file=plot${P}.csv
    if [ -f $file ] ; then
       rm $file
    fi
    echo "Type","N","Time"> plot${P}.csv
    ~/UGP/allocator/src/allocator.out ${P} 8 > log
    
    for N in 256 1024 4096 16384 65536 262144 1048576
    do
        for t in 1 2 3 4 5
        do
            mpiexec -np ${P} -f hosts ./src.x ${N} 50
                                 
        done
    done
    python3 plot.py $P

done
make clean
