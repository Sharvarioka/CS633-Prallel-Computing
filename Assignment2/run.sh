#!/bin/bash
make all
file1=plot_Bcast.csv
file2=plot_Gather.csv
file3=plot_Reduce.csv
file4=plot_Alltoallv.csv
    if [ -f $file1 ] ; then
       rm $file1
    fi
    if [ -f $file2 ] ; then
       rm $file2
    fi
    if [ -f $file3 ] ; then
       rm $file3
    fi
    if [ -f $file4 ] ; then
       rm $file4
    fi
echo "Type","Time"> plot_Bcast.csv 
echo "Type","Time"> plot_Gather.csv
echo "Type","Time"> plot_Reduce.csv
echo "Type","Time"> plot_Alltoallv.csv
for execution in 1 2 3 4 5 6 7 8 9 10
do      
    for P in 4 16
    do
        for ppn in 1 8
        do
         var=`expr $P / 2`
          python script.py 2 $var $ppn
          for D in 16 256 2048
          do
                  mpiexec -np $P*$ppn -f hostfile ./src.x $D $P $ppn
          done
        done
    done
done
python plot.py
make clean

