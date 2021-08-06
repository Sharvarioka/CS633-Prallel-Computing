#!/bin/bash
chmod +x script.sh
file=output.txt
if [ -f $file ] ; then
	rm $file
fi
make all    
for P in 1 2
do
  for ppn in 1 2 4
    do
	    ./script.sh $P
      val=`expr $P \* $ppn`
      mpiexec -np $val -f hostfile.txt -ppn $ppn ./src.x tdata.csv >> output.txt
    done
done
make clean

