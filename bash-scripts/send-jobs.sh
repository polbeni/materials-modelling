#!/bin/bash

for i in {001..040}
do 
cd disp-$i
sbatch run.sh
cd ../
done
