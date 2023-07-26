#!/bin/bash

for i in {001..030}
do 
echo simulation $i
tail -n 1 disp-$i/OUTCAR
done
