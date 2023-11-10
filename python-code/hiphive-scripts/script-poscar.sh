#!/bin/bash

for i in {0001..0010}
do 
mkdir disp-$i
mv POSCAR-$i disp-$i/POSCAR
cp save/* disp-$i
done
