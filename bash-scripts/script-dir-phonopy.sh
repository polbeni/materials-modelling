#!/bin/bash

for i in {001..100}
do 
cp generated_POSCAR/POSCAR-$i simulation-$i/POSCAR
cp save/* simulation-$i
done
