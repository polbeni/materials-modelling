#!/bin/bash

printf -v folders 'dist-%03d ' {1..3}
mkdir $folders