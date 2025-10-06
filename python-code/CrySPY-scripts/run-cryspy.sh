#!/bin/bash

# Number of times to run the command
TIMES=10

# Delay in seconds between iterations
DELAY=30

for i in $(seq 1 $TIMES); do
    cryspy
    sleep $DELAY
done