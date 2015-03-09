#!/bin/bash

FILES="onepulse-qua9co linear-qua9co wet-linear-qua9co"

for f in $FILES; do
    mpirun -np 4 gofem $f
#    GenVtu $f
    go run doplot.go $f
done
