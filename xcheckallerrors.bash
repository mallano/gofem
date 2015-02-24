#!/bin/bash

GOFEM="ana shp inp msolid mconduct mreten mporous fem out"

for p in $GOFEM; do
    echo
    echo
    echo "[1;32m>>> checking $p <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<[0m"
    errcheck github.com/cpmech/gofem/$p
done
