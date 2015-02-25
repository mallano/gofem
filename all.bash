#!/bin/bash

GOFEM="ana shp inp msolid mconduct mreten mporous fem out"

HERE=`pwd`
for p in $GOFEM; do
    echo
    echo
    echo "[1;32m>>> compiling $p <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<[0m"
    cd $p
    touch *.go
    go test
    go install
    cd $HERE
done

echo
echo
echo "[1;32m>>> compiling binaries <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<[0m"
make
