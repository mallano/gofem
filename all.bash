#!/bin/bash

GOFEM="shp inp msolid mreten mporous fem"

HERE=`pwd`
for p in $GOFEM; do
    echo
    echo
    echo "[1;32m>>> compiling $p <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<[0m"
    cd $p
    go test
    go install
    cd $HERE
done

echo
echo
echo "[1;32m>>> compiling binaries <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<[0m"
make
