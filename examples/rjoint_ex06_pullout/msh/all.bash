#!/bin/bash

FILES="o2"

for f in $FILES; do
    echo
    echo "[1;33m>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> $f <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<[0m"
    glabgenmsh $f.inp
    rm $f.blocks "$f"_lines.py
done
