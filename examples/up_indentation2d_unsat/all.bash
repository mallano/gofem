#!/bin/bash

mpirun -np 4 gofem coarse-elast-d2-q9 && GenVtu coarse-elast-d2-q9
go run doplot.go

