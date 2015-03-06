#!/bin/bash

mpirun -np 4 gofem o2elast.sim
go run doplot.go
