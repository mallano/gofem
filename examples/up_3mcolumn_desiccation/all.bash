#!/bin/bash

mpirun -np 4 gofem onepulse-qua9co
GenVtu onepulse-qua9co
go run doplot.go
