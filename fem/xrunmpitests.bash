#!/bin/bash

go build -o /tmp/gofem/t_bh16_main t_bh16_main.go
mpirun -np 2 /tmp/gofem/t_bh16_main
