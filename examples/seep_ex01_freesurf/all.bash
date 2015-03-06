#!/bin/bash

mpirun -np 4 gofem coarse && GenVtu coarse
