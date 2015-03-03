#!/bin/bash

gofem coarse-elast-d2-q9 && GenVtu coarse-elast-d2-q9
go run doplot.go

