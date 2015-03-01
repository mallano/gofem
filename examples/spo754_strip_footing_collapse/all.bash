#!/bin/bash

gofem spo754 && GenVtu spo754
go run doplot.go
