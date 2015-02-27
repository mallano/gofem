#!/bin/bash

gofem spo751 && GenVtu spo751
go run doplot.go
