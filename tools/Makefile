# Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

all: GenVtu ConvertGofemMat MatTable
.PHONY: GenVtu ConvertGofemMat

ConvertGofemMat: ConvertGofemMat.go
	go build -o /tmp/gofem/ConvertGofemMat ConvertGofemMat.go && mv /tmp/gofem/ConvertGofemMat $(GOPATH)/bin/

GenVtu: GenVtu.go
	go build -o /tmp/gofem/GenVtu GenVtu.go && mv /tmp/gofem/GenVtu $(GOPATH)/bin/

MatTable: MatTable.go
	go build -o /tmp/gofem/MatTable MatTable.go && mv /tmp/gofem/MatTable $(GOPATH)/bin/
