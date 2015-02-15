# Copyright 2012 Dorival de Moraes Pedroso. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

all: sg114

debug:
	go test -gcflags "-N -l" -c -o /tmp/test_gofem_fem

beam01: debug
	gdb --args /tmp/test_gofem_fem -test.run=beam01

dom01: debug
	gdb --args /tmp/test_gofem_fem -test.run=dom01

bh16a: debug
	gdb --args /tmp/test_gofem_fem -test.run=bh16a

sg57: debug
	gdb --args /tmp/test_gofem_fem -test.run=sg57

sg114: debug
	gdb --args /tmp/test_gofem_fem -test.run=sg114

p01: debug
	gdb --args /tmp/test_gofem_fem -test.run=p01

bench:
	go test -run=XXX -bench=. -benchtime=5s
