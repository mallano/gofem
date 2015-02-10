# Copyright 2012 Dorival de Moraes Pedroso. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

all: dom01

debug:
	go test -gcflags "-N -l" -c

beam01: debug
	gdb --args ./fem.test -test.run=beam01

dom01: debug
	gdb --args ./fem.test -test.run=dom01

bh16a: debug
	gdb --args ./fem.test -test.run=bh16a

sg57: debug
	gdb --args ./fem.test -test.run=sg57

bench:
	go test -run=XXX -bench=. -benchtime=5s
