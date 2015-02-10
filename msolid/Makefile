# Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

all: dp01

debug:
	go test -gcflags "-N -l" -c

dp01: debug
	gdb --args ./msolid.test -test.run=dp01
