# Copyright 2012 Dorival de Moraes Pedroso. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

all: mdl01

debug:
	go test -gcflags "-N -l" -c -o /tmp/test_gofem_mporous

mdl01: debug
	gdb --args /tmp/test_gofem_mporous -test.run=mdl01
