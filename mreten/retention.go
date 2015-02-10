package mdl

import (
    "github.com/cpmech/gosl/utl"
)

type Retention interface {
    Init ()
}

var rm_factory = map[string]func() Retention {}

func NewRetention(name string) Retention {
    allocator, ok := rm_factory[name]
    if !ok {
        utl.Panic(_retention_err01, name)
    }
    return allocator()
}

// error messages
var (
    _retention_err01 = "retention.go: NewRetention: cannot find retention model named %s"
)
