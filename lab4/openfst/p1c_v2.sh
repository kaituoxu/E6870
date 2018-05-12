#!/bin/bash

fstcompose p1a.fst p1b.fst.binary p1c.fst.binary
fstprint p1c.fst.binary
fstdraw p1c.fst.binary | dot -Tps > p1c.fsm.ps
