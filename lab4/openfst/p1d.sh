#!/bin/bash

fstcompile --isymbols=p1a.isyms --keep_isymbols --acceptor p1d.fsm p1d.fst
fstprint --acceptor p1d.fst
fstdraw --acceptor p1d.fst | dot -Tps > p1d.fsm.ps
