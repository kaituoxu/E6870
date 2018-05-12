#!/bin/bash

fstcompile --isymbols=p1a.isyms --keep_isymbols --acceptor p1a.fsm p1a.fst
fstprint --acceptor p1a.fst
fstdraw --acceptor p1a.fst | dot -Tps > p1a.fsm.ps
