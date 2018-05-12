#!/bin/bash

fstcompile --isymbols=p1a.isyms --osymbols=p1a.osyms --keep_isymbols --keep_osymbols p1a.fst.txt > p1a.fst.binary
fstprint p1a.fst.binary
fstdraw p1a.fst.binary | dot -Tps > p1a.ps
# fstcompile --isymbols=p1a.isyms --osymbols=p1a.osyms p1a.fsm > p1a.fst.binary
# fstprint --isymbols=p1a.isyms --osymbols=p1a.osyms p1a.fst.binary
# fstdraw --isymbols=p1a.isyms --osymbols=p1a.osyms p1a.fst.binary | dot -Tps > p1a.ps
