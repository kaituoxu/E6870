#!/bin/bash

fstcompile --isymbols=p1a.isyms --osymbols=p1b.osyms --keep_isymbols --keep_osymbols p1b.fst.txt > p1b.fst.binary
fstprint p1b.fst.binary
fstdraw p1b.fst.binary | dot -Tps > p1b.ps
