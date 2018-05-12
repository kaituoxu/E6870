#!/bin/bash
# OpenFst example

# compile
fstcompile --isymbols=isyms.txt --osymbols=osyms.txt --keep_isymbols --keep_osymbols text.fst binary.fst

# draw
fstdraw --isymbols=isyms.txt --osymbols=osyms.txt binary.fst binary.dot
dot -Tps binary.dot > binary.ps
