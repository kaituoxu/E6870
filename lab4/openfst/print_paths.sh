#!/bin/bash

DIR=/Users/kaituoxu/Tools/openfst-tools
ISYM_TABLE=$1
FST=$2

$DIR/fstprintpaths $ISYM_TABLE $FST

# example
# $DIR/fstprintpaths p1a.isyms p1a.fst
