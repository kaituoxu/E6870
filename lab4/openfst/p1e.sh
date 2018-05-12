#!/bin/bash

bash p1a_v2.sh > /dev/null # generate FSA p1a.fst
bash p1d.sh  > /dev/null  # generate FSA p1d.fst

fstconcat p1d.fst p1d.fst | fstconcat p1d.fst - p1d_x3.fst
fstdraw --acceptor p1d_x3.fst | dot -Tps > p1d_x3.ps

fstdeterminize p1d_x3.fst | fstdifference - p1a.fst p1e.fst
fstdraw --acceptor p1e.fst | dot -Tps > p1e.ps
