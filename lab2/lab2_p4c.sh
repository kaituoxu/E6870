#!/bin/bash -e


if [[ -e lab2_vit ]] ; then
    binStr="lab2_vit"
elif [[ -e Lab2Vit.class ]] ; then
    binStr="java Lab2Vit"
else
    echo "Couldn't find program to execute."
    exit 1
fi

echo "Decoding ..."
$binStr --gmm p4a.100.gmm --audio_file p018k7t.100.dat \
    --graph_file p018k7.100.5.noloop.fsm --word_syms p018k2.syms \
    --dcd_file p4c.100.dcd

echo "Computing WER ..."
p018h1.calc-wer.sh p4c.100.dcd p018k7t.100.trn temp


