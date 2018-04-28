#!/bin/bash -e



for d in 100 300 1000 ; do

echo "======================================================================"
echo "Training: $d sentences."
echo "======================================================================"

if [[ -e lab2_fb ]] ; then
    binStr="lab2_fb"
elif [[ -e Lab2Fb.class ]] ; then
    binStr="java Lab2Fb"
else
    echo "Couldn't find program to execute."
    exit 1
fi

$binStr --audio_file p018k7.$d.dat --graph_file p018k7.$d.fsm --iters 5 \
    --in_gmm p018k1.gmm.dat --out_gmm p4a.$d.gmm

if [[ -e lab2_vit ]] ; then
    binStr="lab2_vit"
elif [[ -e Lab2Vit.class ]] ; then
    binStr="java Lab2Vit"
else
    echo "Couldn't find program to execute."
    exit 1
fi

echo "Decoding ..."
$binStr --gmm p4a.$d.gmm --audio_file p018k7t.100.dat \
    --graph_file p018k1.noloop.fsm --word_syms p018k2.syms \
    --dcd_file p4a.$d.dcd

echo "Computing WER ..."
p018h1.calc-wer.sh p4a.$d.dcd p018k7t.100.trn temp

done


