#!/bin/bash -e


if [[ -e lab2_fb ]] ; then
    binStr="lab2_fb"
elif [[ -e Lab2Fb.class ]] ; then
    binStr="java Lab2Fb"
else
    echo "Couldn't find program to execute."
    exit 1
fi


$binStr --audio_file p018k7.22.dat --graph_file p018k7.22.fsm --iters 20 \
    --in_gmm p018k1.gmm.dat --out_gmm p3c.gmm


