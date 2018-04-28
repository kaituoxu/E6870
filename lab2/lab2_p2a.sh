#!/bin/bash -e


if [[ -e lab2_train ]] ; then
    binStr="lab2_train"
elif [[ -e Lab2Train.class ]] ; then
    binStr="java Lab2Train"
else
    echo "Couldn't find program to execute."
    exit 1
fi


$binStr --audio_file p018k7.22.dat --align_file p018k7.22.20.align \
    --iters 1 --in_gmm p018k1.gmm.dat --out_gmm p2a.gmm


