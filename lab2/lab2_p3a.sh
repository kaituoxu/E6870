#!/bin/bash -e


if [[ -e lab2_fb ]] ; then
    binStr="lab2_fb"
elif [[ -e Lab2Fb.class ]] ; then
    binStr="java Lab2Fb"
else
    echo "Couldn't find program to execute."
    exit 1
fi


$binStr --audio_file p018k7.1.dat --graph_file p018k7.1.fsm  --iters 1 \
    --in_gmm p018k7.22.2.gmm --out_gmm /dev/null --chart_file p3a_chart.dat


