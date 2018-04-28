#!/bin/bash -e


if [[ -e lab2_vit ]] ; then
    binStr="lab2_vit"
elif [[ -e Lab2Vit.class ]] ; then
    binStr="java Lab2Vit"
else
    echo "Couldn't find program to execute."
    exit 1
fi


$binStr --gmm p3c.gmm --audio_file p018k7t.10.dat \
    --graph_file p018k1.noloop.fsm --word_syms p018k2.syms \
    --dcd_file /dev/null


