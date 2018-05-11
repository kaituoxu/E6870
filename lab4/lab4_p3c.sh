#!/bin/bash -e


dbgStr=""
while [[ ( $# -gt 0 ) && ( "$1" == -?* ) ]] ; do
    strFlag="$1" ; shift
    case "$strFlag" in
        -dbg) dbgStr="dbg.sh" ;;
        *) echo "Unrecognized flag: $strFlag" 1>&2 ; exit 1 ;;
    esac
done

if [[ -e lab4_vit ]] ; then
    binStr="lab4_vit"
elif [[ -e Lab4Vit.class ]] ; then
    binStr="java Lab4Vit"
else
    echo "Couldn't find program to execute."
    exit 1
fi


$dbgStr $binStr --gmm isodigit.gmm --audio_file isodigit.1.dat \
    --graph_file isodigit.eps.fsm --word_syms isodigit.syms \
    --dcd_file /dev/null --chart_file p3c.chart \
    --ac_wgt 0.1 --log10_beam 999 --rank_beam 0


