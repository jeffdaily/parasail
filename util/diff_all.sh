#!/bin/sh

message() {
    if grep -q "saturated=1" $1
    then
        echo "$1 differs, but saturated"
    else
        echo "$1 differs"
    fi
}

for a in nw sg_qb sg_qe sg_qx sg_db sg_de sg_dx sg_qb_de sg_qe_db sg_qb_db sg_qe_de sg sw
do
for f in ${a}_scr_*
do
    if ! diff ${a}_scr_orig_NA_32_32.txt $f > /dev/null
    then
        message $f
    fi
done
done

for a in nw sg_qb sg_qe sg_qx sg_db sg_de sg_dx sg_qb_de sg_qe_db sg_qb_db sg_qe_de sg sw
do
for t in scr mch len
do
for f in ${a}_stats_${t}_*
do
    if ! diff ${a}_stats_${t}_orig_NA_32_32.txt $f > /dev/null
    then
        message $f
    fi
done
done
done

