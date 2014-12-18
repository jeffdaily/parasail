#!/bin/sh

for f in nw_scr_*
do
    if ! diff nw_scr.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

for f in sg_scr_*
do
    if ! diff sg_scr.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

for f in sw_scr_*
do
    if ! diff sw_scr.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

