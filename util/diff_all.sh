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

for f in nw_stats_scr_*
do
    if ! diff nw_stats_scr.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

for f in sg_stats_scr_*
do
    if ! diff sg_stats_scr.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

for f in sw_stats_scr_*
do
    if ! diff sw_stats_scr.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

for f in nw_stats_mch_*
do
    if ! diff nw_stats_mch.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

for f in sg_stats_mch_*
do
    if ! diff sg_stats_mch.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

for f in sw_stats_mch_*
do
    if ! diff sw_stats_mch.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

for f in nw_stats_len_*
do
    if ! diff nw_stats_len.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

for f in sg_stats_len_*
do
    if ! diff sg_stats_len.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

for f in sw_stats_len_*
do
    if ! diff sw_stats_len.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

