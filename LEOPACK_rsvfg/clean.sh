#!/bin/sh
for file in *OUTPUT* rsvfg
do
    if test -r $file
    then
        rm $file
    fi
done
