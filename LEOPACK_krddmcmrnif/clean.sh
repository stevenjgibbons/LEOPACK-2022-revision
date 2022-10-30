#!/bin/sh
for file in *OUTPUT* krddmcmrnif
do
    if test -r $file
    then
        rm $file
    fi
done
