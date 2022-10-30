#!/bin/sh
for file in *OUTPUT* krcmrnif
do
    if test -r $file
    then
        rm $file
    fi
done
