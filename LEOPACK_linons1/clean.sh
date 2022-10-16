#!/bin/sh
for file in *OUTPUT* linons1
do
    if test -r $file
    then
        rm $file
    fi
done
