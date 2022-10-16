#!/bin/sh
for file in *OUTPUT* linons2
do
    if test -r $file
    then
        rm $file
    fi
done
