#!/bin/sh
for file in *OUTPUT* o2ubtctsc2
do
    if test -r $file
    then
        rm $file
    fi
done
