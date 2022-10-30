#!/bin/sh
for file in *OUTPUT* kriepgrf
do
    if test -r $file
    then
        rm $file
    fi
done
