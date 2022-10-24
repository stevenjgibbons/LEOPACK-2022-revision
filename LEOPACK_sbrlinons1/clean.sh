#!/bin/sh
for file in *OUTPUT* sbrlinons1
do
    if test -r $file
    then
        rm $file
    fi
done
