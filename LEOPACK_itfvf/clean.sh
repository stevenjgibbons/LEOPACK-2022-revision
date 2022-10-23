#!/bin/sh
for file in *OUTPUT* itfvf
do
    if test -r $file
    then
        rm $file
    fi
done
