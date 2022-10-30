#!/bin/sh
for file in *OUTPUT* djiepgrf
do
    if test -r $file
    then
        rm $file
    fi
done
