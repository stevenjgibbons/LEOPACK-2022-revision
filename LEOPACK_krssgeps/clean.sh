#!/bin/sh
for file in *OUTPUT* krssgeps
do
    if test -r $file
    then
        rm $file
    fi
done
