#!/bin/sh
for file in *OUTPUT*
do
    if test -r $file
    then
        rm $file
    fi
done
