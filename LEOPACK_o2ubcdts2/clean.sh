#!/bin/sh
for file in *OUTPUT* o2ubcdts2
do
    if test -r $file
    then
        rm $file
    fi
done
