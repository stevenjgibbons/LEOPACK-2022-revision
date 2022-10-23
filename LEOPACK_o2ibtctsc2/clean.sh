#!/bin/sh
for file in *OUTPUT* o2ibtctsc2
do
    if test -r $file
    then
        rm $file
    fi
done
