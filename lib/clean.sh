#!/bin/sh
for file in *.a
do
    if test -r $file
    then
        rm $file
    fi
done
