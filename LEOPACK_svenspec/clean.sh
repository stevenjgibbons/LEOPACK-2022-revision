#!/bin/sh
for file in *OUTPUT* svenspec
do
    if test -r $file
    then
        rm $file
    fi
done
