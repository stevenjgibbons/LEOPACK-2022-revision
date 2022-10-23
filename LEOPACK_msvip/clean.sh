#!/bin/sh
for file in *OUTPUT* msvip
do
    if test -r $file
    then
        rm $file
    fi
done
