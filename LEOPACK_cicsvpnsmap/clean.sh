#!/bin/sh
for file in *OUTPUT* cicsvpnsmap
do
    if test -r $file
    then
        rm $file
    fi
done
