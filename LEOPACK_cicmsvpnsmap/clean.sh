#!/bin/sh
for file in *OUTPUT* cicmsvpnsmap
do
    if test -r $file
    then
        rm $file
    fi
done
