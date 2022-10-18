#!/bin/sh
for file in *OUTPUT* cicmibcdts2
do
    if test -r $file
    then
        rm $file
    fi
done
