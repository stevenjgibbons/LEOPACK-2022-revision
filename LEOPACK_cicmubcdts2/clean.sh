#!/bin/sh
for file in *OUTPUT* cicmubcdts2
do
    if test -r $file
    then
        rm $file
    fi
done
