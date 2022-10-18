#!/bin/sh
for file in *OUTPUT* cicm2ocdisplay
do
    if test -r $file
    then
        rm $file
    fi
done
