#!/bin/sh
for file in *OUTPUT* cicibcdts2
do
    if test -r $file
    then
        rm $file
    fi
done
