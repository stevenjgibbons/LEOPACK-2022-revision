#!/bin/sh
for file in *OUTPUT* cicubcdts2
do
    if test -r $file
    then
        rm $file
    fi
done
