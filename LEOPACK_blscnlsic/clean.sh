#!/bin/sh
for file in *OUTPUT* blscnlsic
do
    if test -r $file
    then
        rm $file
    fi
done
