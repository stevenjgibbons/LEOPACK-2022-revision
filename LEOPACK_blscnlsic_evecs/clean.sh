#!/bin/sh
for file in *OUTPUT* blscnlsic_evecs
do
    if test -r $file
    then
        rm $file
    fi
done
