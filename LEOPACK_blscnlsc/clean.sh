#!/bin/sh
for file in *OUTPUT* blscnlsc
do
    if test -r $file
    then
        rm $file
    fi
done
